import os
import subprocess as sp
import tempfile
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from typing import Dict, List, Tuple, Union

from dxcint.Context import Context, ContextEmpty
from dxcint.RegisteredTest import RegisteredTest


# Token expected by the in-process HTTP server. Mirrors the value used by the
# scalatest integration spec; keep these in sync if it ever needs to change.
_EXPECTED_TOKEN = "secret-token-abc"

# Env variable that dxCompiler reads to pick up per-domain Bearer tokens for
# WDL imports. Mirrored from dx.core.languages.wdl.WdlImportHttpAuth.TokensEnvVar.
_BEARER_TOKENS_ENV_VAR = "DXCOMPILER_WDL_IMPORT_BEARER_TOKENS"

# Filename of the protected WDL served by the local HTTP server. Lives next to
# the main workflow template in the test fixture directory.
_PROTECTED_WDL_FILENAME = "imported.wdl"

# Placeholder in the main WDL fixture that is replaced with the live server
# URL at test time.
_URL_PLACEHOLDER = "{URL}"


class BearerAuthImport(RegisteredTest):
    """Local-only test that verifies dxCompiler handles
    `DXCOMPILER_WDL_IMPORT_BEARER_TOKENS` when following http(s) WDL imports.

    The test brings up an HTTP server that serves a single WDL document behind
    a static Bearer token, renders the main workflow with the server URL, and
    invokes `java -jar dxCompiler.jar compile ... -compileMode IR` three times:
    with no token, a wrong token, and the correct token. All three must fail
    with HTTP 401 because bearer credentials are never attached to plain HTTP
    imports.

    It does not interact with the platform: there is no upload, no DXAnalysis,
    and no messenger. `get_test_result` is overridden so the normal
    compile/run/wait pipeline is bypassed entirely.
    """

    def __init__(
        self,
        src_file: str,
        category: str,
        test_name: str,
        context: Union[Context, ContextEmpty],
    ):
        super().__init__(src_file, category, test_name, context)
        self._fixture_dir = os.path.dirname(src_file)
        with open(src_file, "r") as f:
            self._main_template = f.read()
        protected_path = os.path.join(self._fixture_dir, _PROTECTED_WDL_FILENAME)
        with open(protected_path, "rb") as f:
            self._protected_bytes = f.read()

    # --- the normal RegisteredTest lifecycle is bypassed for this test ----

    def get_test_result(self) -> bool:
        result = self._validate()
        if result["passed"]:
            self._context.logger.info(
                f"Test {self._test_name} successfully PASSED. {result['message']}"
            )
            return True
        self._context.logger.error(
            f"Test {self._test_name} FAILED with message: {result['message']}"
        )
        return False

    def _validate(self) -> Dict:
        server, host, port = self._start_protected_server(
            self._protected_bytes, _EXPECTED_TOKEN
        )
        try:
            url = f"http://{host}:{port}/{_PROTECTED_WDL_FILENAME}"
            scenarios: List[Tuple[str, Union[str, None], bool, List[str]]] = [
                (
                    "no token configured",
                    None,
                    False,
                    [
                        "HTTP 401 Unauthorized",
                        "ensure the credentials are provided",
                    ],
                ),
                (
                    "wrong token configured",
                    f"{host}:wrong-token",
                    False,
                    [
                        "HTTP 401 Unauthorized",
                        "ensure the credentials are provided",
                    ],
                ),
                (
                    "correct token configured over HTTP",
                    f"{host}:{_EXPECTED_TOKEN}",
                    False,
                    [
                        "HTTP 401 Unauthorized",
                        "ensure the credentials are provided",
                    ],
                ),
            ]
            with tempfile.TemporaryDirectory(prefix="dxcint-bearer-auth-") as workdir:
                main_path = os.path.join(workdir, "main.wdl")
                with open(main_path, "w") as f:
                    f.write(self._main_template.replace(_URL_PLACEHOLDER, url))
                failures: List[str] = []
                for label, env_value, expect_success, must_include in scenarios:
                    rc, combined = self._run_compile(main_path, env_value)
                    ok, msg = self._check_outcome(
                        rc, combined, expect_success, must_include
                    )
                    if ok:
                        self._context.logger.info(
                            f"BearerAuthImport: PASS [{label}] ({msg})"
                        )
                    else:
                        self._context.logger.error(
                            f"BearerAuthImport: FAIL [{label}] ({msg}); "
                            f"output: {combined.strip()[:1000]}"
                        )
                        failures.append(label)
                if failures:
                    return {
                        "passed": False,
                        "message": (
                            f"{len(failures)}/{len(scenarios)} bearer-auth "
                            f"scenarios failed: {', '.join(failures)}"
                        ),
                    }
                return {
                    "passed": True,
                    "message": (f"all {len(scenarios)} bearer-auth scenarios passed"),
                }
        finally:
            server.shutdown()
            server.server_close()

    # --- helpers ----------------------------------------------------------

    @staticmethod
    def _start_protected_server(
        body_bytes: bytes, expected_token: str
    ) -> Tuple[HTTPServer, str, int]:
        protected_path = "/" + _PROTECTED_WDL_FILENAME

        class Handler(BaseHTTPRequestHandler):
            def _auth_status(self) -> int:
                # 401 when no credentials are supplied, 403 when credentials are
                # supplied but do not match, 200 when they match.
                auth = self.headers.get("Authorization")
                if auth is None:
                    return 401
                if auth != f"Bearer {expected_token}":
                    return 403
                return 200

            def _serve(self, send_body: bool) -> None:
                if self.path != protected_path:
                    self.send_response(404)
                    self.end_headers()
                    return
                status = self._auth_status()
                if status != 200:
                    self.send_response(status)
                    self.end_headers()
                    return
                self.send_response(200)
                self.send_header("Content-Type", "text/plain; charset=utf-8")
                self.send_header("Content-Length", str(len(body_bytes)))
                self.end_headers()
                if send_body:
                    self.wfile.write(body_bytes)

            def do_GET(self):  # noqa: N802
                self._serve(send_body=True)

            def do_HEAD(self):  # noqa: N802
                self._serve(send_body=False)

            def log_message(self, fmt, *args):  # silence access log
                return

        server = HTTPServer(("127.0.0.1", 0), Handler)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        host, port = server.server_address
        return server, host, port

    def _run_compile(
        self, main_wdl_path: str, env_value: Union[str, None]
    ) -> Tuple[int, str]:
        jar_path = os.path.join(
            self._context.repo_root_dir, f"dxCompiler-{self._context.version}.jar"
        )
        env = dict(os.environ)
        if env_value is None:
            env.pop(_BEARER_TOKENS_ENV_VAR, None)
        else:
            env[_BEARER_TOKENS_ENV_VAR] = env_value
        cmd = [
            "java",
            "-jar",
            jar_path,
            "compile",
            main_wdl_path,
            "-compileMode",
            "IR",
            "-quiet",
        ]
        self._context.logger.info(f"BearerAuthImport: COMPILE COMMAND {' '.join(cmd)}")
        proc = sp.run(cmd, env=env, capture_output=True, text=True)
        return proc.returncode, (proc.stdout or "") + (proc.stderr or "")

    @staticmethod
    def _check_outcome(
        rc: int,
        combined: str,
        expect_success: bool,
        must_include: List[str],
    ) -> Tuple[bool, str]:
        if expect_success:
            if rc != 0:
                return False, f"expected compile success, got exit {rc}"
            return True, "compile succeeded as expected"
        if rc == 0:
            return False, "expected compile failure, got exit 0"
        missing = [s for s in must_include if s not in combined]
        if missing:
            return (
                False,
                f"expected failure output to contain {missing!r}",
            )
        return True, f"compile failed as expected (exit {rc})"
