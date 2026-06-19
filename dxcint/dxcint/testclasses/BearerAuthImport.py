import os
import shutil
import ssl
import subprocess as sp
import tempfile
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Dict, List, Optional, Tuple, Union

from dxcint.Context import Context, ContextEmpty
from dxcint.RegisteredTest import RegisteredTest


# Token expected by the in-process HTTPS server. Mirrors the value used by the
# scalatest integration spec; keep these in sync if it ever needs to change.
_EXPECTED_TOKEN = "secret-token-abc"

# Env variable that dxCompiler reads to pick up per-domain Bearer tokens for
# WDL imports. Mirrored from dx.core.languages.wdl.WdlImportHttpAuth.TokensEnvVar.
_BEARER_TOKENS_ENV_VAR = "DXCOMPILER_WDL_IMPORT_BEARER_TOKENS"

# Filename of the WDL served by the local HTTP(S) servers. Lives next to the
# main workflow template in the test fixture directory.
_PROTECTED_WDL_FILENAME = "imported.wdl"

# Placeholder in the main WDL fixture that is replaced with the live server
# URL at test time.
_URL_PLACEHOLDER = "{URL}"

# Loopback host the servers bind to. Used verbatim in the import URL, in the
# per-domain token map, and in the certificate SAN so JVM hostname
# verification succeeds.
_HOST = "127.0.0.1"

# Password protecting the generated Java truststore.
_TRUSTSTORE_PASSWORD = "changeit"


class BearerAuthImport(RegisteredTest):
    """Local-only test that verifies how dxCompiler follows http(s) WDL imports
    and honours `DXCOMPILER_WDL_IMPORT_BEARER_TOKENS`.

    It exercises two supported flows plus the Bearer-token security gate:

    * a public document over plain HTTP (no credentials) -> compile succeeds;
    * a document over HTTPS behind a Bearer token, served with a self-signed
      certificate that is trusted by the JVM via an injected truststore:
        - no token   -> compile fails with HTTP 401;
        - wrong token -> compile fails with HTTP 403;
        - right token -> compile succeeds.

    Bearer tokens are intentionally only attached to HTTPS requests, so the
    plain-HTTP flow never sends credentials.

    The certificate, server key, and Java truststore are all generated at
    runtime (openssl + keytool); nothing is checked into the repository.

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
        with tempfile.TemporaryDirectory(prefix="dxcint-bearer-auth-") as workdir:
            truststore = self._generate_tls_material(workdir)

            http_server = self._start_server(require_auth=False)
            https_server = self._start_server(
                require_auth=True, ssl_context=self._server_ssl_context()
            )
            try:
                http_url = self._url("http", http_server)
                https_url = self._url("https", https_server)

                # (label, import url, token env value, inject truststore, expect
                #  success, required output fragments on failure)
                scenarios: List[
                    Tuple[str, str, Optional[str], bool, bool, List[str]]
                ] = [
                    (
                        "public import over HTTP (no auth)",
                        http_url,
                        None,
                        False,
                        True,
                        [],
                    ),
                    (
                        "HTTPS import, no token configured",
                        https_url,
                        None,
                        True,
                        False,
                        [
                            "HTTP 401 Unauthorized",
                            "ensure the credentials are provided",
                        ],
                    ),
                    (
                        "HTTPS import, wrong token configured",
                        https_url,
                        f"{_HOST}:wrong-token",
                        True,
                        False,
                        [
                            "HTTP 403 Forbidden",
                            "the provided credentials are invalid or lack the required permissions",
                        ],
                    ),
                    (
                        "HTTPS import, correct token configured",
                        https_url,
                        f"{_HOST}:{_EXPECTED_TOKEN}",
                        True,
                        True,
                        [],
                    ),
                ]

                failures: List[str] = []
                for (
                    label,
                    url,
                    env_value,
                    inject_truststore,
                    expect_success,
                    must_include,
                ) in scenarios:
                    main_path = os.path.join(workdir, "main.wdl")
                    with open(main_path, "w") as f:
                        f.write(self._main_template.replace(_URL_PLACEHOLDER, url))
                    rc, combined = self._run_compile(
                        main_path,
                        env_value,
                        truststore if inject_truststore else None,
                    )
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
                for server in (http_server, https_server):
                    server.shutdown()
                    server.server_close()

    # --- TLS material -----------------------------------------------------

    def _generate_tls_material(self, workdir: str) -> str:
        """Generate a self-signed cert + key (openssl) and a Java truststore
        (keytool) trusting it. Returns the truststore path.

        The cert lives at ``<workdir>/cert.pem`` / ``<workdir>/key.pem`` so the
        HTTPS server can present it; the truststore lets the JVM running
        dxCompiler trust that cert when fetching the import.
        """
        openssl = shutil.which("openssl")
        if openssl is None:
            raise RuntimeError("openssl is required to generate the test certificate")
        keytool = shutil.which("keytool")
        if keytool is None:
            java = shutil.which("java")
            if java:
                keytool = os.path.join(os.path.dirname(java), "keytool")
        if not keytool or not os.path.exists(keytool):
            raise RuntimeError("keytool is required to build the test truststore")

        cert_pem = os.path.join(workdir, "cert.pem")
        key_pem = os.path.join(workdir, "key.pem")
        truststore = os.path.join(workdir, "truststore.p12")

        sp.run(
            [
                openssl,
                "req",
                "-x509",
                "-newkey",
                "rsa:2048",
                "-keyout",
                key_pem,
                "-out",
                cert_pem,
                "-days",
                "1",
                "-nodes",
                "-subj",
                f"/CN={_HOST}",
                "-addext",
                f"subjectAltName=IP:{_HOST},DNS:localhost",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        sp.run(
            [
                keytool,
                "-importcert",
                "-noprompt",
                "-alias",
                "bearer",
                "-file",
                cert_pem,
                "-keystore",
                truststore,
                "-storetype",
                "PKCS12",
                "-storepass",
                _TRUSTSTORE_PASSWORD,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        self._cert_pem = cert_pem
        self._key_pem = key_pem
        return truststore

    def _server_ssl_context(self) -> ssl.SSLContext:
        ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        ctx.load_cert_chain(certfile=self._cert_pem, keyfile=self._key_pem)
        return ctx

    # --- HTTP(S) server ---------------------------------------------------

    def _start_server(
        self, require_auth: bool, ssl_context: Optional[ssl.SSLContext] = None
    ) -> ThreadingHTTPServer:
        body_bytes = self._protected_bytes
        protected_path = "/" + _PROTECTED_WDL_FILENAME
        expected_token = _EXPECTED_TOKEN

        class Handler(BaseHTTPRequestHandler):
            def _auth_status(self) -> int:
                # 401 when no credentials are supplied, 403 when credentials are
                # supplied but do not match, 200 when they match.
                if not require_auth:
                    return 200
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

        server = ThreadingHTTPServer((_HOST, 0), Handler)
        if ssl_context is not None:
            server.socket = ssl_context.wrap_socket(server.socket, server_side=True)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        return server

    @staticmethod
    def _url(scheme: str, server: ThreadingHTTPServer) -> str:
        port = server.server_address[1]
        return f"{scheme}://{_HOST}:{port}/{_PROTECTED_WDL_FILENAME}"

    # --- compile ----------------------------------------------------------

    def _run_compile(
        self,
        main_wdl_path: str,
        env_value: Optional[str],
        truststore: Optional[str],
    ) -> Tuple[int, str]:
        jar_path = os.path.join(
            self._context.repo_root_dir, f"dxCompiler-{self._context.version}.jar"
        )
        env = dict(os.environ)
        if env_value is None:
            env.pop(_BEARER_TOKENS_ENV_VAR, None)
        else:
            env[_BEARER_TOKENS_ENV_VAR] = env_value
        # JVM args must precede -jar. Injecting the truststore makes the JVM
        # trust the self-signed cert; IR compile is local-only, so replacing the
        # default truststore for this invocation is harmless.
        jvm_args: List[str] = []
        if truststore is not None:
            jvm_args = [
                f"-Djavax.net.ssl.trustStore={truststore}",
                "-Djavax.net.ssl.trustStoreType=PKCS12",
                f"-Djavax.net.ssl.trustStorePassword={_TRUSTSTORE_PASSWORD}",
            ]
        cmd = (
            ["java"]
            + jvm_args
            + [
                "-jar",
                jar_path,
                "compile",
                main_wdl_path,
                "-compileMode",
                "IR",
                "-quiet",
            ]
        )
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
