# Bearer-token authenticated WDL imports

End-to-end test for the `DXCOMPILER_WDL_IMPORT_BEARER_TOKENS` env var.

Files:
- `imported.wdl` — served by a local HTTP server behind a static Bearer token.
- `main.wdl.template` — workflow that imports the protected document via an
  HTTP URL. `{URL}` is replaced at test time with the server's address.

Driver: `scripts/run_tests.py --test bearer_auth` (or the equivalent
pseudo-test name `--test bearer_auth_import`) brings up the server, renders the
template, and invokes the dxCompiler JAR (with `-compileMode IR`, no platform
required) under three configurations:

1. no token → expect compile failure with HTTP 401 (no credentials)
2. wrong token → expect compile failure with HTTP 403 (credentials rejected)
3. correct token → expect compile success
