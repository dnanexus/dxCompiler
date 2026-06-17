version 1.0

# Template workflow for the bearer-auth WDL import test.
# The {URL} placeholder is replaced with the address of a local HTTP server by
# the BearerAuthImport dxcint test class before invoking the compiler. This
# file is intentionally not parseable as-is — it is rendered into a temp dir at
# test time. See dxcint/dxcint/testclasses/BearerAuthImport.py.

import "{URL}" as lib

workflow bearer_auth_import {
  input {
    String name = "world"
  }
  call lib.hello { input: name = name }
  output {
    String greeting = hello.greeting
  }
}
