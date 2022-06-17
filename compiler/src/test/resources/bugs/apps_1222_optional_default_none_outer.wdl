version 1.1


import "apps_1222_optional_default_none_inner.wdl" as apps_1222_inner
workflow apps_1222_optional_default_none_outer {
  input {
    Int arg2 = 2
  }
  call apps_1222_inner.apps_1222_optional_default_none_inner {
    input:
      arg1=1
  }

  output {
    String outer_out = "Hello World"
  }
}
