version 1.0
import "apps_1175_nested_inner.wdl" as nested_inner

workflow apps_1175_nested_wf_intermediate_outputs {
  input {
    String s
  }
  Boolean a = true
  if (a) {
    call nested_inner.test_inner1 as call_1 { input: test_in = s }
    call nested_inner.test_inner1 as call_2 { input: test_in = s }
  }
  output {
    File? out_inner = call_1.test_out
  }
}