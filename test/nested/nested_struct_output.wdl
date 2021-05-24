version 1.0

import "nested_wf.wdl"

workflow nested_struct_output {
  input {
    Boolean b
    String x
  }

  if (b) {
    call nested_wf.wf as nwf {
      input: x = x
    }
    String s = nwf.result.s
  }

  output {
    String? out = s
  }
}