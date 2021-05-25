version 1.0

import "nested_wf.wdl"

workflow apps_612 {
  input {
    Boolean b
    String x
  }

  if (b) {
    call nested_wf.nested_wf as nwf {
      input: x = x
    }
    String s = nwf.result.s
  }

  output {
    String? out = s
  }
}