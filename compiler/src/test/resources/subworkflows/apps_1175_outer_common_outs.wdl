version 1.0
import "apps_1175_inner_common_outs.wdl" as nested_inner_co

workflow apps_1175_outs_not_uploaded_outer {
  input {
    String s
  }
  call nested_inner_co.nested_inner_co as ni_co { input: s = "~{s+s}" }
  output {
    File out_inner = ni_co.nested_inner_wf_out
  }
}