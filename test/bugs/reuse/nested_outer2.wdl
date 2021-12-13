version 1.0
import "nested_inner.wdl" as nested_inner

workflow nested_outer2 {
    input {
        String s
    }
    call nested_inner.nested_inner { input: s = "~{s+s}" }
    call nested_inner.nested_inner as nested_inner2 { input: s = "~{s}" }
    output {
        File out_inner = nested_inner.nested_inner_wf_out
        File out_inner2 = nested_inner2.nested_inner_wf_out
    }
}