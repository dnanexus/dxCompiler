version 1.0
import "nested_inner.wdl" as nested_inner

workflow nested_outer4 {
    input {
        String s
    }
    call nested_inner.nested_inner { input: s = "nested_outer4" }
    output {
        File out_inner = nested_inner.nested_inner_wf_out
    }
}