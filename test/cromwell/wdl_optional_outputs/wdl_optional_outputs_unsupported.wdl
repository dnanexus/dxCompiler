version 1.0

struct Files {
    File yes
    File? maybe
}

task unsupported_pairs {
    command {
        touch yes
        # no no touching
        # touch no
    }
    output {
        # Even though the nonexistent file is being assigned to the optional File? in the Pair,
        # Cromwell can't currently work out optionality in Pairs with both optional and non-optional files.
        Pair[File?, File] one_optional = ("no", "yes")
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
}


task unsupported_structs {
    command {
        touch yes
        # no no touching
        # touch no
    }
    output {
        # Cromwell currently cannot work out optionality in structs with both non-optional and optional files.
        Files files = object { yes: "yes", maybe: "no" }
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
}

workflow wdl_optional_outputs_unsupported {
    call unsupported_pairs
    call unsupported_structs
}
