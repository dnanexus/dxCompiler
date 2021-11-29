version 1.0

workflow view_and_count {
    input {
        File bam
    }

    call slice_bam {
        input: bam = bam,
               size = "8k"
    }

    scatter (slice in slice_bam.slices) {
        call slice_bam as slice2 {
            input: bam = slice,
                   size="1k"
        }
        scatter ( slice3 in slice2.slices) {
            call identity {
                input: bam = slice3
            }
        }
    }

    output {
        Array[Array[File]] counts = identity.res
    }
}

task slice_bam {
    input {
        File bam
        String size
    }

    command <<<
        split -b~{size} -a6 --additional-suffix=.bam  ~{bam} bam
    >>>

    output {
        Array[File] slices = glob("bam*")
    }
}

task identity {
    input {
        File bam
    }
    command <<<
        cd .
    >>>

    output {
        File res = bam
    }
}