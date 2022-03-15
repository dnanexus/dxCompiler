version 1.0

workflow bam_chrom_counter {
    input {
        File bam
    }

    call slice_bam {
        input : bam = bam
    }
    scatter (slice in slice_bam.slices) {
        call count_bam {
            input: bam = slice
        }
    }
    output {
        File bai = slice_bam.bai
        Array[Int] count = count_bam.count
    }
}

task slice_bam {
    input {
        File bam
        Int num_chrom = 22
    }
    command <<<
    set -ex
    samtools index ~{bam}
    mkdir slices/
    for i in `seq ~{num_chrom}`; do
        samtools view -b ~{bam} -o slices/$i.bam chr$i
    done
    >>>
    runtime {
        docker: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }
    output {
        File bai = "~{bam}.bai"
        Array[File] slices = glob("slices/*.bam")
    }
}

task count_bam {
    input {
        File bam
    }

    command <<<
        samtools view -c ~{bam}
    >>>
    runtime {
        docker: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }
    output {
        Int count = read_int(stdout())
    }
}