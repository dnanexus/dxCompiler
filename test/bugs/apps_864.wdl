version 1.0
workflow apps_864 {

    input {
        String sample
        String old_fastq_folder
        String R1_or_R2
    }

    call find_fastq_from_previous {
        input:
            sample = sample,
            old_fastq_folder = old_fastq_folder,
            R1_or_R2 = R1_or_R2
    }

    output {
        Array[File] fq = find_fastq_from_previous.fq
    }
}

task find_fastq_from_previous {
    input {
        String sample
        String old_fastq_folder
        String R1_or_R2
    }
    command <<<
        project=$DX_PROJECT_CONTEXT_ID
        folder=$project://~{old_fastq_folder}

        fq=$(dx find data --name "~{sample}~{R1_or_R2}.fastq.gz" --path ${folder} --norecurse --brief)
        for file_id in $fq
            do
                echo "dx://${file_id}"
            done
    >>>
    output {
        Array[File] fq = read_lines(stdout())
    }
}