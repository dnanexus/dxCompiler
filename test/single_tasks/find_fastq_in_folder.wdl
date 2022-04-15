version 1.0
task find_fastq_in_folder {
    input {
        String fastq_folder
    }
    command<<<
        project=$DX_PROJECT_CONTEXT_ID
        folder=$project://~{fastq_folder}
        fq=$(dx find data --name "*.fastq.gz" --path ${folder} --norecurse --brief)
        for file_id in $fq
        do
        echo "dx://${file_id}"
        done
    >>>
    output {
        Array[File] fq = read_lines(stdout())
    }
}