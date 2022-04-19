version 1.0
task unzip_files {
    input {
        Array[File] zipped_files
    }
    command <<<
        for f in ~{sep=' ' zipped_files}; do
            gunzip "${f}"
            f_unzipped="${f%%.gz}"
            file_id=$(dx upload "${f_unzipped}" --brief --tag unzipped --project ${DX_WORKSPACE_ID} --wait)
            # clean up the unzipped file to save local disk space
            rm "${f_unzipped}"
            echo "dx://${DX_WORKSPACE_ID}:${file_id}"
        done
    >>>
    output {
        Array[File] files_out = read_lines(stdout())
    }
}