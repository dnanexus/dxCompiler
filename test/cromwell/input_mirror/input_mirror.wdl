task mkFile {
    command {
        echo "The sixth sheikh's sheep is sick"
    }
    output {
        File sick_sheep = stdout()
    }
    runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

task input_mirror {
    File inFile
    command {
        # NO-OP
    }
    output {
        File outFile = inFile
        String out = read_string(outFile)
    }
    runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

workflow mirror_mirror_on_the_wall {
    call mkFile
    call input_mirror as mirror1 { input: inFile = mkFile.sick_sheep }
    call input_mirror as mirror2 { input: inFile = mirror1.outFile }
}
