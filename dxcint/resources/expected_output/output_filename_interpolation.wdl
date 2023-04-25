task output_filename_interpolation_t01 {
    String outputPrefix = "Oops"
    String outputMain = "Upside"
    String outputSuffix = "YourHead"

    command {
        mkdir ${outputPrefix}
        echo a > ${outputMain}
        echo b > ${outputMain + "." + outputSuffix}
        echo c > ${outputPrefix}/${outputMain}.${outputSuffix}
    }

    output {
        String a = read_string(outputMain)
        String b = read_string("${outputMain}.${outputSuffix}")
        String c = read_string(outputPrefix + "/" + outputMain + "." + outputSuffix)
    }

    runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

workflow output_filename_interpolation {
    call output_filename_interpolation_t01
}
