task singleFile {
  command {
    echo hello
  }
  output {
    File out = stdout()
  }
  runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

task listFiles {
  Array[File] manyIn
  command {
    cat ${sep=" " manyIn}
  }
  output {
    String result = read_string(stdout())
  }
  runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

workflow single_to_array_conversion {
  call singleFile
  call listFiles { input: manyIn = [ singleFile.out ] }
}
