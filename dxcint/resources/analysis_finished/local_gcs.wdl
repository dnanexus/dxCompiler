task makeGcsFileInJes {
  command {
    echo "I am writing stuff" > out
  }

  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }

  output {
    File outputFile = "out"
  }
}

task readGcsFileFromLocal {
  File gcsFile

  command {
    cat ${gcsFile}
  }

  runtime {
    backend: "Local"
  }

  output {
    String result = read_string(stdout())
  }
}

workflow local_gcs {
  call makeGcsFileInJes
  call readGcsFileFromLocal { input: gcsFile = makeGcsFileInJes.outputFile }
}
