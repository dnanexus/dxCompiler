task subdirTask {
  command {
    mkdir subdir
    cd subdir
    echo "I'm in a subdirectory !" > subFile
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
  output {
    File outputFile = "subdir/subFile"
  }
}

workflow subdirWorkflow {
  call subdirTask
}
