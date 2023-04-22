task exitTask {
  command {
    exit 5
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    returnCodes: 5
  }
}

workflow exitWorkflow {
  call exitTask
}
