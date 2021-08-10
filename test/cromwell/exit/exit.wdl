task exitTask {
  command {
    exit 5
  }
  runtime {
    docker: "ubuntu:latest"
    returnCodes: 5
  }
}

workflow exitWorkflow {
  call exitTask
}
