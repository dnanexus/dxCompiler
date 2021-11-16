task noop {
  command {}
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow workflow_type_and_version_cwl {
  call noop
}
