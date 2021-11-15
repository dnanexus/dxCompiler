task t {
  command {
    echo hi
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow monitoring_script_localization_failure {
  call t
}
