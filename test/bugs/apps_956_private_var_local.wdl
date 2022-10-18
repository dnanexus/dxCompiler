version 1.1

task apps_956_private_var_local_t1 {
  input {
    File file_in
  }
  File file_var = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G0XfZ8Q0yzZV60XJJ4x4pg5z"

  command <<<
  set -euxo pipefail
  echo "File from input should be localized."
  cat ~{file_in}
  echo "File from private var should be localized, also."
  cat ~{file_var}
  >>>

  runtime {
    container: "ubuntu:20.04"
  }
}

workflow apps_956_private_var_local {
  call apps_956_private_var_local_t1 {
    input:
      file_in = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV6180yzZyvZ124KB0q46P"
  }
}
