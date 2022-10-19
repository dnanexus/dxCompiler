version 1.1

import "global_wf_from_wdl_sub.wdl" as lib

task global_wf_from_wdl_t1 {
  command <<<
  set -euxo pipefail
  echo "Task 1"
  >>>

  runtime {
    container: "ubuntu:20.04"
  }
}

task global_wf_from_wdl_t2 {
  input {
    String str
  }

  command <<<
  set -euxo pipefail
  echo "Task 2"
  echo "Input: ~{str}"
  >>>

  runtime {
    container: "ubuntu:20.04"
  }
}

task global_wf_from_wdl_t3 {
  command <<<
  set -euxo pipefail
  echo "Docker image file as container"
  cat /etc/os-release
  >>>

  runtime {
    # Bundled dependency due to Docker image file
    # test_data/ubuntu_20_04.tar.gz
    container: "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G51b4F8016p804xB309zVFpB"
  }
}

task global_wf_from_wdl_t4 {
  input {
    Array[File]+ files
    String? output_filename
  }

  command <<< >>>

  output {
    File out = "placeholder.txt"
  }

  meta {
    type: "native"
    id: "app-BZ9ZQzQ02VP5gpkP54b96pYY"
  }
}

task global_wf_from_wdl_t5 {
  command <<<
  set -euxo pipefail
  echo "Task 5"
  >>>

  runtime {
    docker: "ubuntu:20.04"
    dx_instance_type: "mem1_ssd2_x4"
  }
}

workflow global_wf_from_wdl {
  Boolean t = true
  Int i = 1

  # Inner applet should be bundled by outer frag applet
  if (t) {
    call global_wf_from_wdl_t1
  }

  # Inner applet should be bundled by outer frag applet
  call global_wf_from_wdl_t2 {
    input:
      str = "~{if t then '${1 + i}' else 0}"
  }

  # Docker image file should be bundled
  call global_wf_from_wdl_t3

  # Native app should be in description
  call global_wf_from_wdl_t4 {
    input:
      files = [
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-GJ87YJ00yzZq4KJ51KF966fB",
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-GJ7VKG00yzZgfZzb1JgY2K2Q"
      ]
  }

  # Hard-coded dx instance type should be in description
  call global_wf_from_wdl_t5

  call lib.global_wf_from_wdl_sub
}
