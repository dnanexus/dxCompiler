version 1.1

task apps_623_t1 {
  command <<<
  set -euxo pipefail
  echo "Task 1"
  >>>

  runtime {
    container: "ubuntu:20.04"
  }
}

task apps_623_t2 {
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

task apps_623_t3 {
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

workflow apps_623_wf {
  Boolean t = true
  Int i = 1

  # Inner applet should be bundled by outer frag applet
  if (t) {
    call apps_623_t1
  }

  # Inner applet should be bundled by outer frag applet
  call apps_623_t2 {
    input:
      str = "~{if t then '${1 + i}' else 0}"
  }

  # Docker image file should be bundled
  call apps_623_t3
}
