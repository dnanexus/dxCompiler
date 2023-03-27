version 1.1

task t2_1 {
  command <<<
  set -euxo pipefail
  echo "Task 2.1"
  >>>

  runtime {
    container: "ubuntu:20.04"
  }
}

task t2_2 {
  input {
    String str
  }

  command <<<
  set -euxo pipefail
  echo "Task 2.2"
  echo "Input: ~{str}"
  >>>

  runtime {
    container: "ubuntu:20.04"
  }
}

task t2_3 {
  input {
    String str
  }

  command <<<
  set -euxo pipefail
  echo "Task 2.3"
  echo "Input: ~{str}"
  >>>

  runtime {
    container: "ubuntu:20.04"
  }
}

workflow apps_867_w2 {
  Boolean t = true
  Int i = 1

  if (t) {
    call t2_1
  }

  Array[String] arr = ["A1", "A2", "A3"]
  scatter (s in arr) {
    call t2_2 {
      input:
        str = s
    }
  }

  call t2_3 {
    input:
      str = "hello"
  }

  call t2_3 as t2_3_eval {
    input:
      str = "~{if t then '${1 + i}' else 0}"
  }
}
