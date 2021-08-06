version 1.0
task diskspace_exhauster {
  input {
    Int count
  }
  command <<<
    set -euo pipefail
    rc=0
    mkdir tmp_folder
    for i in $(seq 1 ~{count}); do
      dd if=/dev/zero bs=1G count=1 of=tmp_folder/${i} || rc=1
      df -h /
      if [ "$rc" == 1 ]; then
        yes >tmp_folder/yes
        df -h /
        exit 1
      fi
    done
  >>>
  output {
    Array[File] outFiles = glob("tmp_folder/*")
  }
  runtime {
    docker: "ubuntu:20.04"
    dx_instance_type: "mem3_ssd1_x2"
  }
}
