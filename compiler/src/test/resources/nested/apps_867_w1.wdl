version 1.1

import "apps_867_w2.wdl" as lib

task t1_1 {
  command <<<
  set -euxo pipefail
  echo "Docker image file in same project"
  cat /etc/os-release
  >>>
}

workflow w1 {
  call t1_1
  call lib.apps_867_w2
}
