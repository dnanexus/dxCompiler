version 1.1

# Input files are stored in dxCompiler_playground:/unit_tests/dependency_report.
# Their contents are not important; they only need to exist as platform files.

import "dependency_report_wf2.wdl" as lib

task dependency_report_t3 {
  File body3 = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV6100yzZz17bx4JkQkybb"

  command <<<
  set -euxo pipefail

  cat ~{body3} >> output.txt
  >>>

  output {
    File out = "output.txt"
  }

  runtime {
    docker: "alpine:3.14"
    dx_instance_type: "mem1_ssd2_x4"
  }
}

workflow dependency_report_wf1 {
  input {
    File default = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV6180yzZyvZ124KB0q46P"
    File pmeta
  }

  call lib.dependency_report_wf2 as sub {
    input:
      default = default,
      pmeta = pmeta
  }

  call dependency_report_t3

  output {
    File dependency_report_t1_out = sub.dependency_report_t1_out
    File dependency_report_t2_out = sub.dependency_report_t2_out
    File dependency_report_t3_out = dependency_report_t3.out
  }
}