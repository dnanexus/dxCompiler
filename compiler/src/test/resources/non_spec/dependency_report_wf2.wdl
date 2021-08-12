version 1.1

# Imported in dependency_report_wf1.wdl
# Input files are stored in dxCompiler_playground:/test_data/dependency_report.
# Their contents are not important; they only need to exist as platform files.

task dependency_report_t1 {
  input {
    File default = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV6180yzZyvZ124KB0q46P"
    File pmeta
    File body1
  }

  parameter_meta {
    pmeta: {
      suggestions: [
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV61Q0yzZq60Jj4K5vfG92"
      ],
      choices: [
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV61Q0yzZq60Jj4K5vfG92",
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV61j0yzZgf6JQKxP4gQ3Y"
      ]
    }
  }

  File body2 = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV6200yzZq60Jj4K5vfG94"

  command <<<
  set -euxo pipefail

  cat ~{default} > output.txt
  cat ~{pmeta} >> output.txt
  cat ~{body1} >> output.txt
  cat ~{body2} >> output.txt
  >>>

  output {
    File out = "output.txt"
  }

  runtime {
    docker: "ubuntu:20.04"
    dx_instance_type: "mem1_ssd2_x4"
  }
}

task dependency_report_t2 {
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

workflow dependency_report_wf2 {
  input {
    File default = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV6180yzZyvZ124KB0q46P"
    File pmeta
  }

  parameter_meta {
    pmeta: {
      suggestions: [
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV61Q0yzZq60Jj4K5vfG92"
      ],
      choices: [
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV61Q0yzZq60Jj4K5vfG92",
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV61j0yzZgf6JQKxP4gQ3Y"
      ]
    }
  }

  File body1 = "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV6280yzZq60Jj4K5vfG96"

  call dependency_report_t1 {
    input:
      default = default,
      pmeta = pmeta,
      body1 = body1
  }

  call dependency_report_t2 {
    input:
      files = [
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV62Q0yzZjqf4F4KJvv1gG",
        "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:file-G4BV62j0yzZyvZ124KB0q46V"
      ]
  }

  output {
    File dependency_report_t1_out = dependency_report_t1.out
    File dependency_report_t2_out = dependency_report_t2.out
  }
}