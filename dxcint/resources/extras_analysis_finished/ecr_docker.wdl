version 1.0

task ecr_docker {
  input {
    String s
  }
  command <<<
  echo ~{s}
  >>>
  runtime {
    docker: "263799606133.dkr.ecr.us-east-1.amazonaws.com/dxcompiler_integration"
  }
  output {
    String sout = s
  }
}
