version 1.0

workflow basename_fields {
  input {
    File tool
  }

  call echo_file_root as root {
    input: f = basename(tool)
  }

  call echo_file_ext as ext {
    input: f = basename(tool)
  }

  output {
    File rootFile = root.out
    File extFile = ext.out
  }
}

task echo_file_root {
  input {
    String f
    String? name
  }

  command <<<
    filename=~{f}
    echo "${filename%.*}"
  >>>

  output {
    File out = stdout()
  }
}

task echo_file_ext {
  input {
    String f
    String? name
  }

  command <<<
    filename=~{f}
    echo "${filename##*.}"
  >>>

  output {
    File out = stdout()
  }
}