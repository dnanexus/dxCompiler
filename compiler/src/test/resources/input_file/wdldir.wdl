version development

task foldertest {
  input {
    Directory WorkingDir
  }

  command {
    set -euxo pipefail
    ls -l ~{WorkingDir}
  }
  hints {}
  output {}
}

workflow folderrun {
  input {
    Directory WorkingDir
  }

  call foldertest {
    input: 
      WorkingDir=WorkingDir
  }
}