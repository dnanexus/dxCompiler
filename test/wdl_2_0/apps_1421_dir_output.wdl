version development

task apps_1421_task_01 {
  input {}

  command <<<
    mkdir folderoutput
    echo hello > folderoutput/hello.txt
  >>>
  output {
    Directory outdir = "folderoutput/"
  }
}

task apps_1421_task_02 {
  input {
    Directory indir
  }

  command <<<
    mkdir folderoutput_02
    cp ~{indir}/* folderoutput_02/
  >>>
  output {
    Directory outdir = "folderoutput_02/"
  }
}

workflow apps_1421_dir_output {
  input {}
  call apps_1421_task_01 {}
  call apps_1421_task_02 {
    input: indir = apps_1421_task_01.outdir
  }
  output {
    Directory outdir = apps_1421_task_02.outdir
  }
}
