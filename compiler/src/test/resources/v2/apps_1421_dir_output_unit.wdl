version development

task apps_1421_unit_task_01 {
  input {}

  command <<<
    mkdir folderoutput
    echo hello > folderoutput/hello.txt
  >>>
  output {
    Directory outdir = 'folderoutput/'
  }
}

workflow apps_1421_dir_output_unit {
  input {}
  call apps_1421_unit_task_01 {}
  output {
    Directory outdir = apps_1421_unit_task_01.outdir
  }
}
