version 1.0

workflow apps_1858_inner {
    input {
        String run_name
    }

    call apps_1858_task_01_inner  {input: num = 1}
    output {
        String wf_out_01 = apps_1858_task_01_inner.out
    }
}

task apps_1858_task_01_inner {
    input {
        Int num
    }
    command <<<
        echo "Hello World"
    >>>
    output {
        String out = "Hello"
    }
}

