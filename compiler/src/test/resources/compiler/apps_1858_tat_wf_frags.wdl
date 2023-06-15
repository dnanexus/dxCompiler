version 1.0


import "apps_1858_inner.wdl" as inner

workflow apps_1858_tat_wf_frags {
    input {
        String run_name
    }

    if (true) {
        call apps_1858_task_01  {input: num = 1}
    }
    call apps_1858_task_02  {input: num = 2}

    if (true) {
        call inner.apps_1858_inner {input: run_name = run_name}
    }

    output {
        String wf_out_01 = apps_1858_task_01.out
        String wf_out_02 = apps_1858_task_02.out
    }
}

task apps_1858_task_01 {
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

task apps_1858_task_02 {
    input {
        Int num
    }
    command <<<
        echo "Hello World 2"
    >>>
    output {
        String out = "Hello 2"
    }
}
