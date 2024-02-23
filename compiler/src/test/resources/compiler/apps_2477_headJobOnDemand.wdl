version 1.0

workflow apps_2477_on_demand_wf {
    input {
        String run_name
    }

    call apps_2477_task_01  {input: num = 1}
    call apps_2477_task_02  {input: num = 2}

    output {
        String wf_out_01 = apps_2477_task_01.out
        String wf_out_02 = apps_2477_task_02.out
    }
}

task apps_2477_task_01 {
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

task apps_2477_task_02 {
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
