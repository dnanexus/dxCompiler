version 1.0

workflow apps_1647_dx_restart {
    input {
        String one_input = "foo"
    }
    if (true) {
        call apps_1647_dx_restart_t01 as t01_conditional { input: test_in = one_input }
    }
    output {
        String? t01_conditional_result = t01_conditional.test_out
    }
}


task apps_1647_dx_restart_t01 {
    input {
        String test_in
    }
    command <<<
        echo "Hello ${test_in}"
    >>>
    output {
        String test_out = "Hello ${test_in}"
    }
    runtime {
        dx_restart: object {
         max: 5,
         errors: object {
             UnresponsiveWorker: 3,
             ExecutionError: 3,
            }
        }
    }
}
