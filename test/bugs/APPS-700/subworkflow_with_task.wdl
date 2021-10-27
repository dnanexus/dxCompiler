version 1.0
import "test_task.wdl" as test_task

workflow subworkflow {

    input {
        String string_input = "test-string-2"
        Boolean bool_input = true
        Int int_default = 2222
        Int int_not_default
    }
    String test_string = string_input
    call test_task.test_task as subtask {
        input:
            string = test_string,
            boolean = bool_input,
            integer1 = int_default,
            integer2 = int_not_default
    }


    output {
        Int res_int1 = subtask.integer1_res
        Int res_int2 = subtask.integer2_res
        String res_string = subtask.string_res
        Boolean res_boolean = subtask.boolean_res
        String res_string_subtask_input = string_input
    }
}