version 1.0

import "subworkflow_with_task.wdl" as subworkflow

workflow apps_700 {

    input {
        Boolean input_bool
        String input_string
        Int input_int = 1111
        Int input_int_not_default_in_subworkflow2 = 11111
    }

    call subworkflow.subworkflow as subwf {
        input:
            bool_input = input_bool,
            string_input = input_string,
            int_default = input_int,
            int_not_default = input_int_not_default_in_subworkflow2
    }

    output {
        Int integer_result = subwf.res_int1
        Int integer_result_not_default = subwf.res_int2
        String string_result = subwf.res_string
        Boolean boolean_result = subwf.res_boolean
        String subworkflow_string_input = subwf.res_string_subtask_input
    }
}