version 1.0

import "subworkflow_with_task.wdl" as subworkflow

workflow apps_700 {

    input {
        Boolean input_bool
        String input_string
        Int input_int = 1111
        Int input_int_not_default_in_subworkflow2 = 11111
    }

    call subworkflow.subworkflow_with_task as subwf700 {
        input:
            bool_input = input_bool,
            string_input = input_string,
            int_default = input_int,
            int_not_default = input_int_not_default_in_subworkflow2
    }

    output {
        Int integer_result = subwf700.res_int1
        Int integer_result_not_default = subwf700.res_int2
        String string_result = subwf700.res_string
        Boolean boolean_result = subwf700.res_boolean
        String subworkflow_string_input = subwf700.res_string_subtask_input
    }
}