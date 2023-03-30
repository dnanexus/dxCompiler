version 1.1

task apps_1341_t01 {
    input {
        String inp_00
        String? inp_01
    }

    command <<<
    echo "Hello World"
    >>>

    output {
        String report_html = "placeholder_00"
        String stats_txt = "placeholder_01"
    }
}

workflow apps_1341 {
    input {
        String wf_inp_00
        String? wf_inp_01
    }

    String? new_wf_inp_01 = if defined(wf_inp_01) then wf_inp_01 else None
    call apps_1341_t01 as t01_0 {
        input:
            inp_00 = wf_inp_00,
            inp_01 = wf_inp_01
      }

    call apps_1341_t01 as t01_1 {
        input:
            inp_00 = if defined(wf_inp_00) then wf_inp_00 else wf_inp_00,
            inp_01 = wf_inp_01
    }

    call apps_1341_t01 as t01_2 {
        input:
            inp_00 = wf_inp_00,
            inp_01 = if defined(wf_inp_01) then wf_inp_01 else None
    }

    call apps_1341_t01 as t01_3 {
        input:
            inp_00 = wf_inp_00,
            inp_01 = new_wf_inp_01
    }
    
    output {
        String report0 = t01_0.report_html
        String report1 = t01_1.report_html
        String report2 = t01_2.report_html
        String report3 = t01_3.report_html
    }
}