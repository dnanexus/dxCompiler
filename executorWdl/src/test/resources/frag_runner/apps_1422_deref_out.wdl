version 1.0

workflow apps_1422_deref {
    input {
        Array[String] ina1 = ["Hello", "Halo"]
        Array[String] ina2 = ["World", "Welt"]
    }
    Array[Pair[String, String]] pairs = zip(ina1, ina2)
    call apps_1422_deref_task_01 as t01_1 { input: t01_in = ina1[0]}
    String xx = pairs[0].left
    call apps_1422_deref_task_01 as t01_2 { input: t01_in = xx}
    output {
        String apps_1422_deref_out_01 = t01_1.t01_out
        String apps_1422_deref_out_02 = t01_2.t01_out
    }
}

task apps_1422_deref_task_01 {
    input {
        String t01_in
    }

    command <<<>>>
    output {
        String t01_out = t01_in
    }
}