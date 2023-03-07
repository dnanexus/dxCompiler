version 1.0

workflow apps_1473_scatter_empty_arr {
    Array[String] a = []
    scatter (item in a) {
        call apps_1473_scatter_empty_arr_t01
    }
    output {
        Array[Array[String]] somethings = apps_1473_scatter_empty_arr_t01.something
    }
}

task apps_1473_scatter_empty_arr_t01 {
    command <<<
        echo "Hello World!"
    >>>
    output {
        Array[String] something = ["thing 1", "thing 2"]
    }
}