version 1.0

workflow apps_1351 {
    call aa
    output {
        Array[Pair[String, String]] zipped = zip(aa.a, aa.b)
    }
}

task aa {
    command <<<
        echo a > a.txt
        echo b > b.txt
    >>>
    output {
        Array[String] a = read_lines("a.txt")
        Array[String] b = read_lines("b.txt")
    }
}