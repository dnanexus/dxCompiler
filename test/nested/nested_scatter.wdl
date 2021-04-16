version 1.0

workflow nested_scatter {
    Array[Int] ints1 = [1, 2, 3]
    Array[Int] ints2 = [10, 11]

    scatter (x in ints1) {
        scatter (y in ints2) {
            call xxx_add{ input: a = x, b = y }

            Int i = ints1[x]
            if (i < 3) {
                String s = "hello ${x} ${y}"
            }
        }
    }

    output {
        Array[Array[Int]] result = xxx_add.result
        Array[String] strings = select_all(flatten(s))
    }
}


task xxx_add {
    input {
        Int a
        Int b
    }
    command {}
    output {
        Int result = a + b
    }
}
