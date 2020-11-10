version 1.0

workflow nested_closure {
    input {
        Int i
    }

    # 1 - input=[i], output=[indexes]
    # input=[i], output=[indexes]
    Array[Int] indexes = range(i)

    #2 - input=[i, x], output=[indexes]
    # input=[i, x], output=[indexes]
    scatter (x in indexes) {
        # 3 - input=[i, x], output=[indexes, s]
        # input=[i], output=[indexes, s]
        String s = "item ~{x}"
        # 5 - input=[i, x], output=[indexes, s, foo.out, t]
        # input=[i, foo.out], output=[indexes]
        String t = foo.out

        # 4 - input=[i, x], output=[indexes, s, foo.out]
        call foo {
            input: s = s
        }

        # 10 - input=[i, x, z, j], output=[indexes, s, foo.out, t, string, n, out]
        Array[String] out = n

        # 6 - input=[i, x, z], output=[indexes, s, foo.out, t]
        scatter (z in [s, t]) {
            # 9 - input=[i, x, z, j], output=[indexes, s, foo.out, t, string, n]
            String n = "~{sep=',' string}"

            # 7 - input=[i, x, z, j], output=[indexes, s, foo.out, t]
            scatter (j in [1,2]) {
                # 8 - input=[i, x, z, j], output=[indexes, s, foo.out, t, string]
                String string = "~{j} ~{t}"
            }
        }
    }

    output {
        Array[String] out2 = out
    }
}