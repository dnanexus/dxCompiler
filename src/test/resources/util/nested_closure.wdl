version 1.0

# This workflow demonstrates the need to sort elements by dependency
# order, rather than process them in order of appearance, when computing
# the inputs and outputs for each block of statements. For example, if
# we processed elements in order of appearance, `foo.out` would be
# classified as an input rather than an output, since it is referenced
# prior to its definition.
workflow nested_closure {
    input {
        Int i
    }

    # 1 - input=[i], output=[indexes]
    Array[Int] indexes = range(i)

    #2 - input=[i, x], output=[indexes]
    scatter (x in indexes) {
        # 3 - input=[i, x], output=[indexes, s]
        String s = "item ~{x}"
        # 5 - input=[i, x], output=[indexes, s, foo.out, t]
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
        Array[Array[String]] out2 = out
    }
}

task foo {
    input {
        String s
    }
    command {}
    output {
        String out = s
    }
}