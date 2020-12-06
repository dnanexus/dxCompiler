version 1.0

# This is a version of the nested_closure workflow that doesn't have forward references.
workflow nested_closure_no_fw_ref {
  input {
    Int i
  }

  # 1 - input=[i], output=[indexes]
  Array[Int]+ indexes = [i, i + 1]

  #2 - input=[i, x], output=[indexes]
  scatter (x in indexes) {
    # 3 - input=[i, x], output=[indexes, s]
    String s = "item ~{x}"

    # 4 - input=[i, x], output=[indexes, s, foo.out]
    call foo {
      input: s = s
    }

    # 5 - input=[i, x], output=[indexes, s, foo.out, t]
    String t = foo.out

    # 6 - input=[i, x], output=[indexes, s, foo.out, t]
    # note that 'z' is never referenced, so it is not added
    # to the inputs
    scatter (z in [s, t]) {
      # 7 - input=[i, x, z, j], output=[indexes, s, foo.out, t]
      scatter (j in [1,2]) {
        # 8 - input=[i, x, z, j], output=[indexes, s, foo.out, t, string]
        String string = "~{j} ~{t}"
      }

      # 9 - input=[i, x, z, j], output=[indexes, s, foo.out, t, string, n]
      String n = "~{sep=',' string}"
    }

    # 10 - input=[i, x, z, j], output=[indexes, s, foo.out, t, string, n, out]
    Array[String]+ out = n
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