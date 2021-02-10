version 1.0

workflow nested_pairs {
  input {
    Array[Pair[String, Pair[String, String]]] manifest
  }

  scatter (item in manifest) {
    call testtask {
      input:
        s1 = item.right.left,
        s2 = item.right.right,
    }

    String name = item.left
  }

  output {
    Array[String] messages = testtask.message
    Array[String] names = name
  }
}

task testtask {
  input {
    String s1
    String s2
  }

  command {
    echo "hello bam"
  }

  output {
    String message = "~{s1} ~{s2}"
  }
}
