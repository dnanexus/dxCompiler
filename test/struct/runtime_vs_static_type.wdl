version 1.0

struct WordInfo {
  String word
  Int len
}

workflow runtime_vs_static_type {
    call rvs_opt_int { input: x = 14 }

    call rvs_opt_array { input: xa = [14,15,20] }

    WordStruct manitoba = object {
        word: "Manitoba",
        len: 8
    }

    call rvs_opt_struct { input : ao = [manitoba] }

    output {
        Int result = rvs_opt_int.result
        String result2 = rvs_opt_array.numbers
        String result3 = rvs_opt_struct.w
    }
}

task rvs_opt_int {
    input {
        Int? x
    }
    command {
        echo $(( ~{x} + 10 ))
    }
    output {
        Int result = read_int(stdout())
    }
}

task rvs_int {
  input {
    Int? x
  }
  command {
    echo $(( ~{x} + 10 ))
  }
  output {
    Int result = read_int(stdout())
  }
}

task rvs_opt_array {
    input {
        Array[Int?] xa
    }
    Array[Int] a = select_all(xa)
    command {
        echo ~{sep=',' a}
    }
    output {
        String numbers = read_string(stdout())
    }
}

task rvs_opt_struct {
    input {
        Array[WordStruct?]? ao
    }
    Array[WordStruct?] a = select_first([ao])
    Array[WordStruct] a2 = select_all(a)
    command{}
    output {
        String w = a2[0].word
    }
}
