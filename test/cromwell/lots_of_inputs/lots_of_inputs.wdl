task do_nothing {
  Array[File] f

  command {
    echo "no-op"
  }
  output {
    String o = read_string(stdout())
  }
  runtime {
    docker: "dx://file-GJ7q5KQ0yzZv7GQzJjX2x9Pb"
  }
}

task make_array {
  Int n
  command {
    python <<CODE
    for i in range(${n}):
      filename = 'file-' + str(i)
      with open(filename, 'w') as fp:
        fp.write(filename)
      print(filename)
    CODE
  }
  output {
    Array[File] a = glob("file-*")
  }
  runtime {
    docker: "dx://file-GJ7q5KQ0yzZv7GQzJjX2x9Pb"
  }
}

workflow lots_of_inputs {
  Int how_many_is_lots
  call make_array { input: n = how_many_is_lots }
  call do_nothing { input: f = make_array.a }

  output {
    Int out_count = length(make_array.a)
    String nothing_out = do_nothing.o
  }
}
