task write_map {
  Map[String, String] file_to_name
  command {
    cat ${write_map(file_to_name)}
  }
  output {
    String contents = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

task read_map {
  command <<<
    python <<CODE
    map = {'x': 500, 'y': 600, 'z': 700}
    print("\\n".join(["{}\\t{}".format(k,v) for k,v in map.items()]))
    CODE
  >>>
  output {
    Map[String, Int] out_map = read_map(stdout())
  }
  runtime {
    docker: "dx://file-GJ941b80yzZvGbK68zxQzB0B"
  }
}

workflow read_write_map {
  Map[String, String] map = {"f1": "alice", "f2": "bob", "f3": "chuck"}
  call write_map {input: file_to_name = map}
  call read_map
  output {
     Map[String, Int] out = read_map.out_map
     String contents = write_map.contents
     }
}
