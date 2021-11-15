version 1.0

task make_files {
  command {
    for f in f1 f2 f3; do touch $f; done
  }
  output {
    Array[File] files = ["f1", "f2", "f3"]
  }
  runtime {
    docker: "dx://file-G66qz3Q0yzZfy6pg5q3yK3Kz"
  }
}

task write_map {
  input {
    Map[File, String] file_to_name
  }
  command {
    cat ${write_map(file_to_name)}
  }
  output {
    String contents = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qz3Q0yzZfy6pg5q3yK3Kz"
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
    docker: "dx://file-G66qz3Q0yzZfy6pg5q3yK3Kz"
  }
}
workflow wf {
  call make_files
  Map[File, String] map = {make_files.files[0]: "alice", make_files.files[1]: "bob", make_files.files[2]: "chuck"}
  call write_map {input: file_to_name = map}
  call read_map
}
