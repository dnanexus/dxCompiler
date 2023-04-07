version 1.0

workflow nested_optional {
  input {
    String? s
  }

  Array[String?] opt_strings = [s]

  # make sure `file` is not a nested optional
  if (defined(s)) {
    String? str = select_first(opt_strings)
  }

  call nested_opt_task {
    input: strings = select_all(opt_strings)
  }

  output {
    String? outstr = str
    String? outstr2 = nested_opt_task.sout
  }
}

task nested_opt_task {
  input {
    Array[String] strings
    String? s1
  }

  String? s2 = if length(strings) > 0 then strings[0] else s1

  command <<< >>>

  output {
    String? sout = s2
  }
}