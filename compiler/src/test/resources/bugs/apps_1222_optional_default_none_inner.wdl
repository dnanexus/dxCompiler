version 1.1

workflow apps_1222_optional_default_none_inner {
  input {
    Int arg1 = 1
    String? caller_to_emulate = None
    Int? maybe_none_int = None
    File? maybe_none_file = None
  }
  call none_inner_t1 {input: x = arg1}
  output {
    String outer_out = none_inner_t1.rg_out
  }
}

task none_inner_t1 {
  input {
    Int x
  }

  command {
    echo ~{x}
  }
  output {
    String rg_out = "~{x}"+"~{x}"
  }
}