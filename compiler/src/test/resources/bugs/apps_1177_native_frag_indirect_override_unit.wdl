version 1.0

workflow override {

    # Memory
    call mem_int {}

    # Conditional mem
    Boolean a = true
    if (a) {
        call mem_int as conditional_mem {}
    }


}


task mem_int {
  input {
    String? mock_input = "asdl"
  }

  command <<< >>>

  output {
    String mock_out = ""
  }

  runtime {
    memory: "30 GB"
    dx_app: object {
              type: "app",
              id: "app-G6G0jX80g1FZX1Z57z3zbG6v",
              name: "mock_app_sciprodbuild/0.0.6"
            }
  }
}