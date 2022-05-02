version 1.0

workflow apps_1177_native_frag_indirect_override_unit {

    # Memory
    call apps_1177_mem_int {}

    # Conditional mem
    Boolean a = true
    if (a) {
        call apps_1177_mem_int as conditional_mem {}
    }


}


task apps_1177_mem_int {
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
                type: "applet",
                id: "applet-G9VZBF00yzZgj97K42BQKJ7Z",
                name: "apps_1177"
            }
  }
}