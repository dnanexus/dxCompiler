version 1.0

workflow apps_1197_native_frag_default_unit {

    # Memory
    call apps_1197_default_isntance {}

    # Conditional mem
    Boolean a = true
    if (a) {
        call apps_1197_default_isntance as conditional_mem {}
    }


}


task apps_1197_default_isntance {
  input {
    String? mock_input = "asdl"
  }

  command <<< >>>

  output {
    String mock_out = ""
  }

  runtime {
    dx_app: object {
                type: "applet",
                id: "applet-G9j8B680yzZVZ8X77XQ6ypXJ",
                name: "apps_1177"
            }
  }
}