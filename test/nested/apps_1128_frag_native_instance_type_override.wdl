version 1.1

workflow apps_1128_frag_native_instance_type_override {

  call apps_1128_default_instance {
    input:
      mock_input = "apps_1128_default_instance"
  }

  call apps_1128_override_instance_name {
    input:
      mock_input = "apps_1128_override_instance_name"
  }

  # Conditional, single call
  Boolean n = true
  if (n) {
    call apps_1128_override_instance_name as conditional_instance_type {
      input:
        mock_input = "conditional_instance_type"
    }
  }

  # scatter, single call
  scatter (i in [1,2]) {
    call apps_1128_override_instance_name as scatter_instance_type {
      input:
        mock_input = "scatter_instance_type~{i}"
    }
  }

  # Conditional, multiple calls
  Boolean m = true
  if (m) {
    call apps_1128_override_instance_name as conditional_instance_type_1 {
      input:
        mock_input = "conditional_instance_type_block1"
    }
    call apps_1128_override_instance_name_2 as conditional_instance_type_2 {
      input:
        mock_input = "conditional_instance_type_block2"
    }
  }

  output {
    String outs_apps_1128_default_instance = apps_1128_default_instance.instance_name
    String outs_apps_1128_override_instance_name = apps_1128_override_instance_name.instance_name
    String? outs_conditional_instance_type = conditional_instance_type.instance_name
    Array[String]? outs_scatter_instance_type = scatter_instance_type.instance_name
    String? outs_conditional_instance_type_1 = conditional_instance_type_1.instance_name
    String? outs_conditional_instance_type_2 = conditional_instance_type_2.instance_name
  }
}

task apps_1128_default_instance {
  input {
    String? mock_input
  }

  command <<< >>>

  output {
    String instance_name = ""
  }

  runtime {
    dx_app: object {
              type: "applet",
              id: "applet-G9j8B680yzZVZ8X77XQ6ypXJ",
              name: "apps_1177"
            }
  }
}
task apps_1128_override_instance_name {
  input {
    String? mock_input
  }

  command <<< >>>

  output {
    String instance_name = ""
  }

  runtime {
    dx_instance_type: "mem1_ssd1_v2_x2"
    dx_app: object {
              type: "applet",
              id: "applet-G9j8B680yzZVZ8X77XQ6ypXJ",
              name: "apps_1177"
            }
  }
}

task apps_1128_override_instance_name_2 {
  input {
    String? mock_input
  }

  command <<< >>>

  output {
    String instance_name = ""
  }

  runtime {
    dx_instance_type: "mem1_ssd1_v2_x4"
    dx_app: object {
              type: "applet",
              id: "applet-G9j8B680yzZVZ8X77XQ6ypXJ",
              name: "apps_1177"
            }
  }
}
