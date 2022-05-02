version 1.0

workflow apps_1177_native_indirect_override {

  call apps_1177_default_instance {
    input:
      mock_input = "apps_1177_default_instance"
  }

  call apps_1177_mem_int {
    input:
      mock_input = "apps_1177_mem_int"
  }

  # Conditional, single call
  Boolean n = true
  if (n) {
    call apps_1177_mem_int as conditional_mem_int {
      input:
        mock_input = "conditional_mem_int"
    }
  }

  # scatter, single call
  scatter (i in [1,2]) {
    call apps_1177_mem_int as scatter_mem_int {
      input:
        mock_input = "scatter_mem_int~{i}"
    }
  }

  # Conditional, multiple calls
  Boolean m = true
  if (m) {
    call apps_1177_mem_int as conditional_mem_int_1 {
      input:
        mock_input = "conditional_instance_type_block1"
    }
    call apps_1177_cpu_int as conditional_cpu_int_2 {
      input:
        mock_input = "conditional_instance_type_block2"
    }
  }

  output {
    String outs_apps_1177_default_instance = apps_1177_default_instance.instance_name
    String outs_apps_1177_mem_int = apps_1177_mem_int.instance_name
    String? outs_conditional_mem_int = conditional_mem_int.instance_name
    Array[String]? outs_scatter_mem_int = scatter_mem_int.instance_name
    String? outs_conditional_mem_int_1 = conditional_mem_int_1.instance_name
    String? outs_conditional_cpu_int_2 = conditional_cpu_int_2.instance_name
  }
}

task apps_1177_default_instance {
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
              id: "applet-G9VZBF00yzZgj97K42BQKJ7Z",
              name: "apps_1177"
            }
  }
}
task apps_1177_mem_int {
  input {
    String? mock_input
  }

  command <<< >>>

  output {
    String instance_name = ""
  }

  runtime {
    memory: "16 GB"
    dx_app: object {
              type: "applet",
              id: "applet-G9VZBF00yzZgj97K42BQKJ7Z",
              name: "apps_1177"
            }
  }
}

task apps_1177_cpu_int {
  input {
    String? mock_input
  }

  command <<< >>>

  output {
    String instance_name = ""
  }

  runtime {
    cpu: 4
    dx_app: object {
              type: "applet",
              id: "applet-G9VZBF00yzZgj97K42BQKJ7Z",
              name: "apps_1177"
            }
  }
}
