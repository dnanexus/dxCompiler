version 1.0

workflow override {

  # Default
  call default {}

  # Memory
  call mem_int {}

  # CPU
  call cpu_int {}
}

task default {
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
                    id: "applet-G9VZBF00yzZgj97K42BQKJ7Z",
                    name: "apps_1177"
                }
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
                type: "applet",
                id: "applet-G9VZBF00yzZgj97K42BQKJ7Z",
                name: "apps_1177"
            }
  }
}

task cpu_int {
  input {
    String? mock_input = "asdl"
  }

  command <<< >>>

  output {
    String mock_out = ""
  }

  runtime {
    cpu: 8
    dx_app: object {
                type: "applet",
                id: "applet-G9VZBF00yzZgj97K42BQKJ7Z",
                name: "apps_1177"
            }
  }
}