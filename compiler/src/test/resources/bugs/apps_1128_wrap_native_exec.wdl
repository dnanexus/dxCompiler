version 1.0

workflow apps_1128_native_instance {

    Boolean n = true
    if (n) {
        call apps_1128_override_instance_name as conditional_instance_type {
            input:
                mock_input = "conditional_instance_type"
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
        dx_instance_type: "mem1_ssd1_v2_x8"
        dx_app: object {
                    type: "applet",
                    id: "applet-G9j8B680yzZVZ8X77XQ6ypXJ",
                    name: "apps_1177"
                }
    }
}
