version 1.0

workflow apps_1197_native_frag_default {

  call apps_1197_default_instance {
    input:
      mock_input = "apps_1197_default_instance"
  }

  # default in frag
  Boolean a = true
  if (a) {
    call apps_1197_default_instance as frag_default {
      input:
        mock_input = "frag_default"
    }
  }

  # meta section
  call apps_1197_default_instance_meta {
    input:
      mock_input = "apps_1197_default_instance_meta"
  }


  # default in frag with meta
  Boolean b = true
  if (a) {
    call apps_1197_default_instance_meta as frag_default_meta {
      input:
        mock_input = "frag_default_meta"
    }
  }

  output {
    String outs_apps_1197_default_instance = apps_1197_default_instance.instance_name
    String? outs_frag_default = frag_default.instance_name
    String outs_apps_1197_default_instance_meta = apps_1197_default_instance_meta.instance_name
    String? outs_frag_default_meta = frag_default_meta.instance_name
  }
}

task apps_1197_default_instance {
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


task apps_1197_default_instance_meta {
  input {
    String? mock_input
  }

  command <<< >>>

  output {
    String instance_name = ""
  }

  meta {
    type : "native"
    id : "applet-G9j8B680yzZVZ8X77XQ6ypXJ"
  }
}
