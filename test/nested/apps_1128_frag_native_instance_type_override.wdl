version 1.1

workflow apps_1128_frag_native_instance_type_override {
  call apps_1128_create_input

  call apps_1128_default_instance {
    input: files = apps_1128_create_input.inputs
  }

  call apps_1128_override_instance_name {
    input: files = apps_1128_create_input.inputs
  }

  # Conditional, single call
  Boolean n = true
  if (n) {
    call apps_1128_override_instance_name as conditional_isntance_type {
      input: files = apps_1128_create_input.inputs
    }
  }

  # scatter, single call
  scatter (i in [1,2]) {
    call apps_1128_override_instance_name as scatter_isntance_type {
      input: files = apps_1128_create_input.inputs
    }
  }

  # Conditional, multiple calls
  Boolean m = true
  if (m) {
    call apps_1128_override_instance_name as conditional_isntance_type_1 {
      input: files = apps_1128_create_input.inputs
    }
    call apps_1128_override_instance_name_2 as conditional_isntance_type_2 {
      input: files = apps_1128_create_input.inputs
    }
  }
}

task apps_1128_create_input {
  command <<<
    echo "a" > input_a.txt
    echo "b" > input_b.txt
  >>>
  output {
    Array[File] inputs = ["input_a.txt", "input_b.txt"]
  }
}

task apps_1128_default_instance {
  input {
    Array[File]+ files
    String? output_filename
  }

  command <<< >>>

  output {
    File file = "placeholder.txt"
  }

  runtime {
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}
task apps_1128_override_instance_name {
  input {
    Array[File]+ files
    String? output_filename
  }

  command <<< >>>

  output {
    File file = "placeholder.txt"
  }

  runtime {
    dx_instance_type: "mem1_ssd1_v2_x2"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}

task apps_1128_override_instance_name_2 {
  input {
    Array[File]+ files
    String? output_filename
  }

  command <<< >>>

  output {
    File file = "placeholder.txt"
  }

  runtime {
    dx_instance_type: "mem1_ssd1_v2_x4"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}
