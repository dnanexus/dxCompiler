version 1.1

workflow override {
  call create_input

  call default {
    input: files = create_input.inputs
  }

  call instance_type {
    input: files = create_input.inputs
  }

  # Conditional, single call
  Boolean n = true
  if (n) {
    call instance_type as conditional_isntance_type {
      input: files = create_input.inputs
    }
  }

  # scatter, single call
  scatter (i in [1,2]) {
    call instance_type as scatter_isntance_type {
      input: files = create_input.inputs
    }
  }

  # Conditional, multiple calls
  Boolean m = true
  if (m) {
    call instance_type as conditional_isntance_type_1 {
      input: files = create_input.inputs
    }
    call instance_type_2 as conditional_isntance_type_2 {
      input: files = create_input.inputs
    }
  }
}

task create_input {
  command <<<
    echo "a" > input_a.txt
    echo "b" > input_b.txt
  >>>
  output {
    Array[File] inputs = ["input_a.txt", "input_b.txt"]
  }
}

task default {
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
task instance_type {
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

task instance_type_2 {
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
