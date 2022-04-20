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

  # In-task instance eval (with static instance name)
  Boolean l = true
  if (l) {
    call instance_type_dynamic as dynamic_instance {
      input: files = create_input.inputs
    }
  }

  # In-task instance eval (with RAM)
  Boolean x = true
  if (x) {
    call mem_int as memory_static {
      input: files = create_input.inputs
    }
  }

  # In-task instance eval (with RAM calculation)
  Boolean y = true
  if (y) {
    call mem_calc as memory_calculated {
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
    dx_instance_type: "mem1_ssd1_v2_x16"
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
    dx_instance_type: "mem1_ssd1_v2_x8"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}

task instance_type_dynamic {
  input {
    Array[File]+ files
    String? output_filename
  }
  Int ver = max(8,16)
  command <<< >>>

  output {
    File file = "placeholder.txt"
  }

  runtime {
    dx_instance_type: "mem1_ssd1_v2_x~{ver}"
    dx_app: object {
              type: "app",
              id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
              name: "file_concatenator/2.0.0"
            }
  }
}

task mem_int {
  input {
    Array[File]+ files
    String? output_filename
  }

  command <<< >>>

  output {
    File file = "placeholder.txt"
  }

  runtime {
    memory: "30 GB"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}
task mem_calc {
  input {
    Array[File]+ files
    String? output_filename
  }
  Int mem = max(30,15)

  command <<< >>>

  output {
    File file = "placeholder.txt"
  }

  runtime {
    memory: "~{mem} GB"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}