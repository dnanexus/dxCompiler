version 1.1

workflow override {
    call create_input {}

    call default {
      input: files = create_input.inputs
    }

    call instance_type {
      input: files = create_input.inputs
    }

    Boolean n = true
    if (n) {
      call instance_type as conditional_isntance_type {
        input: files = create_input.inputs
      }
    }

    call mem_int {
      input: files = create_input.inputs
    }

    call mem_calc {
      input: files = create_input.inputs
    }

	Boolean t = true
	if (t) {
      call mem_calc as conditional_mem {
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
    memory: "${mem} GB"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}