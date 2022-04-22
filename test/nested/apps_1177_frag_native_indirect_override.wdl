version 1.0

workflow override {
  call create_input

  # Direct Call
  call mem_int {
    input: files = create_input.inputs
  }

  # Conditional, single call
  Boolean n = true
  if (n) {
    call mem_int as conditional_mem {
      input: files = create_input.inputs
    }
  }

  # Conditional, multiple calls
  Boolean m = true
  if (m) {
    call mem_int as conditional_mem_multi {
      input: files = create_input.inputs
    }
    call cpu_int as conditional_cpu_multi {
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
  runtime {
    memory: "30 GB"
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

task cpu_int {
  input {
    Array[File]+ files
    String? output_filename
  }

  command <<< >>>

  output {
    File file = "placeholder.txt"
  }

  runtime {
    cpu: 8
    dx_app: object {
              type: "app",
              id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
              name: "file_concatenator/2.0.0"
            }
  }
}