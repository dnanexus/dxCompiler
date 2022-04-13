version 1.1


task mem_calc {
  input {
    Array[Int] nums = [1, 2]
  }

  command <<< >>>

  output {
    String out = "Hello World"
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