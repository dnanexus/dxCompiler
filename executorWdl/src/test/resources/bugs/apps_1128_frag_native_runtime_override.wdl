version 1.1

workflow override {

    Boolean a = true
    if (a) {
      call instance_type
    }

    Boolean b = true
    if (b) {
      call mem_int
    }

	Boolean c = true
	if (c) {
      call mem_calc
	}
}



task instance_type {
  command <<< >>>
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
  command <<< >>>
  output {
    String hello = "Hello World With Memory"
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
  Int mem = max(30,15)

  command <<< >>>

  output {
    String hello = "Hello World With Memory Dynamic"
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