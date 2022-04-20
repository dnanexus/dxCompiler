version 1.1

workflow override {

    Boolean a = true
    if (a) {
      call instance_type
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

}