version 1.0

workflow frag_default {

    Boolean a = true
    if (a) {
      call instance_type
    }
}

task instance_type {
  command <<< >>>
  runtime {
    dx_instance_type: "mem1_ssd1_v2_x16"
  }
}
