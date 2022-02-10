version 1.1

task global_wf_from_wdl_multiply {
  input {
    Int multiply_first
    Int multiply_second
  }

  command <<< >>>

  output {
    Int product = 0
  }

  # Applet source code test/multi_user/global_wf_from_wdl_multiply/
  runtime {
    dx_app: object {
      type: "applet",
      id: "applet-G6ggFX80yzZykfx7JVP20BGK",
      name: "global_wf_from_wdl_multiply"
    }
  }
}

workflow global_wf_from_wdl_sub {

  # Native applet should be cloned, runnable
  call global_wf_from_wdl_multiply {
    input:
      multiply_first = 2,
      multiply_second = 4
  }
}
