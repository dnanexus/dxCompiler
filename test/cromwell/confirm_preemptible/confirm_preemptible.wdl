task check_preemption {
  Int set_preemptible

  command {
    curl "http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible" -H "Metadata-Flavor: Google"
  }

  output {
    String out = read_string(stdout())
  }

  runtime {
    preemptible: "${set_preemptible}"
    # includes curl
    docker: "dx://file-G66qz3Q0yzZfy6pg5q3yK3Kz"
  }
}


workflow confirm_preemptible {
  call check_preemption as yes { input: set_preemptible = 5 }
  call check_preemption as no  { input: set_preemptible = 0 }

  output {
    String should_be_true  = yes.out
    String should_be_false = no.out
  }
}
