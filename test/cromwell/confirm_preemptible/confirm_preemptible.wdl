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
    docker: "dx://file-GJ7q5KQ0yzZv7GQzJjX2x9Pb"
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
