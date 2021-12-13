version 1.0

task add {
  input {
    Int a
    Int c
  }
  command {
    echo $((${a} + ${c}))
  }
  meta {
    summary: "Adds two int together"
  }

  output {
    Int result = read_int(stdout())
  }
}
