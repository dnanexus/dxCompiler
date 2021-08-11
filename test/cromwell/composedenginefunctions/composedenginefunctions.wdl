##
# Check that we can:
# - Compose engine functions together.
##

task composeEngineFunctions {

  command {
    echo "Hello, I am a small test string"
    echo 2 >&2
  }
  output {
    File blah = stdout()
    String x = read_string(stdout())
    String y = read_int(stderr()) + x + read_string(blah)
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow composedenginefunctions {
  call composeEngineFunctions

  output {
    String x_out = composeEngineFunctions.x
    String y_out = composeEngineFunctions.y
  }
}
