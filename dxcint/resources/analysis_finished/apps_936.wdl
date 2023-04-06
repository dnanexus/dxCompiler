version 1.0

struct Foo {
  String docker_image
}

task apps_936 {
  input {
    Foo foo
  }

  command <<<
  echo 'hello'
  >>>

  output {
    String docker_image = foo.docker_image
  }

  runtime {
    docker: foo.docker_image
  }
}