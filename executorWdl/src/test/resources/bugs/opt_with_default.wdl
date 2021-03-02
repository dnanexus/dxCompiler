version 1.1

workflow Test {
  input {
    Int? opt
    File file
  }

  Int ref_size = ceil(size(file, "GB"))

  call Foo {
    input:
      opt=opt,
      non_opt=ref_size
  }
}

task Foo {
  input {
    Int? opt = 20
    Int non_opt
  }
  command <<<
  echo ~{default="20" opt}
  >>>
}
