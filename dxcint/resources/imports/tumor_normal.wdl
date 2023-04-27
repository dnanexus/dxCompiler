version 1.1

workflow wf {
  input {
    String id
    File file
  }

  call varcall {
    input:
      id = id,
      file = file
  }

  output {
    File out = varcall.out
  }
}

task varcall {
  input {
    String id
    File file
  }

  command <<<
  head ~{file} > ~{id}.vcf
  >>>

  output {
    File out = "${id}.vcf"
  }
}