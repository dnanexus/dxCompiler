version 1.0

task string_int_concat {
  Int disk_size = 100

  command <<<
    echo "hello"
  >>>

  runtime {
    disks: "local-disk " + disk_size + " HDD"
  }
}