task printInt {
  Int? int

  command { echo "${int}" > out.txt }
  output { Int out = read_int("out.txt") }
  runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

workflow arrays_scatters_ifs {

  Array[Array[Int]] table = [[0,0], [1,1]]
  scatter (row in table) {

    if (length(row) == 2) {
      Int int = row[1]
    }

    call printInt {input: int=int }
  }
  output {
    Array[Int] ints = printInt.out
  }
}
