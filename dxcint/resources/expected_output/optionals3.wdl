version 1.0

workflow optionals3 {
  input {
    Int? number_to_scatter
  }

  Int defined_number_to_scatter = select_first([number_to_scatter, 2])

  scatter (i in range(defined_number_to_scatter)) {
    Int j = i + 1
  }

  output {
    Array[Int] k = j
  }
}