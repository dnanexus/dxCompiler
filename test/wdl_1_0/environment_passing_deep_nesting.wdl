version 1.0

workflow environment_passing_deep_nesting {
    input {
        Boolean verify
    }

    call bear
    call bull

    if (verify) {
        call compare { input:
            a = bear.result,
            b = bull.result
        }

        call assert
    }

    output {
        Int? result = compare.result
    }
}


task bear {
    command {}
    output {
        Int result = 10
    }
}

task bull {
  command {}
  output {
      Int result = 3
  }
}

task assert {
    command {}
}


task compare {
    input {
        Int a
        Int b
    }
    command {}
    output {
        Int result = a - b
  }
}
