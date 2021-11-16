task countTo {
    Int value
    command {
        seq 0 1 ${value}
    }
    runtime {
          docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
      }
    output {
        File range = stdout()
    }
}

task filterEvens {
    File numbers
    command {
        grep '[02468]$' ${numbers} > evens
    }
    runtime {
          docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
      }
    output {
        File evens = "evens"
    }
}

workflow countEvens {
    Int max = 10
    
    call countTo { input: value = max }
    call filterEvens { input: numbers = countTo.range }
    output {
        String someStringOutput = "I'm an output"
        File evenFile = filterEvens.evens
    }
}
