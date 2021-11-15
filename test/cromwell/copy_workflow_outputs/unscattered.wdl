
task A {
  command {
    echo "Enfin un peu de francais pour contrer ce raz-de-marée anglais !" > out
    echo "Jacques Chirac fait du jetski sur la Seine en costume traditionnel russe" > out2
  }
  output {
    File out = "out"
    File out2 = "out2"
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

task B {
  command {
     echo "Je contre avec un bonnet peruvien et tire une carte chance" > out
     echo "Kamoulox !" > out2
  }
  output {
     Array[File] outs = ["out", "out2"]
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

task C {
  command {
    cat > out <<END
    (_)
     !_________________________________________
     !*  *  *  *  * |##########################|
     ! *  *  *  *  *|                          |
     !*  *  *  *  * |##########################|
     ! *  *  *  *  *|                          |
     !*  *  *  *  * |##########################|
     ! *  *  *  *  *|                          |
     !*  *  *  *  * |##########################|
     !~~~~~~~~~~~~~~~                          |
     !#########################################|
     !                                         |
     !#########################################|
     !                                         |
     !###################################JGS###|
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !
     !
     !
     !
     !
     !
     !
    END
  }
  output {
      File out = "out"
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow wfoutputs {
  call A
  call B
  call C
  output {
    Array[File] a_out = A.*
    Array[File] b_out = B.outs
  }
}
