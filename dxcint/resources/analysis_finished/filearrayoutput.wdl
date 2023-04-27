task createFileArray {
  command <<<
    mkdir out
    echo "hullo" > out/hello.txt
    echo "buh-bye" > out/ciao.txt
  >>>
  output {
    Array[File] out = [ "out/hello.txt", "out/ciao.txt" ]
  }
  runtime {docker:"ubuntu:latest"}
}

task combiner {
  Array[File] in_file
  command <<<
    cat ${sep=' ' in_file}
  >>>
  output {
    String result = read_string(stdout())
  }
  runtime {docker:"ubuntu:latest"}
}

workflow filearrayoutput {
    call createFileArray
    call combiner { input: in_file = createFileArray.out }
    output {
        String combiner_out = combiner.result
    }
}
