import "hello_world.wdl" as hello
import "fjoin.wdl"
import "sub_interactions.wdl" as interactions

task join {
  Int grepCount
  Int wcCount
  command {
    expr ${wcCount} / ${grepCount}
  }
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow recursive_imports {
  call hello.sub.wf_hello {}
  call fjoin.join { input: wcCount=15, grepCount=10 }
  call interactions.counter.countTo { input: value = 3 }

  output {
    Int join_out = join.proportion
  }
}
