# When a task is dependent on a failing task, it should not start. When a task is not dependent
# on a failing task, it should run and succeed.

task boundToFail {
  command {
    cd fake/file/path
  }
  output {
   String badOutput = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

task shouldNotStart {
  String str
    command {
     echo ${str}
    }
    runtime {
       docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
}

task shouldSucceed {
  String str
    command {
     sleep 100
     echo ${str}
    }
    runtime {
       docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
     String stalling = read_string(stdout())
    }
}

task delayedTask {
  String str_2
    command {
     echo ${str_2}
    }
    runtime {
       docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
     String notUsed = read_string(stdout())
    }
}

workflow cont_while_possible {
  call boundToFail
  call shouldNotStart { input: str = boundToFail.badOutput }
  call shouldSucceed { input: str = "echoing nonsense" }
  call delayedTask  { input: str_2 = shouldSucceed.stalling }
}
