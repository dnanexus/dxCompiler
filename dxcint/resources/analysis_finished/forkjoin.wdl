##
# Checks a simple branch and join operation.
# We start with a task, branch into two parallel executions, and then rejoin to calculate the result.
##

task mkFile {
  command {
    for i in `seq 1 1000`
    do
      echo $i
    done
  }
  output {
    File numbers = stdout()
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

task grep {
  String pattern = "10"
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

task join {
  Int grepCount
  Int wcCount
  command {
    expr ${wcCount} / ${grepCount}
  }
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

workflow forkjoin {
  call mkFile
  call grep { input: in_file = mkFile.numbers }
  call wc { input: in_file=mkFile.numbers }
  call join { input: wcCount = wc.count, grepCount = grep.count }
  output {
    Int join_out = join.proportion
  }
}
