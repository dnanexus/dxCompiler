task make_files {
  command {
    echo "content" > file1
    echo "content" > file2
  }

  output {
    Array[File] files = ["file1", "file2"]
  }
}

task localize {
    Array[File] array
    command {
        ls -1 "$(dirname ${array[0]})" | wc -l | tr -d '[[:space:]]'
    }

    output {
        String ls = read_string(stdout())
    }
}

task echo_int {
  Int int
  command {
    echo ${int} > out
  }
  output {File out = "out"}
}

task localize_with_docker {
    Array[File] array
    command {
        ls -1 "$(dirname ${array[0]})" | wc -l | tr -d '[[:space:]]'
    }
    output {
        String ls = read_string(stdout())
    }
    runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

task echo_int_with_docker {
  Int int
  command {
    echo ${int} > out
  }
  output {File out = "out"}
  runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

workflow wf {
    Array[Int] ints = [1,2]

   call make_files

   scatter(i in ints) {
    call echo_int {
      input: int = i
    }
   }

   scatter(i in ints) {
      call echo_int_with_docker {
        input: int = i
      }
   }

  call localize as fromDifferentDirectories { input: array = echo_int.out }
  call localize as fromSameDirectory { input: array = make_files.files }

  call localize_with_docker as fromDifferentDirectories_with_docker { input: array = echo_int_with_docker.out }
  call localize_with_docker as fromSameDirectory_with_docker { input: array = make_files.files }
}
