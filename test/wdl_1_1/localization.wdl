version 1.1

task t1 {
  input {
    Int label
  }
  Array[Int] a = [1,2,3]
  command <<<
  for i in ~{sep(' ', a)}; do
    echo "${i}_~{label}" > ${i}
  done
  >>>
  output {
    Array[File] out = glob("[0-9]")
  }
}

task t2 {
  input {
    Pair[File,File] p
  }
  command <<<
  cat ~{p.right} ~{p.left} > out
  >>>
  output {
    File out = "out"
  }
}

workflow zip {
  call t1 as c1 {
     input:
        label=1
  }
  call t1 as c2 {
     input:
        label=2
  }
  scatter (p in zip(c1.out, c2.out)) {
     call t2 as t2_scattered {
       input: p=p
    }
  }
  output {
    Array[File] af=t2_scattered.out
  }
}
