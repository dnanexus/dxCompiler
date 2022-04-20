version 1.0

workflow simple_wf {
    input {}
    call simple_task
    output {
      Array[File] out = simple_task.out
    }
}

task simple_task {

  input {}

  command <<<
    echo "out12" > out.vcf
    echo "out2" > out.bam
  >>>

  output {
     Array[File] out= ["out.vcf", "out.bam"]
  }
}