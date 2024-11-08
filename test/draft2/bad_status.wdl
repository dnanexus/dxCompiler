task BadCommand {
      command {
          ls /xx/yyy
      }
      runtime {
          docker: "broadinstitute/genomes-in-the-cloud:2.2.5-1485277291"
      }
      output {
          Int rc = 1
      }
}

workflow bad_status {
    call BadCommand
    output {
        Int rc = BadCommand.rc
    }
}
