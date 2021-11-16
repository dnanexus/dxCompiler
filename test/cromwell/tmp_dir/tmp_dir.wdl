workflow tmp_dir {
    call mkTmpFile
    call writeToTmpDir

  output {
      String out = mkTmpFile.out
      String out2 = writeToTmpDir.tmpDir
  }
}

task mkTmpFile {
    command {
        echo "tmp_dir test wdl" > tmp
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
        String out = read_string("tmp")
    }
}

task writeToTmpDir {
   command {
        echo "tmp_dir test wdl 2" > $TMPDIR/tmp
        cat $TMPDIR/tmp
   }
   runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
   }
   output {
        String tmpDir = read_string(stdout())
   }
}
