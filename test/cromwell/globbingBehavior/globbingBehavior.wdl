version 1.1
task writeFiles {
  input {
    String dir
  }
  command <<<
    mkdir ${dir}
    echo "hullo" > ${dir}/hello.txt

    echo "what" > what.txt
    echo "how" > how.txt
    echo "why" > why.txt

    cd ${dir}
    echo "buh-bye" > ciao.txt
  >>>
  output {
    Array[File] strInterpolation =  glob("${dir}/*.txt")
    Array[File] callDirectory = glob("*.txt")
    Array[File] specificDirectory = glob(dir + "/*.txt")
  }
  runtime {
    docker:"ubuntu:latest"
  }
}

task unexpectedReturnCode {

  command {
    echo "easiest glob job ever!" > tiny.txt
    exit 3
  }

  output {
    Array[File] easyGlob = glob("*.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    returnCodes: 3
  }
}

task writeLinks {
  input {
    String dir
  }
  command <<<
    mkdir ${dir}
    cd ${dir}

    echo "hi file" > orig_file.txt

    ln orig_file.txt hard_link.txt

    ln -s hard_link.txt sym_link.txt
    ln -s sym_link.txt sym_link_2.txt

    ln sym_link_2.txt hard_link_to_sym.txt

    echo "updated file" > orig_file.txt
  >>>
  output {
    Array[File] origFileContents =  glob("${dir}/orig_file*.txt")
    Array[File] symLinkContents =  glob("${dir}/sym_link*.txt")
    Array[File] hardLinkContents =  glob("${dir}/hard_link*.txt")
  }
  runtime {
    docker:"ubuntu:latest"
  }
}

workflow globbingBehavior {
  input {
    String outputDir = "out"
  }

  call writeFiles { input: dir = outputDir }

  call unexpectedReturnCode {}

  call writeLinks { input: dir = outputDir }

  scatter (f in writeFiles.strInterpolation) {
    String strInterpolation_contents = read_string(f)
  }

  scatter (f in writeFiles.callDirectory) {
    String callDirectory_contents = read_string(f)
  }

  scatter (f in writeFiles.specificDirectory) {
    String specificDirectory_contents = read_string(f)
  }

  scatter (f in unexpectedReturnCode.easyGlob) {
    String easyGlob_contents = read_string(f)
  }

  scatter (f in writeLinks.origFileContents) {
    String origFileContents_contents = read_string(f)
  }

  scatter (f in writeLinks.symLinkContents) {
    String symLinkContents_contents = read_string(f)
  }

  scatter (f in writeLinks.hardLinkContents) {
    String hardLinkContents_contents = read_string(f)
  }

  output {
    Array[String] strInterpolation_contents_out = strInterpolation_contents
    Array[String] callDirectory_contents_out = callDirectory_contents
    Array[String] specificDirectory_contents_out = specificDirectory_contents
    Array[String] easyGlob_contents_out = easyGlob_contents
    Array[String] origFileContents_contents_out = origFileContents_contents
    Array[String] symLinkContents_contents_out = symLinkContents_contents
    Array[String] hardLinkContents_contents_out = hardLinkContents_contents
 }
}
