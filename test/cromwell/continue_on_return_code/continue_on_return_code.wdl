version 1.1

task continueOnRC1 {
    command <<<
        echo "echo 'OH NO!'" > script.sh
        echo "exit 1" >> script.sh
        chmod +x script.sh
        ./script.sh
    >>>
    output {
        File ohno = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
        returnCodes: [1]
    }
}

task continueOnRC2 {
    command <<<
        echo "echo 'OH NO!'" > script.sh
        echo "exit 12" >> script.sh
        chmod +x script.sh
        ./script.sh
    >>>
    output {
        String ohno = read_string(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
        returnCodes: [12]
    }
}

task continueOnRC3 {
    command <<<
        echo "echo 'OH NO!'" > script.sh
        echo "exit 123" >> script.sh
        chmod +x script.sh
        ./script.sh
    >>>
    output {
        String ohno = read_string(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
        returnCodes: [1, 12, 123, 23, 3]
    }
}

task finisher {
    input {
        File in1
        String in2
        String in3
    }
    command <<<
       cat ${in1} && echo ${in2} && wc -l <<< "${in3}"
    >>>
    output {
        String finished = read_string(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
        returnCodes: [0]
    }
}

workflow runtime_continueOnRC {
    call continueOnRC1
    call continueOnRC2
    call continueOnRC3
    call finisher { input:
        in1=continueOnRC1.ohno,
        in2=continueOnRC2.ohno,
        in3=continueOnRC3.ohno
    }

    output {
        String finisher_out = finisher.finished
    }
}
