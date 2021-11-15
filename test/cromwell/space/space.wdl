task hello {
    command {
        echo "Hello world" > 'out with space.txt'
        echo "Hello world" > 'out%20with%20.uxu'
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
        File singleSpace = "out with space.txt"
        File singlePercent = "out%20with%20.uxu"
        Array[File] globbedSpace = glob("*.txt")
        Array[File] globbedPercent = glob("*.uxu")
    }
}

task goodbye {
    File f
    File f2
    Array[File] files
    Array[File] files2
    command {
        cat "${f}"
        cat "${f2}"
        cat "${sep=" " files}"
        cat "${sep=" " files2}"
    }
    runtime {
            docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
        Array[String] out = read_lines(stdout())
    }
}

workflow space {
    call hello
    String s1 = read_string(hello.singleSpace)
    String s2 = read_string(hello.singlePercent)
    
    call goodbye { input: 
     f = hello.singleSpace,
     f2 = hello.singlePercent, 
     files = hello.globbedSpace,
     files2 = hello.globbedPercent 
    }
    
    output {
        String o1 = s1
        String o2 = s2
        Array[String] o3 = goodbye.out
    }
}
