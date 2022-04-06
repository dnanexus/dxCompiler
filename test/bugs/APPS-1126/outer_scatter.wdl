version 1.0

task nf {
    input {
        File name
    }
    command {
        echo "something new" > filename2.out
    }

    output {
        File out = "filename2.out"
    }
}

workflow outer_scatter {
    input {
        Array[String] files
    }


    output {
        Array[File] foo = nf.out
    }
}
