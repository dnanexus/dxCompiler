version 1.0

task file_outputs_from_input {
    input {
        String name
    }

    command {
        echo "foo" > ${name}
    }
    output {
        File out = name
    }
}

task new_file {
    input {
        File name
    }
    command {
        echo "something new" > filename.out
    }

    output {
        File out = "filename.out"
    }
}

workflow conditional {
input {
    Array[String] files
    Boolean makeCall
}

    if (makeCall) {
        call file_outputs_from_input as fo {
            input: name = select_first(files)
        }
        call new_file as nf {
            input: name = fo.out
        }
        call new_file as nf2 {
            input: name = nf.out
        }
    }
output {
    File? foo = nf2.out
}
}
