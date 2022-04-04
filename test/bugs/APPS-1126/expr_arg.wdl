version 1.0
import "inner.wdl" as inner

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

workflow expr_arg {
input {
    String files
}

        call file_outputs_from_input as fo {
            input: name = "~{files + files}"
        }
        call new_file as nf {
            input: name = select_first([fo.out])
        }
        call new_file as nf2 {
            input: name = select_first([nf.out])
        }
output {
    File? foo = nf2.out
}
}
