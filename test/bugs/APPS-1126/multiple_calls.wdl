task file_outputs_from_input {
    String outputName1
    String outputName2
    String outputName3

    command {
        echo "foo" > ${outputName1}
        echo "bar" > ${outputName2}
        echo "baz" > ${outputName3}.txt
    }
    output {
        File foo = outputName1
        File bar = "${outputName2}"
        File baz = "${outputName3}.txt"
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
}

task idempotent {
    File name1
    File name2
    File name3

    command {
        echo "running idempotent..."
    }
    output {
        File foo = name1
        File bar = name2
        File baz = name3
    }
}

task new_file {
    File name

    command {
        echo "something new" > filename.out
    }

    output {
        File out = "filename.out"
    }
}

task new_file2 {
    File name

    command {
        echo "something new" > filename.out
    }

    output {
        File out = "filename.txt"
    }
}

workflow file_outputs_from_input_wf2 {
    call file_outputs_from_input { input: outputName1 = "___something___", outputName2 = "___anything___", outputName3 = "___nothing___" }
    call idempotent as i1 { input: name1 = file_outputs_from_input.foo, name2 = file_outputs_from_input.bar, name3 = file_outputs_from_input.baz }
    call idempotent as i2 { input: name1 = i1.foo, name2 = i1.bar, name3 = i1.baz }
    call idempotent as i3 { input: name1 = i2.foo, name2 = i2.bar, name3 = i2.baz }
    call new_file as n1 {input: name = i3.foo}
    call new_file as n1a {input: name = n1.out}
    call new_file as n1b {input: name = n1a.out}
    call new_file2 as n1c {input: name = n1b.out}
    call new_file2 as n1d {input: name = n1c.out}
    call new_file2 as n1e {input: name = n1d.out}

    call new_file as n2 {input: name = i3.bar}
    call new_file as n2a {input: name = n2.out}

    output {
        File foo = n1e.out
        File foo2 = n2a.out
        String bar = read_string(i3.bar)
        String baz = read_string(i3.baz)
        File out = i3.foo
    }
}
