# Tips and tricks

This page contains examples of coding patterns in WDL. Sometimes, the correct WDL way of
doing things is not immediately obvious, and here we provide some worked examples.

## Optional argument with a default

Examine a situation where you want to write a task that takes an
optional argument with a default. For example, the task `hello` takes
an optional username argument. If this is not supplied, the default is `dave.jones`.

```wdl
workflow foo {
    String? username
    hello { input: who = username }
}

task hello {
    String? who
    String who_actual = select_first([who, "dave.jones"])

    command {
        echo "Hello, ${who_actual}"
    }
    output {
        String result = read_string(stdout())
    }
}
```

## Assertions

Currently, WDL does not have a way to throw exceptions, or report errors. A way around
this limitation is to write an `assert` task (suggested by [@mlin](https://github.com/mlin)).

```wdl
task assert {
    Boolean value
    String msg

    command {
        if [ "${value}" == "false" ]; then
            echo $msg
            exit 1
        fi
        exit 0
    }
}
```

You can call it from a workflow when you want to check a condition.

```wdl
workflow foo {
    call assert { input: value= CONDITION_TO_CHECK,
                         msg = ERROR_MESSAGE }
}
```


A method that works directly with the platform error reporting mechanism, suggested by [@jtratner](https://github.com/jtratner), is:

```wdl
workflow show_error {
    call make_error {}
}

task make_error {
    command <<<
        echo '{"error": {"type": "AppError", "message": "x must be at least 2' $(hostname)'"}}' > job_error.json
        exit 1
    >>>
}
```

(including the hostname in the error message allows the message to have the root
job ID, facilitating debugging in complex workflows that may be deeply nested.)


## Large number of input/output files
When passing a large number of files from stage to stage, we recommend doing that in the form of a `tar.gz` archive to 
improve the I/O of the upload and download transactions with the DNAnexus platform storage. This will also make it easier 
to update input and output specifications of a running job. Of note, this suggestion is applicable for the scenarios when 
outputs of a stage A are mapped 1:1 to the inputs of a stage B, but not applicable when the outputs of a stage A are 
to be used for a parallel execution, i.e. a scatter.  
Since the `File` type declaration is represented by a string value, please follow the official [WDL documentation](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#strings) 
regarding the format and allowed characters.  
#### Examples
**Scenario 1.** Stage A creates `N` files to be used only by Stage B, and there is no need for a separate archiving task.
In this case the archiving can be done in the Stage A `task` section, and un-archiving should be done in the `task` 
section of the Stage B.
```wdl
task stage_a {
...
    command <<<
        # code that creates N files
        tar czf out.tgz ${DIRECTORY_WITH_N_FILES} 
    >>>
    output {
        File tarball = "out.tgz"
    }
}


task stage_b {
    input {
        File input_tarball
    }
    command <<<
        mkdir -p unpack
        tar xzvf ~{input_tarball} -C unpack 1>tarball_listing
        # continue your analysis here
    >>>
...
}

workflow wf_tar {
...
    call stage_a {...}
    call stage_b {input: input_tarball=stage_a.tarball}
...
}
```
**Scenario 2.** Stage A creates `N` files to be used by Stage B, C, ... m, and it makes sense to refactor the archiving 
step into a stand-alone stage.
```wdl
task stage_a {
    ...
    command <<<
        # Creates N files e.g. with the pattern `file*.out`
    >>>
      output {
         Array[File] generated_files = glob("file*.out")
      }
}


task stage_b_c_through_m {
    input {
        File input_tarball
    }
    command <<<
        mkdir -p unpack
        tar xzvf ~{input_tarball} -C unpack 1>tarball_listing
        # continue your analysis here
    >>>
    ...
}
        
task archive {
    input {
        Array[File] files_to_tar
    }
    command <<<
        tar czf out.tgz -C / ~{sep=' ' files_to_tar} 
    >>>
    output {
        File tarball = "out.tgz"
    }
}

workflow wf_tar {
...
    call stage_a {...}
    call archive {input: stage_a.generated_files}
    call stage_b_c_through_m {input: input_tarball=archive.tarball}
...
}
```