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
When passing a large number of files (50+) from task to task, we recommend doing that in the form of a `tar.gz` archive to
improve the speed of (de)localizing those files between cloud storage and VMs that execute the WDL tasks. This will also make it easier 
to update input and output specifications of a running job. Of note, this suggestion is applicable for the scenarios when 
outputs of a task A are mapped 1:1 to the inputs of a task B, but not applicable when the outputs of a task A are 
to be used for a parallel execution, i.e. a scatter.
### Example
```wdl
task task_a {
...
    output {
        Array[File] files_out = [...]
    }
}
        
task task_b {
    input {
        Array[File] files_in
    }
...
}
        
task task_c {
    input {
        File single_file
    }
...
}

workflow wf_tar {
...
    call task_a {...}
    # APPLICABLE BELOW: use tar to optimize
    call task_b {input: files_in=task_a.files_out}
    # NOT APPLICABLE BELOW: need to preserve an array of output files, not to reduce them to a single archive
    scatter (s in  task_a.files_out) {
        call task_c {
          input:
          single_file = s
        }
    }
...
```
Since the `File` type declaration is represented by a string value, please follow the official [WDL documentation](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#strings) 
regarding the format and allowed characters.  
### Examples
**Scenario 1.** Task A creates `N` files to be used only by Task B, and there is no need for a separate archiving task.
In this case the archiving can be done in the Task A `task` section, and un-archiving should be done in the `task` 
section of the Task B.
```wdl
task task_a {
...
    command <<<
        # code that creates N files
        tar czf out.tgz ${DIRECTORY_WITH_N_FILES} 
    >>>
    output {
        File tarball = "out.tgz"
    }
}


task task_b {
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
    call task_a {...}
    call task_b {input: input_tarball=task_a.tarball}
...
}
```
**Scenario 2.** Task A creates `N` files to be used by multiple tasks B, C, ...M, here represented by 
`task_b_c_through_m`, it makes sense to refactor the archiving
```wdl
task task_a {
    ...
    command <<<
        # Creates N files e.g. with the pattern `file*.out`
    >>>
      output {
         Array[File] generated_files = glob("file*.out")
      }
}


task task_b_c_through_m {
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
    call task_a {...}
    call archive {input: task_a.generated_files}
    call task_b_c_through_m {input: input_tarball=archive.tarball}
...
}
```

**Scenario 3.** Task A creates 1..N files and is wrapped in a scatter with cardinality M. Then, all tarballs are to be combined. 
```wdl
task task_a {
    ...
    command <<<
        # Creates N files e.g. with the pattern `file*.out`
    >>>
      output {
         Array[File] generated_files = glob("file*.out")
      }
}
                
task archive {
    input {
        Array[File] files_to_tar
        Int archive_number = 0
    }
    command <<<
        archive_dir="archive_~{archive_number}"
        mkdir ${archive_dir}
        mv ~{sep=' ' files_to_tar} ${archive_dir}
        tar czf out_~{archive_number}.tgz ${archive_dir} 
    >>>
    output {
        File tarball = "out_~{archive_number}.tgz"
    }
}
        
task combine {
    input {
        Array[File] archives_to_combine
    }
    command <<<
        cat ~{sep=' ' archives_to_combine} >> out.tgz
    >>>
    output {
        File final_tarball = "out.tgz"
    }
}

workflow wf_tar {
...
    scatter (s in  range(M)) {
        call task_a {...}
        call archive {input: files_to_tar = task_a.generated_files, archive_number = s}
    }
    call combine {input: archives_to_combine = archive.tarball}
...
}
```
