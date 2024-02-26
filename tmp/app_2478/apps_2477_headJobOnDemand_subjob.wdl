version 1.0


task apps_2477_headJobOnDemand_subjob_t01 {
    input {
        File fruits
    }
    # provide dynamical size of the disk
    Int disk_req_gb = ceil(size(fruits, "GB")) + 50

  command <<<
    lines=$(df -t btrfs | grep dev)
    size_kb=$(echo $lines | cut -d ' ' -f 2)
    let "size_gb= $size_kb / (1024 * 1024)"
    if [[ $size_gb -ge disk_req_gb ]]; then
       echo "true"
    else
       echo "false"
   fi
 >>>

    runtime {
        disks: "local-disk ${disk_req_gb} HDD"
    }
    output {
        String retval = read_string(stdout())
    }
}

task apps_2477_headJobOnDemand_subjob_t02 {
    input {
        Int num = 1
    }
    # provide dynamical size of the disk

  command <<<
    echo 1
 >>>

    output {
        String t02_retval = "Hello"
    }
}
