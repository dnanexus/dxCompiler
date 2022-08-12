version 1.0
workflow apps_1269_1270_unqualified_ids_manifest {
    input {
        File missing_file_qualified
        File cloned_file_unqualified
        File cloned_file_qualified
    }

    Boolean k = true
    Boolean l = true
    Boolean m = true

    if (k) {
        call apps_1269_missing_file_qualified {
            input:
                task_missing_file = missing_file_qualified
        }
    }

    if (l) {
        call apps_1270_cloned_file_unqualified {
            input:
                task_cloned_file_unqualified = cloned_file_unqualified
        }
    }

    if (m) {
        call apps_1270_cloned_file_qualified {
            input:
                task_cloned_file_qualified = cloned_file_qualified
        }
    }


    output {
        Int? o1 = apps_1269_missing_file_qualified.out1
        Int? o2 = apps_1270_cloned_file_unqualified.out2
        Int? o3 = apps_1270_cloned_file_qualified.out3
    }
}

task apps_1269_missing_file_qualified {
    input {
        File task_missing_file
    }
    command <<<
        cat ~{task_missing_file}
    >>>
    output {
        Int out1 = 1
    }
}

task apps_1270_cloned_file_unqualified {
    input {
        File task_cloned_file_unqualified
    }
    command <<<
        cat ~{task_cloned_file_unqualified}
    >>>
    output {
        Int out2 = 2
    }
}

task apps_1270_cloned_file_qualified {
    input {
        File task_cloned_file_qualified
    }
    command <<<
        cat ~{task_cloned_file_qualified}
    >>>
    output {
        Int out3 = 3
    }
}