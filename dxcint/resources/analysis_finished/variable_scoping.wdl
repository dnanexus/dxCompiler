task scoping_task {
    Int number1
    Int number2 = 20

    command {
        echo "${number1} and ${number2}"
    }

    output { String out = read_string(stdout()) }
    runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
}

workflow scoping_wf {
    Int number1 = 5
    Int number2 = 10

    # Assign Workflow's number2 variable to Task's number1:
    call scoping_task { input: number1 = number2 }
}
