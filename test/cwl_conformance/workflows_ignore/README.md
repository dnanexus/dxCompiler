These are workflow tests that we are currently ignoring.

## Open questions about whether the test is valid:

* cond-wf-002*: it seems the workflow output should be optional
* count-lines11-extra-step-wf-noET: the workflow input is optional but cat-tool input is not
* The expression tool returns `null` but the return type is `Any`:
    * count-lines11-null-step-wf
    * count-lines11-null-step-wf-noET
* input_schema-def: missing a.bam file referenced in input

## Non-required tests that we currently can't support

* Requires a step input to be linked to a workflow input and also to have a default value that is used if the workflow input is unspecified.
    * dynresreq-workflow-stepdefault
    * count-lines11-wf1 (inputs/results renamed with 'ignore_' prefix because the wf is still used for other tests)

## Tests that can't be run due to bugs in cwltool or cwljava

* parse error
    * record-in-secondaryFiles-wf
    * schemadef-wf
* anonymous processes with no ID: the auto-generated ID is not stable ([issue](https://github.com/common-workflow-language/cwltool/issues/1520))
    * count-lines2-wf
    * io-file-default-wf
    * io-int-default-tool-and-wf
    * io-int-default-wf
    * io-int-optional-wf
    * io-int-wf
    * io-union-input-default-wf
    * no-inputs-wf
    * no-outputs-wf
    * output-arrays-file-wf
    * output-arrays-int-wf
    * scatter-valuefrom-inputs-wf1
    * scatter-valuefrom-wf1
    * scatter-valuefrom-wf5
    * scatter-wf1.cwl.json step-valuefrom5-wf
    * steplevel-resreq
    * sum-wf-noET
    * sum-wf
    * timelimit2-wf
    * timelimit3-wf
    * wf-loadContents
    * wf-loadContents2
    * wf-loadContents3
    * wf-loadContents4