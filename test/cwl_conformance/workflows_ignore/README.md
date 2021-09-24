These are workflow tests that we are currently ignoring.

## Open questions about whether the test is valid:

* cond-wf-002*: it seems the workflow output should be optional
* count-lines11-extra-step-wf-noET: the workflow input is optional but cat-tool input is not
* The expression tool returns `null` but the return type is `Any`:
    * count-lines11-null-step-wf
    * count-lines11-null-step-wf-noET
* input_schema-def: missing a.bam file referenced in input
* cond-wf-011* output type is `array[array[array[string]]]`, but the third-level array can contain nulls and thus should be optional

## Tests that can't be run due to bugs in cwltool or cwljava

* parse error
    * record-in-secondaryFiles-wf
    * schemadef-wf
* cwltool `--single-step` does not work for steps in nested workflows ([issue](https://github.com/common-workflow-language/cwltool/issues/1530))
    * count-lines8-wf
    * count-lines8-wf-noET
* cwltool `--single-step` does not apply inherited requirements and hints
    * env-wf2
    * env-wf3
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
    * iwdr_with_nested_dirs