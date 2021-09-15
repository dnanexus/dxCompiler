These are workflow tests that we are currently ignoring.

## Open questions about whether the test is valid:

* cond-wf-002*: it seems the workflow output should be optional
* count-lines11-extra-step-wf-noET: the workflow input is optional but cat-tool input is not

## Tests that can't be run due to bugs in cwltool or cwljava

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