# Synopsis
Here are some recommendations for debugging workflows and individual stages compiled and run with dxCompiler on the 
DNAnexus platform.


## Debugging workflow stages
General approach follows the same strategy as the debugging of any jobs on the platform. Consult the official [DNAnexus documentation](https://documentation.dnanexus.com/developer/apps/execution-environment#debugging-and-connecting-to-jobs-via-ssh) 
on the matter. Below is a step-by-step guide
### Use case: debug a failing stage (job)
1. Find the job which failed. Use platform UI (`View Failure Cause` button) or write your own program to parse the execution tree.
2. 