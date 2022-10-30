# Synopsis
Here are some recommendations for debugging the workflows and individual stages compiled and run with dxCompiler on the 
DNAnexus platform.


## Debugging workflow stages
General approach follows a similar strategy as the debugging of any job on the platform. Consult the official [DNAnexus documentation](https://documentation.dnanexus.com/developer/apps/execution-environment#debugging-and-connecting-to-jobs-via-ssh) 
on the matter. Below is a step-by-step guide:
### Use case: debug a stage (job)
1. Find the job which failed. Use platform UI (`View Failure Cause` button) or write a program to parse the execution tree.
2. Run the stage by cloning the failed job:
```bash
dx run --clone job-xxxx --priority=high -y --debug-on All --ssh
```
A new `job-yyyy` will start and with an attempt to login to the worker. If an attempt fails, ssh directly to the `job-yyyy` 
```bash
dx ssh job-yyyy
```

3. Inside the running job, if you expect the job to succeed, but still want to debug it:
```bash
sudo touch /.dx-hold
```
However, if the job successfully finished, it gets to the `terminated` state without an option to get access to the worker.

4. To restart whatever processing was done in the container if it already finished, first switch to the root user
while preserving job environment:
```bash
sudo -E bash
```

5. The job script (from the `command` section of the workflow task/process) is located in `/home/dnanexus/meta/commandScript` for debugging.   Add a delay to this script (e.g. `sleep 3600`) if the containerized execution inside this script finishes too quickly to debug.  Then invoke the job script:
 ```bash
 /home/dnanexus/meta/containerRunScript
 ```

6. If Docker is used as an execution environment in your task, you can debug inside the Docker container while it's running, by using GNU screens:
   * create a new screen by pressing `Ctrl a` followed by `c`
   * in the new screen:
   ```bash
   sudo -E bash
   docker exec -it $(docker ps -q | head -n1) /bin/bash
   ```
