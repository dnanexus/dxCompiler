# Developer Setup for FG

Sandbox project: project-Gz4gZ1Q0KPY0kbfGKFG2y5kX
Environment: prod
Region: aws:us-east-1

Let me know if you want to change the region; it should be in the same region as another project where you want to run test workflows.

I will add you as a contributor of this project.

From the dxCompiler repo folder, run

```bash
./scripts/clean_build.sh
```

This will build dxCompiler, including creating a copy of the runtime assets in the sandbox project. You can find the built dxCompiler jar file in the dxCompiler repo folder. Each time you run the script, the previous assets and local dxCompiler jar will be replaced.
