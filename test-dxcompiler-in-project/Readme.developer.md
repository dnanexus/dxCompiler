# test-dxcompiler-in-project Developer Readme

## Applet Requirements

Applet `test-dxcompiler-in-project` can be built and deployed in customer-provided projects to facilitate testing on their WDL workflows. It requires
- Java 8, to run `java -jar dxCompiler.jar compile <workflow>`
- `CONTRIBUTE` access to the project, to run the test workflows
- `VIEW` access to all projects, in case workflows reference objects in other projects
- Regional dxCompiler assets to already exist on the platform
