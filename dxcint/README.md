# dxCint

## Test suite architecture
The following architecture will be implemented to support integration tests for the compiler.
1. **CLI** orchestrates the flow 
2. **TestDiscovery** scans for config files, registers the tests by creating an array of RegisteredTest objects 
3. **Dependency** class sets up necessary environment with immutable state for **RegisteredTest**’s 
4. **RegisteredTest** has methods to compile test, run test and output results 
5. **Messenger** is attached to a **RegisteredTest** to monitor the completion of the test 
6. **Validator** contains methods to validate test results if **Messenger** shows that the test is finished 
7. **PlatformUtil** and **JsonUtil** have methods for platform interaction and JSON parsing, respectively


### Requirements
* New architecture for a test suite
* Common logic for internal and customer integration tests
* TDD for the suite functionality

### Non-Requirements
* Keep in mind possibilities for a unified CI/CD and possible interactions with sciprodbuild framework
