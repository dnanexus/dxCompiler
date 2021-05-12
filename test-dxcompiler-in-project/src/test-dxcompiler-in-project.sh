#!/bin/bash

main() {
    # Download and run the test script from the project
    dx download test/test_wdl.sh || dx-jobutil-report-error "Test script expected at /test/test_wdl.sh"
    chmod +x test_wdl.sh
    ./test_wdl.sh
}
