#!/bin/bash

main() {
    # Set dx project context
    unset DX_WORKSPACE_ID
    dx cd $DX_PROJECT_CONTEXT_ID:

    # Download and run the test script from the project
    dx download test/test_dxcompiler.sh || \
        dx-jobutil-report-error "Test script expected at /test/test_dxcompiler.sh"
    chmod +x test_dxcompiler.sh
    ./test_dxcompiler.sh
}
