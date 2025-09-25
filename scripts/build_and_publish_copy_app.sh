#!/bin/bash -e

# Global variables
top_dir=""
token=""
should_publish_app=""

# https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
# Get the source directory of the distribution
function get_top_dir {
    SCRIPT=$(realpath "$0")
    SCRIPTPATH=$(dirname "$SCRIPT")
    top_dir=$(realpath "$SCRIPTPATH/..")
    echo "The dxCompiler top directory is: $top_dir"
}

function build {
    # build the release on staging
    echo "building staging release"
    dx login --staging --token $staging_token --noprojects
    $top_dir/scripts/build_release.py --multi-region $build_flags

    ## test that it actually works
    echo "running multi region tests on staging"
    $top_dir/scripts/multi_region_tests.py
    #$top_dir/scripts/proxy_test.py

    echo "leave staging"
    dx clearenv

    ## build on production
    echo "building on production"
    dx login --token $production_token --noprojects
    $top_dir/scripts/build_release.py --multi-region $build_flags
}

function parse_cmd_line {
    while [[ $# -ge 1 ]]
    do
        case "$1" in
            --token)
                token=$2
                shift
                ;;
            --publish-app)
                should_publish_app=$2
                shift
                ;;
            *)
                echo "unknown argument $1"
                echo "arguments: "
                echo "  --token: dx auth token"
                echo "  --should_publish_app: either 'true' or 'false'"
                exit 1
        esac
        shift
    done

    if [[ $token == "" ]]; then
        echo "dx token is missing"
        exit 1
    fi

    if [[ $should_publish_app == "" ]]; then
        echo "should_publish_app is missing"
        exit 1
    fi

    if [[ $should_publish_app != "true" && $should_publish_app != "false" ]]; then
        echo "should_publish_app must be 'true' or 'false'"
        exit 1
    fi
}

# main program
parse_cmd_line $@
get_top_dir
build
