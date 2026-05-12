#!/bin/bash -e

# Global variables
top_dir=""
version=""
dry_run=""
build_flags=""
staging_token=""
production_token=""
skip_clone_regions=""

# https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
# Get the source directory of the distribution
function get_top_dir {
    SCRIPT=$(realpath "$0")
    SCRIPTPATH=$(dirname "$SCRIPT")
    top_dir=$(realpath "$SCRIPTPATH/..")
    echo "The dxCompiler top directory is: $top_dir"
}

function get_version {
    # figure out the release tag
    local config=$top_dir/core/src/main/resources/application.conf
    version=$(grep version "${config}" | cut --delimiter='"' --fields=2)
    if [ -z "$version" ]; then
        echo "could not figure out the dxCompiler release version"
        exit 1
    fi
    echo "dxCompiler version is $version"
}

function basic_checks {
    # make sure dx is in our path
    local path_to_dx=$(which dx)
    if [ -z "$path_to_dx" ] ; then
        echo "Could not find the dx CLI"
        exit 1
    fi
    echo "Found the dx CLI: $path_to_dx"

    local branch=$(git symbolic-ref --short HEAD)
    if [[ $branch != "$target_branch" ]]; then
        echo "This isn't $target_branch branch, please do 'git checkout $target_branch'"
        exit 1
    fi

    echo "making sure $target_branch is up to date"
    git pull
}

function build {
    local skip_regions_flag=""
    if [[ -n "$skip_clone_regions" ]]; then
        skip_regions_flag="--skip-clone-regions $skip_clone_regions"
    fi

    # build the release on staging
    echo "building staging release"
    dx login --staging --token $staging_token --noprojects
    $top_dir/scripts/build_release.py --multi-region $build_flags $skip_regions_flag

    ## test that it actually works
    echo "running multi region tests on staging"
    if [[ -n "$skip_clone_regions" ]]; then
        $top_dir/scripts/multi_region_tests.py --skip-regions $skip_clone_regions
    else
        $top_dir/scripts/multi_region_tests.py
    fi
    #$top_dir/scripts/proxy_test.py

    echo "leave staging"
    dx clearenv

    ## build on production
    echo "building on production"
    dx login --token $production_token --noprojects
    $top_dir/scripts/build_release.py --multi-region $build_flags $skip_regions_flag
}

function usage_die
{
    echo "arguments: "
    echo "  --force: remove existing build artifacts, and build new ones"
    echo "  --dry-run: don't actually run anything"
    echo "  --staging-token <string>: an auth token for the staging environment"
    echo "  --production-token <string>: an auth token for the production environment"
    echo "  --branch <string>: branch to build from (default=main)"
    echo "  --skip-clone-regions <string>: comma-separated regions to skip asset cloning (e.g. aws:eu-west-2)"
    exit 1
}

function parse_cmd_line {
    while [[ $# -ge 1 ]]
    do
        case "$1" in
            --force)
                build_flags="$build_flags --force"
                ;;
            --dry-run|--dry_run|--dryrun)
                dry_run=1
                build_flags="$build_flags --dry-run"
                ;;
            --staging-token)
                staging_token=$2
                shift
                ;;
            --production-token)
                production_token=$2
                shift
                ;;
            --branch)
                target_branch=$2
                shift
                ;;
            --skip-clone-regions)
                skip_clone_regions=$2
                shift
                ;;
            *)
                echo "unknown argument $1"
                usage_die
        esac
        shift
    done

    if [[ $dry_run == "1" ]]; then
        return
    fi

    if [[ $staging_token == "" ]]; then
        echo "staging token is missing"
        exit 1
    fi
    if [[ $production_token == "" ]]; then
        echo "production token is missing"
        exit 1
    fi
    if [[ $target_branch == "" ]]; then
        target_branch="main"
    fi
}

# main program
parse_cmd_line $@
basic_checks
get_top_dir
get_version
build
