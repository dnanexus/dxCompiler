#!/bin/bash -e

# Global variables
top_dir=""
version=""
dry_run=""
build_flags=""
staging_token=""
production_token=""
docker_password=""
docker_user=""

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


# Create a public docker image for dxCompiler that allows a simple command line
# invocation
function build_docker_image {
    cd $top_dir/scripts
    ln $top_dir/dxCompiler-${version}.jar .

    echo "building a docker image"
    sudo docker build --build-arg VERSION=${version} -t dnanexus/dxcompiler:${version} .

    echo "tagging as latest"
    sudo docker tag dnanexus/dxcompiler:${version} dnanexus/dxcompiler:latest

    echo "For the next steps to work you need to:"
    echo "(1) be logged into docker.io"
    echo "(2) have permissions to create a repository for dnanexus"
    echo $docker_password | sudo docker login -u $docker_user --password-stdin

    echo "pushing to docker hub"
    sudo docker push dnanexus/dxcompiler:${version}
    sudo docker push dnanexus/dxcompiler:latest
}

function usage_die
{
    echo "arguments: "
    echo "  --force: remove existing build artifacts, and build new ones"
    echo "  --dry-run: don't actually run anything"
    echo "  --staging-token <string>: an auth token for the staging environment"
    echo "  --production-token <string>: an auth token for the production environment"
    echo "  --docker-user <string>: docker user name"
    echo "  --docker-password <string>: docker password"
    echo "  --branch <string>: branch to build from (default=main)"
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
            --docker-user)
                docker_user=$2
                shift
                ;;
            --docker-password)
                docker_password=$2
                shift
                ;;
            --branch)
                target_branch=$2
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
    if [[ $docker_user == "" ]]; then
        echo "docker user name is missing"
        exit 1
    fi
    if [[ $docker_password == "" ]]; then
        echo "docker password is missing"
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
build_docker_image
