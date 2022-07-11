import click
import logging
from typing import Optional

from dxcint.Context import Context
from dxcint.Terraform import Terraform
from dxcint.TestDiscovery import TestDiscovery, TestDiscoveryError
from dxcint.LogFormatter import LogFormatter


@click.group()
@click.option(
    "--verbosity",
    default="info",
    type=click.Choice(["error", "warning", "info", "debug"]),
)
def dxcint(verbosity: str = "info") -> None:
    level = {"error": 40, "warning": 30, "info": 20, "debug": 10}[verbosity]
    log_format = '%(asctime)s | %(levelname)8s | %(message)s'
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(LogFormatter(log_format))
    logger.addHandler(ch)
    # TODO check build.sbt in the CWD


@dxcint.command()
@click.argument('dxc_repository_root', required=True)
@click.option(
    "-t",
    "--test_name",
    default=None,
    type=str,
    help="Name of a single test or suite of tests. If not provided - builds just the core compiler library "
         "(e.g. dxCompiler-VERSION.jar)",
)
def integration(
        dxc_repository_root: str,
        test_name: Optional[str]
) -> None:
    test_context = Context(project="dxCompiler_playground", repo_root=dxc_repository_root)
    test_discovery = TestDiscovery(test_context)
    if test_name:
        try:
            registered_tests = test_discovery.discover(test_name)
        except TestDiscoveryError:
            registered_tests = test_discovery.discover_single_test(test_name)
    else:
        registered_tests = []
    terraform = Terraform(
        languages={x.language for x in registered_tests},
        context=test_context,
        dependencies=test_discovery.discover_dependencies()
    )
    _ = terraform.build()
    for test in registered_tests:
        passed = test.test_result


@dxcint.command()
@click.argument('dxc_repository_root', required=True)
@click.option(
    "-d",
    "--directory",
    default=".",
    type=str,
    help="Extension of the test workflow script files. Allowed are 'cwl', 'cwl.json', 'wdl'",
)
@click.option(
    "-x",
    "--extension",
    default="wdl",
    type=str,
    help="Extension of the test workflow script files. Allowed are 'cwl', 'cwl.json', 'wdl'",
)
@click.option(
    "-s",
    "--suite",
    type=str,
    help="Test suite name. Usually a team-defined group of tests to be run in each CI/CD step",
)
@click.option(
    "-c",
    "--category",
    type=str,
    help="Test category name. Usually a team-defined category that reflects a type of the test",
)
def add(
        dxc_repository_root: str,
        directory: str,
        extension: str,
        suite: str,
        category: str
) -> None:
    test_context = Context(project="dxCompiler_playground", repo_root=dxc_repository_root)
    test_discovery = TestDiscovery(test_context)
    _ = test_discovery.add_tests(directory, extension, suite, category)
