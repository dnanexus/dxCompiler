import click
import logging
import sys
from typing import Optional

from dxcint.Context import Context
from dxcint.Terraform import Terraform
from dxcint.TestDiscovery import TestDiscovery, TestDiscoveryError


@click.group()
@click.option(
    "--verbosity",
    default="info",
    type=click.Choice(["error", "warning", "info", "debug"]),
)
def dxcint(verbosity: str = "info") -> None:
    level = {"error": 40, "warning": 30, "info": 20, "debug": 10}[verbosity]
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=level)
    logging.captureWarnings(True)
    # TODO check build.sbt in the CWD
    pass


@dxcint.command()
@click.option(
    "-t",
    "--test_name",
    default=None,
    type=str,
    help="Name of a single test or suite of tests. If not provided - builds just the core compiler library "
         "(e.g. dxCompiler-VERSION.jar)",
)
def integration(
        test_name: Optional[str]
) -> None:
    test_context = Context(project="dxCompiler_playground")
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
        test.validate()


