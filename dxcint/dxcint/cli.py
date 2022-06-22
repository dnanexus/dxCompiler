import click
import logging
import sys
from typing import Optional

from dxcint.Context import Context


@click.group()
@click.option(
    "--verbosity",
    default="info",
    type=click.Choice(["error", "warning", "info", "debug"]),
)
def dxcint(verbosity: str) -> None:
    # TODO check build.sbt in the CWD
    pass


@dxcint.command()
@click.argument("location", required=True)
@click.option(
    "-t",
    "--test_name",
    default=None,
    type=str,
    help="Name of a single test or a suite of tests. If not provided - only builds a core compiler library "
         "(e.g. dxCompiler-VERSION.jar)",
)
def integration(
        location: str,
        test_name: Optional[str]
) -> None:
    # test_discovery = TestDiscovery(location)
    cont = Context(project="dxCompiler_playground", folder="/gvaihir/dxcint_testing")
    print(cont.repo_root_dir)


@dxcint.command()
@click.option(
    "-t",
    "--test_name",
    default=None,
    type=str,
    help="test name. If not provided - only builds a core compiler library (e.g. dxCompiler-VERSION.jar)",
)
def client(test_name: Optional[str]) -> None:
    pass
