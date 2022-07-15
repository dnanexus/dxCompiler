import click
from typing import Optional

from dxcint.Context import Context, ContextEmpty
from dxcint.Terraform import Terraform
from dxcint.TestDiscovery import TestDiscovery, TestDiscoveryError


class CliContext(object):
    def __init__(self):
        self.verbosity = None


@click.group()
@click.pass_context
@click.option(
    "--verbosity",
    default="info",
    type=click.Choice(["error", "warning", "info", "debug"]),
)
def dxcint(ctx, verbosity: str = "info") -> None:
    ctx.obj = CliContext()
    ctx.obj.verbosity = verbosity


@dxcint.command()
@click.pass_context
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
        ctx,
        dxc_repository_root: str,
        test_name: Optional[str]
) -> None:
    test_context = Context(
        project="dxCompiler_playground",
        repo_root=dxc_repository_root,
        logger_verbosity=ctx.obj.verbosity
    )
    test_discovery = TestDiscovery(test_context)
    if test_name:
        try:
            registered_tests = test_discovery.discover(test_name)
        except TestDiscoveryError:
            registered_tests = test_discovery.discover_single_test(test_name)
    else:
        registered_tests = []
        test_context.logger.info(
            "CLI: No test/suite name provided. dxCint will only build the core dxCompiler executable"
        )
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
@click.argument("suite", required=True)
@click.argument("category", required=True)
@click.option(
    "-d",
    "--directory",
    default=".",
    type=str,
    help="Directory containing workflow fixtures to be added. DEFAULT: current WD",
)
@click.option(
    "-x",
    "--extension",
    default="wdl",
    type=str,
    help="Extension of the test workflow script files. Allowed are 'cwl', 'cwl.json', 'wdl'. DEFAULT: wdl",
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


@dxcint.command()
def suites() -> None:
    empty_context = ContextEmpty()
    test_discovery = TestDiscovery(empty_context)
    print(
        f"Available suites are:\n"
        f"{test_discovery.suites}"
    )


if __name__ == "__main__":
    dxcint()
