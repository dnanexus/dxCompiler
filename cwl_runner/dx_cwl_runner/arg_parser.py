import argparse
import os

from dx_cwl_runner.dx import Dx


def parse_args():
    parser = argparse.ArgumentParser(description="""This tool implements the cwl-runner interface to run a 
                                                    Common Workflow Language tool or workflow from the
                                                    command line. This is primarily used for testing conformance 
                                                    to the standard.""")
    parser.add_argument("--outdir", default=os.getcwd(), type=str,
                        help=f"Output directory, defaults {os.getcwd()}")
    parser.add_argument("--quiet", default=False, action="store_true",
                        help="No diagnostic output")
    parser.add_argument("--version", action="version", version=f"%(prog)s {Dx.compiler_version}",
                        help="Report the name & version, then quit without further processing")
    parser.add_argument("--basedir", default=os.getcwd(), type=str,
                        help=f"Directory to which relative paths will be resolved; defaults to {os.getcwd()}")
    parser.add_argument("--dryrun", action="store_true", default=False,
                        help="Print out commands and input files but don't actually do anything")
    parser.add_argument("processfile", nargs="?",
                        help="""The CommandLineTool, ExpressionTool, or Workflow description
                                to run. Optional if the jobfile has a `cwl:tool` field to indicate 
                                which process description to run.""")
    parser.add_argument("jobfile", help="""The input job document.""")
    return parser.parse_args()
