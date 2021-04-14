import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""This tool implements generic interface to run a 
                                                    Common Workflow Language tool or workflow from the
                                                    command line. To be implemented by each CWL compliant
                                                    execution platform for testing conformance to the 
                                                    standard and optionally for use by users.""")
    parser.add_argument('--outdir', default="./", type=str,
                        help='Output directory, defaults to the current directory')
    parser.add_argument('--quiet', default=False, action='store_true',
                        help='no diagnostic output')
    parser.add_argument('--version', action='version', version='%(prog)s 1.2',
                        # TODO: what version should be displayed?
                        help='report the name & version, then quit without further processing')
    parser.add_argument('processfile', help="""The CommandLineTool, ExpressionTool, or Workflow description
                                               to run. Optional if the jobfile has a `cwl:tool` field to indicate 
                                               which process description to run.""")
    parser.add_argument('jobfile', help="""The input job document.""")

    return parser.parse_args()
