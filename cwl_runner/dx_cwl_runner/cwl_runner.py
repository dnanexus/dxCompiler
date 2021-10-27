#!/usr/bin/env python3
from dx_cwl_runner import arg_parser
from dx_cwl_runner.compiler import CwlCompiler
from dx_cwl_runner.dx import Dx
from dx_cwl_runner.utils import Log


def main():
    args = arg_parser.parse_args()
    log = Log(verbose=not args.quiet, dryrun=args.dryrun)
    compiler = CwlCompiler(Dx(log))
    compiler.run_test(args.processfile, args.jobfile, args.basedir, args.outdir)


if __name__ == "__main__":
    main()
