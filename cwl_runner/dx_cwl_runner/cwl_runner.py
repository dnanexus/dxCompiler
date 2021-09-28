#!/usr/bin/env python3
import os

from dx_cwl_runner import arg_parser, input_utils, output_utils
from dx_cwl_runner.dx import Dx
from dx_cwl_runner.utils import Log


def main():
    args = arg_parser.parse_args()
    log = Log.init(verbose=not args.quiet, dryrun=args.dryrun)
    dx = Dx(log)
    dx.check_outdir(args.outdir, create=not args.dryrun)

    # create a tempdir to write updated CWL and input files
    with dx.tempdir() as tmpdir:
        # Convert the jobfile into a dx inputs file; upload any local files/directories.
        # Check whether there are any hard-coded file/directory paths in the CWL; if so, 
        # upload them and replace with URIs.
        log.debug("Creating DNAnexus inputs...")
        process_file, dx_input_file = input_utils.create_dx_input(
            args.processfile, args.jobfile, args.basedir, tmpdir, dx, log
        )

        # Compile and run the CWL
        log.debug("Compiling and running process file...")
        execution, execution_log = dx.run_cwl(process_file, dx_input_file)
        # if execution_log is not None:
        #    log.debug()

    # Convert dx to CWL outputs
    if not args.dryrun:
        log.debug("Creating outputs...")
        results = output_utils.create_dx_output(
            execution, process_file, args.outdir, dx
        )
        output_utils.write_output(
            f"{os.path.splitext(process_file)[0]}_results.json", results
        )

    log.debug("Finished successfully!")


if __name__ == "__main__":
    main()
