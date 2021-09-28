#!/usr/bin/env python3
from contextlib import contextmanager
import os
import shutil
import sys
import tempfile
from typing import Optional, Tuple

from dx_cwl_runner import arg_parser, input_utils, output_utils, utils
from dx_cwl_runner.dx import Dx
from dx_cwl_runner.utils import Log


def run_cwl(
    dx_compiler: str, processfile: str, dx_input_file: str, dx: Dx
) -> Tuple[Optional[str], Optional[str]]:
    cmd = (
        f"java -jar {dx_compiler} compile {processfile} -force -folder {dx.test_folder} "
        f"-project {dx.current_dx_project} -locked -inputs {dx_input_file}"
    )
    if dx.log.dryrun:
        dx.log.log(f"compile command: {cmd}")
        return None, None
    else:
        executable = utils.run_cmd(cmd)
        new_dx_input = input_utils.get_new_dx_input(dx_input_file)

        dx.log.debug(f"Running {executable} with {new_dx_input} input.", dx.log.verbose)
        job_id = utils.run_cmd(f"dx run {executable} -f {new_dx_input} -y --brief")

        dx.log.debug(f"Waiting for {job_id} to finish...", dx.log.verbose)
        utils.run_cmd(f"dx wait {job_id}")
        return job_id, utils.run_cmd(f"dx watch {job_id} --quiet")


@contextmanager
def tempdir():
    path = tempfile.mkdtemp()
    try:
        yield path
    finally:
        try:
            shutil.rmtree(path)
        except IOError:
            sys.stderr.write(f"Failed to clean up temp dir {path}")


def main():
    args = arg_parser.parse_args()
    log = Log(verbose=not args.quiet, dryrun=args.dryrun)
    dx = Dx(log)
    dx.check_outdir(args.outdir, create=not args.dryrun)

    # create a tempdir to write updated CWL and input files
    with tempdir() as tmpdir:
        # Convert the jobfile into a dx inputs file; upload any local files/directories.
        # Check whether there are any hard-coded file/directory paths in the CWL; if so, 
        # upload them and replace with URIs.
        log.debug("Creating DNAnexus inputs...")
        process_file, dx_input_file = input_utils.create_dx_input(
            args.processfile, args.jobfile, args.basedir, tmpdir, dx, log
        )

        # Compile and run the CWL
        log.debug("Compiling and running process file...")
        execution, execution_log = run_cwl(dx.compiler_jar, process_file, dx_input_file, dx)
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
