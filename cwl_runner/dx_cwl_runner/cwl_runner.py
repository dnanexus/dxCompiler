#!/usr/bin/env python3
import os

from dx_cwl_runner import arg_parser
from dx_cwl_runner import utils
from dx_cwl_runner import input_utils
from dx_cwl_runner import output_utils
from dx_cwl_runner.dx import Dx
from dx_cwl_runner.utils import Log


def run_cwl(dx_compiler, processfile, dx_input, dx: Dx):
    cmd = f"java -jar {dx_compiler} compile {processfile} -force -folder {dx.test_folder} " \
          f"-project {dx.current_dx_project} -locked -inputs {dx_input}"
    if dx.log.dryrun:
        dx.log.log(f"compile command: {cmd}")
        return None, None
    else:
        executable = utils.run_cmd(cmd)
        new_dx_input = input_utils.get_new_dx_input(dx_input)

        dx.log.debug(f"Running {executable} with {new_dx_input} input.", dx.log.verbose)
        job = utils.run_cmd(f"dx run {executable} -f {new_dx_input} -y --brief")

        dx.log.debug(f"Waiting for {job} to finish...", dx.log.verbose)
        utils.run_cmd(f"dx wait {job}")
        return job, utils.run_cmd(f"dx watch {job} --quiet")


def main():
    args = arg_parser.parse_args()
    log = Log(verbose=not args.quiet, dryrun=args.dryrun)
    dx = Dx(log)
    dx.check_outdir(args.outdir, create=not log.dryrun)

    log.debug("Creating DX inputs...", not args.quiet)
    inputs, input_processfile = input_utils.create_dx_input(args.jobfile, args.basedir, dx)
    process_file = args.processfile or input_processfile
    dx_input = input_utils.write_dx_input(inputs, args.jobfile, log)

    log.debug("Running process file...", not args.quiet)
    execution, execution_log = run_cwl(dx.compiler_jar, process_file, dx_input, dx)
    if execution_log is not None:
        log.debug(log)

    if not log.dryrun:
        log.debug("Creating outputs...", not args.quiet)
        results = output_utils.create_dx_output(execution, process_file, args.outdir, dx)
        output_utils.write_output(f"{os.path.splitext(process_file)[0]}_results.json", results)

    log.debug("Finished successfully!", not args.quiet)


if __name__ == "__main__":
    main()
