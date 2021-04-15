#!/usr/bin/env python3
import os
import sys
import arg_parser
import utils
import input_utils
import output_utils
import constants as C


def run_cwl(dxCompiler, processfile, dx_input, verbose=True):
    applet = utils.run_cmd(f"java -jar {dxCompiler} compile {processfile} -force -folder {C.test_dx_folder} -project {C.current_dx_project} -locked -inputs {dx_input}", returnOutput=True)
    new_dx_input = input_utils.get_new_dx_input(dx_input)

    utils.print_if_verbose(f"Running {applet} with {new_dx_input} input.", verbose)
    job = utils.run_cmd(f"dx run {applet} -f {new_dx_input} -y --brief", returnOutput=True)

    utils.print_if_verbose(f"Waiting for {job} to finish...", verbose)
    utils.run_cmd(f"dx wait {job}", returnOutput=True)
    return job, utils.run_cmd(f"dx watch {job} --quiet", returnOutput=True)


args = arg_parser.parse_args()
utils.check_outdir(args.outdir)

utils.print_if_verbose("Creating DX inputs...", not args.quiet)
inputs = input_utils.create_dx_input(args.jobfile)
dx_input = input_utils.write_dx_input(inputs, args.jobfile)

utils.print_if_verbose("Running process file...", not args.quiet)
# TODO: how to find dxcompiler name?
job, log = run_cwl("dxCompiler.jar", args.processfile, dx_input, verbose=not args.quiet)
utils.print_if_verbose(log, not args.quiet, file=sys.stderr)

utils.print_if_verbose("Creating outputs...", not args.quiet)
results = output_utils.create_dx_output(job, args.processfile, args.outdir)
output_utils.write_output(f"{os.path.splitext(args.processfile)[0]}_results.json", results)

utils.print_if_verbose("Finished successfully!", not args.quiet)
