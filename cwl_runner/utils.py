import subprocess
import os
import sys
import hashlib
import dxpy


def run_cmd(cmd: str, returnOutput: bool = False, print_output=False):
    output = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    if returnOutput:
        if print_output:
            print(output.strip())
        return output.strip().decode('utf-8')


def get_checksum(file_path):
    with open(file_path, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)
    return file_hash.hexdigest()


def get_file_name(outdir, file_id):
    platform_file_name = dxpy.DXFile(file_id).describe(fields={"name": True}).get("name")

    file_name = os.path.join(outdir, platform_file_name)
    counter = 0
    while os.path.exists(file_name):
        counter += 1
        file_name = f"{os.path.join(outdir, platform_file_name)}({counter})"
    return file_name


def upload_file(file, platform_location):
    new_file: dxpy.DXFile = dxpy.upload_local_file(file, folder=platform_location)
    return f"dx://{new_file.get_proj_id()}:{new_file.get_id()}"


def print_if_verbose(message, verbose=True, file=sys.stdout):
    if verbose:
        print(message, file=file)


def check_outdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
    elif not os.path.isdir(dir):
        raise Exception(f"Outdir {dir} is not a directory!")
    elif not os.access(dir, os.W_OK):
        raise Exception(f"You need write access to the outdir repository ({dir})!")
