import hashlib
import subprocess
import sys
from typing import TextIO

class Log:
    @staticmethod
    def log(message: str, file=sys.stderr):
        print(message, file=file)

    def __init__(self, verbose: bool = True, dryrun: bool = False):
        self.verbose = verbose
        self.dryrun = dryrun

    def debug(self, message: str, file: TextIO = sys.stderr):
        if self.verbose or self.dryrun:
            Log.log(message, file)


def run_cmd(cmd: str, verbose: bool = False) -> str:
    output = subprocess.check_output(cmd, shell=True, executable='/bin/bash').strip().decode('utf-8')
    if verbose:
        Log.log(output)
    return output


def get_checksum(file_path: str) -> str:
    with open(file_path, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)
    return file_hash.hexdigest()
