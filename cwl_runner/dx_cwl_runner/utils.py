import hashlib
import logging
import subprocess
import sys


class Log:
    _log = None

    @staticmethod
    def init(file=None, verbose: bool = True, dryrun: bool = False):
        Log._log = Log(file, verbose, dryrun)
        return Log._log

    @staticmethod
    def get():
        if Log._log:
            return Log._log
        else:
            raise Exception("must call Log.init")

    def __init__(self, file=None, verbose: bool = True, dryrun: bool = False):
        if file:
            config = {"filename": file}
        else:
            config = {"stream": sys.stderr}
        if verbose or dryrun:
            config["level"] = logging.DEBUG
        else:
            config["level"] = logging.INFO
        logging.basicConfig(**config)
        self.logger = logging.getLogger()
        self.verbose = verbose
        self.dryrun = dryrun

    def __getattr__(self, item):
        getattr(self.logger, item)


def run_cmd(cmd: str, verbose: bool = False) -> str:
    output = (
        subprocess.check_output(cmd, shell=True, executable="/bin/bash")
        .strip()
        .decode("utf-8")
    )
    if verbose:
        Log.log(output)
    return output


def get_checksum(file_path: str) -> str:
    with open(file_path, "rb") as f:
        file_hash = hashlib.md5()
        while True:
            chunk = f.read(8192)
            if chunk:
                file_hash.update(chunk)
            else:
                break
    return file_hash.hexdigest()
