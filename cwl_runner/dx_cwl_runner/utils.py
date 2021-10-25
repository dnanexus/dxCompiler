import hashlib
import json
import logging
import re
import subprocess
import sys
import urllib.request

try:
    from packaging.version import parse as parse_version
except:
    parse_version = None


class Log:
    _log = None

    @staticmethod
    def init(file=None, verbose: bool = True, dryrun: bool = False):
        Log._log = Log(file, verbose, dryrun)
        return Log._log

    @staticmethod
    def get():
        if Log._log is None:
            Log._log = Log()
        return Log._log

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


def check_tool(tool, version_re, repo, log):
    try:
        if version_re is not None:
            current_version_str = re.fullmatch(
                version_re, run_cmd(f"{tool} --version")
            ).group(1)
        else:
            run_cmd(f"{tool} -h")
            return
    except:
        raise Exception(f"'{tool}' not found - please add it to your $PATH")

    try:
        with urllib.request.urlopen(
                f"https://api.github.com/repos/{repo}/releases/latest"
        ) as url:
            latest_release = json.loads(url.read().decode())
        latest_version_str = latest_release["tag"]
    except Exception as ex:
        log.warning(f"Could not determine latest {tool} version", ex)
        return

    if parse_version:
        latest_version = parse_version(latest_version_str)
        current_version = parse_version(current_version_str)
        if current_version < latest_version:
            log.warning(
                f"A newer version of {tool} ({latest_version_str}) is available; please upgrade "
                f"with 'pip install --upgrade {tool}'"
            )
    elif latest_version_str != current_version_str:
        log.warning(
            f"The installed version of {tool} ({current_version_str}) does not match the latest "
            f"version ({latest_version_str}); you may want to upgrade with "
            f"'pip install --upgrade {tool}'"
        )


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
