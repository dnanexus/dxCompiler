import hashlib
import json
import logging
import os
import re
import subprocess
import sys
import traceback
from typing import Optional
import urllib.request

try:
    from packaging.version import parse as parse_version
except:
    parse_version = None


class Log:
    def __init__(self, file=None, verbose: bool = True, dryrun: bool = False):
        if file:
            config = {"filename": file}
        else:
            config = {"stream": sys.stderr}
        if verbose or dryrun:
            config["level"] = logging.DEBUG
        else:
            config["level"] = logging.WARNING
        logging.basicConfig(**config)
        self.logger = logging.getLogger()
        self.verbose = verbose
        self.dryrun = dryrun

    def __getattr__(self, item):
        return getattr(self.logger, item)


class ErrorLog:
    _log = None
    verbose = False

    @staticmethod
    def get():
        if ErrorLog._log is None:
            ErrorLog._log = ErrorLog()
        return ErrorLog._log

    @staticmethod
    def log(msg, *args, **kwargs):
        sys.stderr.write(msg)
        exc_info = kwargs.get("exc_info")
        if exc_info:
            if isinstance(exc_info, BaseException):
                exc_info = (type(exc_info), exc_info, exc_info.__traceback__)
            elif not isinstance(exc_info, tuple):
                exc_info = sys.exc_info()
            traceback.print_exception(*exc_info)

    def __getattr__(self, item):
        if item == "error":
            return self.log
        else:
            raise NotImplementedError(item)


def check_tool(tool, version_re, repo, log, lib=None):
    try:
        if version_re is not None:
            current_version_str = re.fullmatch(
                version_re, run_cmd(f"{tool} --version", log)
            ).group(1)
        else:
            run_cmd(f"{tool} -h", log)
            return
    except:
        raise Exception(f"'{tool}' not found - please add it to your $PATH")

    # try getting the latest release first, and check tags if that fails
    try:
        with urllib.request.urlopen(
                f"https://api.github.com/repos/{repo}/releases/latest"
        ) as url:
            latest_release = json.loads(url.read().decode())
            latest_version_str = latest_release["tag_name"]
    except:
        try:
            with urllib.request.urlopen(
                    f"https://api.github.com/repos/{repo}/tags"
            ) as url:
                tags = json.loads(url.read().decode())
                latest_version_str = tags[0]["name"]
        except:
            log.warning(f"Could not determine latest {lib or tool} version", exc_info=True)
            return

    if parse_version:
        latest_version = parse_version(latest_version_str)
        current_version = parse_version(current_version_str)
        if current_version < latest_version:
            log.warning(
                f"A newer version of {lib or tool} ({latest_version_str}) is available; please "
                f"upgrade with 'pip install --upgrade {lib or tool}'"
            )
    elif latest_version_str != current_version_str:
        log.warning(
            f"The installed version of {lib or tool} ({current_version_str}) does not match the "
            f"latest version ({latest_version_str}); you may want to upgrade with "
            f"'pip install --upgrade {lib or tool}'"
        )


def run_cmd(cmd: str, log: Optional[Log] = None, verbose: Optional[bool] = False, cwd=None) -> str:
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        shell=True,
        executable="/bin/bash",
        check=False,
        capture_output=True,
        text=True
    )
    if proc.returncode == 0:
        logger = log.debug if log and log.verbose else print if verbose else None
        if logger:
            if proc.stdout:
                logger(proc.stdout)
            if proc.stderr:
                logger(proc.stderr)
        return proc.stdout.strip()
    else:
        raise Exception(
            f"Command '{cmd}' returned non-zero exit status {proc.returncode}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )


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


def get_unique_path(dir, filename):
    file_name = os.path.join(dir, filename)
    counter = 0
    while os.path.exists(file_name):
        counter += 1
        file_name = os.path.join(dir, f"{filename}.{counter}")
    return file_name
