import os
from dxcint.utils import rm_prefix


class RegisteredTestError(Exception):
    """
    Class to handle RegisteredTest errors
    """


class RegisteredTest(object):
    def __init__(
            self,
            src_file: str,
            category: str,
            test_name: str
    ):
        self._src_file = src_file
        self._category = category
        self._test_name = test_name
        self._language = rm_prefix(os.path.basename(src_file), test_name)
        self._exec_id = None
        self._job_id = None

    @property
    def category(self):
        return self._category

    @property
    def name(self):
        return self._test_name

    @property
    def language(self):
        return self._language
