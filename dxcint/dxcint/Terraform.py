from typing import List
from dxcint.RegisteredTest import RegisteredTest
from dxcint.Context import Context


class Terraform(object):
    def __init__(self, registered_tests: List[RegisteredTest], context: Context):
        self._languages = set(x.language for x in registered_tests)
# TODO STOPPED HERE