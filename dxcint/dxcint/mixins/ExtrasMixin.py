from dxcint.RegisteredTest import RegisteredTest
import os


class ExtrasMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._extras_suffix = "_extras.json"

    def _compile_executable(self, *args, **kwargs) -> str:
        extras_path = os.path.join(
            os.path.dirname(self._src_file), self._test_name + self._extras_suffix
        )
        kwargs["additional_compiler_flags"] = kwargs.get(
            "additional_compiler_flags", []
        ) + ["--extras", extras_path]
        return super()._compile_executable(*args, **kwargs)
