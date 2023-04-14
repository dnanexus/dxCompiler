from dxcint.RegisteredTest import RegisteredTest
import os


class ExtrasMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._extras_suffix = "_extras.json"
        self._extras = self._get_extras_src()

    def _get_extras_src(self) -> str:
        extras_basename = f"{self._test_name}{self._extras_suffix}"
        extras_src = os.path.join(os.path.dirname(self._src_file), extras_basename)
        return extras_src

    def _compile_executable(self, *args, **kwargs) -> str:
        super_kwargs = {
            **kwargs.get("additional_compiler_flags", {}),
            **{"-extras": self._extras},
        }
        return super()._compile_executable(
            *args, **{"additional_compiler_flags": super_kwargs}
        )
