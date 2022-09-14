from dxcint.RegisteredTest import RegisteredTest


class LockedMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _compile_executable(self, *args, **kwargs) -> str:
        kwargs["additional_compiler_flags"] = kwargs.get(
            "additional_compiler_flags", []
        ) + ["-locked"]
        return super()._compile_executable(*args, **kwargs)
