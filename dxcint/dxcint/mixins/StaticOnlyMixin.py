from dxcint.RegisteredTest import RegisteredTest


class StaticOnlyMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _compile_executable(*args, **kwargs) -> str:
        super()._compile_executable(
            additional_compiler_flags=["-instanceTypeSelection", "static"],
        )
