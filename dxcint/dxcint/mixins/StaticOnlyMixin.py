from dxcint.RegisteredTest import RegisteredTest


class StaticOnlyMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _compile_executable(self, *args, **kwargs) -> str:
        super_kwargs = {
            **kwargs.get("additional_compiler_flags", {}),
            **{"-instanceTypeSelection": "static"},
        }
        return super()._compile_executable(
            *args, **{"additional_compiler_flags": super_kwargs}
        )
