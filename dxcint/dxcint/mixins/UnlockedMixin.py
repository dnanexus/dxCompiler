from dxcint.RegisteredTest import RegisteredTest


class UnlockedMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _compile_executable(self, *args, **kwargs) -> str:
        super_kwargs = kwargs.get("additional_compiler_flags", {})
        if "-locked" in super_kwargs:
            _ = super_kwargs.pop("-locked")
        return super()._compile_executable(
            *args, **{"additional_compiler_flags": super_kwargs}
        )
