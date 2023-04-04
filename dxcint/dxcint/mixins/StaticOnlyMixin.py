from dxcint.RegisteredTest import RegisteredTest


class StaticOnlyMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def exec_id(self) -> str:
        if not self._exec_id:
            self._exec_id = self._compile_executable(
                {"-instanceTypeSelection": "static"}
            )
        return self._exec_id
