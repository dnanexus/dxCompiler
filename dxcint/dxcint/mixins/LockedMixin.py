from dxcint.RegisteredTest import RegisteredTest


class LockedMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._locked = True
