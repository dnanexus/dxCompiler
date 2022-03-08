class RegisteredTest(object):
    def __init__(self, category: str, test_name: str):
        self._category = category
        self._test_name = test_name

    @property
    def category(self):
        return self._category

    @property
    def name(self):
        return self._test_name
