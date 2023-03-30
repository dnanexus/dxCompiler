import pytest
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.Context import ContextEmpty


@pytest.fixture
def init_ExpectedOutput(mock_messenger):
    yield ExpectedOutput(test_name="eotest", src_file="", category="eo", context=ContextEmpty)

