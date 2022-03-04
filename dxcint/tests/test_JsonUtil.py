import pytest
from dxcint.JsonUtil import JsonUtil


@pytest.fixture(scope="session")
def json_util_init():
    yield JsonUtil()


def test_(json_util_init):
    assert False
