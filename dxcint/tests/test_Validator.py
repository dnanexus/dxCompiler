import pytest
from dxcint.Validator import Validator


@pytest.fixture(scope="session")
def validator_init(some_registered_test):
    yield Validator(some_registered_test)

# TODO how to connect to messenger
