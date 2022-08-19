import pytest
from dxcint.Messenger import Messenger


@pytest.fixture(scope="session")
def messenger_init():
    yield Messenger()

# TODO how to poll messages
