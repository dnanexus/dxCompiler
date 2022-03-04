import os.path

from dxcint.TestDiscovery import TestDiscovery
from dxcint.RegisteredTest import RegisteredTest


def test_discover(fixtures_dir):
    test_discovery = TestDiscovery(os.path.join(fixtures_dir, "mock_tests"))
    discovered_tests = test_discovery.discover()
    assert len(discovered_tests) == 2
    assert all([isinstance(x, RegisteredTest) for x in discovered_tests])
