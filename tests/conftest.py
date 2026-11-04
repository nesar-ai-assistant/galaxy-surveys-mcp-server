"""pytest configuration: add custom markers."""
import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "network: needs internet (deselect with -m 'not network')")
