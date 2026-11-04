"""pytest configuration: add custom markers."""


def pytest_configure(config):
    config.addinivalue_line("markers", "network: needs internet (deselect with -m 'not network')")
