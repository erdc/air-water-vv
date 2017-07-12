import pytest

def pytest_addoption(parser):
    parser.addoption("--runfast", action="store_true", help="run fast tests")
