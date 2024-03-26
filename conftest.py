import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--printout", action="store", default="less", help="Debug printout option value: more or less"
    )


