from __future__ import annotations

import scikit_build_example as m


def test_version():
    assert m.__version__ == "0.0.1"


def test_add():
    assert m.add(1, 2) == 3


def test_sub():
    assert m.subtract(1, 2) == -1


def main():
    x={1,0,3,2,3,4,4};
    y={2,3,4,1,6,4,6};
    m.plot(x,y);
    return 0;