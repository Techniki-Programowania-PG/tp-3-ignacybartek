from __future__ import annotations

import scikit_build_example as m


def test_version():
    assert m.__version__ == "0.0.1"


def test_add():
    assert m.add(1, 2) == 3


def test_sub():
    assert m.subtract(1, 2) == -1

def test_filtr():
    print("Test apply_filter:")
    print(_core.apply_filter([0, 0, 1, 1, 0, 0], [0.25, 0.5, 0.25]))



def main():
   print(m.add(1,2))
   return 0
