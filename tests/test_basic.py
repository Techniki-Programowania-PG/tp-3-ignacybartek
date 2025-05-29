from __future__ import annotations

import scikit_build_example as m


def test_version():
    assert m.__version__ == "0.0.1"


def test_add():
    assert m.add(1, 2) == 3


def test_sub():
    assert m.subtract(1, 2) == -1

def test_filtr1D():
    print(m.apply_filter([0, 0, 1, 1, 0, 0], [0.25, 0.5, 0.25]))

def test_filtr2D():
    matrix = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ]
    kernel = [
        [0, -1, 0],
        [-1, 5, -1],
        [0, -1, 0]
    ]
    result = m.apply_filter_2D(matrix, kernel)
    for row in result:
        print(row)


def main():
   print(m.add(1,2))
   return 0
