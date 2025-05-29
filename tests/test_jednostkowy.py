import sys
sys.path.append("build/Release")
import _core
import unittest
from cmath import isclose

class TestSignals(unittest.TestCase):
    def test_add(self):
        self.assertEqual(_core.add(2, 3), 5)

    def test_subtract(self):
        self.assertEqual(_core.subtract(10, 4), 6)

    def test_sin_signal(self):
        s = _core.sin_signal(1, 0, 1, 100)
        self.assertAlmostEqual(s[0], 0.0628, places=2)

    def test_cos_signal(self):
        s = _core.cos_signal(1, 0, 1, 100)
        self.assertAlmostEqual(s[0], 1.0, places=2)

    def test_square_signal(self):
        s = _core.square_signal(1, 0, 1, 100)
        self.assertIn(s[0], [0.0, 1.0])

    def test_sawtooth_signal(self):
        s = _core.sawtooth_signal(1, 0, 1, 100)
        self.assertTrue(0.0 <= s[0] <= 1.0)

    def test_dft_idft(self):
        signal = [complex(1, 0), complex(0, 0), complex(-1, 0), complex(0, 0)]
        dft = _core.DFT(signal)
        idft = _core.IDFT(dft)
        for a, b in zip(signal, idft):
            self.assertTrue(isclose(a.real, b.real, abs_tol=1e-5))

if __name__ == "__main__":
    unittest.main()
