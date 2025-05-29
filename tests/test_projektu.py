import sys
sys.path.append("build/Release")
import _core
import time
from cmath import sin, pi

# Parametry sygnału
freq = 1
t_start = 0
t_end = 1
samples = 100

# 3. Sygnał sinusoidalny
sinus = _core.sin_signal(freq, t_start, t_end, samples)
print("sin_signal:", sinus[:5])
_core.plot_signal(sinus)
time.sleep(1)

# 4. Cosinusoidalny
cosinus = _core.cos_signal(freq, t_start, t_end, samples)
print("cos_signal:", cosinus[:5])
_core.plot_signal(cosinus)
time.sleep(1)

# 5. Sygnał prostokątny
square = _core.square_signal(freq, t_start, t_end, samples)
print("square_signal:", square[:5])
_core.plot_signal(square)
time.sleep(1)

# 6. Piłozębny
saw = _core.sawtooth_signal(freq, t_start, t_end, samples)
print("sawtooth_signal:", saw[:5])
_core.plot_signal(saw)
time.sleep(1)

# 7. DFT
# użyjmy sinusa z 4 próbek jako przykład
simple_signal = [complex(s, 0) for s in _core.sin_signal(freq, 0, 1, 4)]
dft_result = _core.DFT(simple_signal)
print("DFT:", dft_result)

# 8. IDFT
idft_result = _core.IDFT(dft_result)
print("IDFT:", idft_result)

# 9. Wykres IDFT
_core.plot_signal([x.real for x in idft_result])
