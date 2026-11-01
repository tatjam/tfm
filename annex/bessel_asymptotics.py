import sympy as sp
import numpy as np
from scipy.special import i0e, i1e, ive
import matplotlib.pyplot as plt

import matplotlib

matplotlib.use("kitcat")

κ = sp.symbols("kappa", positive=True)


def I_asym(ν, order):
    """
    Asymptotic series of I_ν(κ) as κ -> ∞
    Source: https://dlmf.nist.gov/10.17#E1 and https://dlmf.nist.gov/10.40
    """
    prefactor = sp.exp(κ) / sp.sqrt(2 * sp.pi * κ)
    series_sum = 1
    for k in range(1, order + 1):
        # 10.17.1 aₖ(ν), inner k iteration named j
        prod = 1
        for j in range(1, k + 1):
            prod *= 4 * ν**2 - (2 * j - 1) ** 2
        prod /= sp.factorial(k) * (8**k)
        term = (-1) ** k * prod / (κ**k)
        series_sum += term
    return prefactor * series_sum


# If we go less it just expands to 0, which is sort of true but only for
# extremely large κ due to the square root.
order = 2
ratio = I_asym(1, order) / I_asym(0, order)

ratio_series = ratio.series(κ, sp.oo, order + 1).removeO().simplify()
print("I1/I0 expansion:")
sp.pprint(ratio_series)

r = ratio_series
b1 = 1 - r
# recurrence formula from wikipedia: I2/I0 = 1 - (2/kappa) * I1/I0
b2 = 1 - (1 - (2 / κ) * r)

arg = b2 / (2 * b1) - 1

arg_series = arg.series(κ, sp.oo, order).removeO().simplify()
final = sp.acos(arg_series).series(κ, sp.oo, 1).removeO()
print("\nFinal arccos expansion:")
sp.pprint(final)


# Visualization to check
def real(κ: np.ndarray) -> np.ndarray:
    i0 = i0e(κ)
    i1 = i1e(κ)
    i2 = ive(2, κ)
    b1 = 1.0 - i1 / i0
    b2 = 1.0 - i2 / i0
    return np.arccos(b2 / (2.0 * b1) - 1.0)


def approx(κ: np.ndarray) -> np.ndarray:
    return np.sqrt(3 / κ)


κs = np.linspace(1, 100, 200)
a = real(κs)
b = approx(κs)

plt.style.use("dark_background")
plt.plot(κs, a)
plt.plot(κs, b)
plt.legend(["Real", "Approx"])
plt.show()

plt.plot(κs, np.abs(a - b))
plt.yscale("log")
plt.legend(["|Real - Approx|"])
plt.show()

κs_error = np.linspace(1, 100000, 1000)
error = np.abs(real(κs_error) - approx(κs_error))
thres = 1e-9
valid = error < thres
kmin = κs_error[np.min(np.nonzero(valid))]

print(f"Error < {thres} for κ > {kmin}")
