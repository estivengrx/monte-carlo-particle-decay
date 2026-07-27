"""
Testing of the differential decay rate dΓ/dx and the 
radiative correction function f(x) in the muon Monte Carlo simulation,
to see how it behaves in the limits x -> 0+ and x -> 1-, and 
in general to see the domain of the decay rate.
"""
import numpy as np
import sympy as sp

from sys import path
from os.path import dirname, abspath

path.append(dirname(dirname(abspath(__file__))))
from src.python.muonMonteCarlo import muonMonteCarlo

# Callin the muonMonteCarlo class to access its methods for testing
muon_mc = muonMonteCarlo(n_events=100_000_000, include_radiative=True)

# This is the complete function for the differential decay rate dΓ/dx, including the x² factor and the radiative correction term f(x).
def dGamma_dx_real(x_val: float) -> float:
    """dΓ/dx tal cual la usa el Monte Carlo (ya trae su propio x² adentro)."""
    return muon_mc._dGamma_dx_differential_decay_rate(np.array([x_val]))[0]

# This is only the radiative correction term f(x), multiplied by x², isolating the pure radiative term.
def x2_times_f_radiative(x_val: float) -> float:
    """x² · f(x), aislando el término radiativo puro."""
    f_val = muon_mc._f_radiative_corrections(np.array([x_val]))[0]
    return (x_val**2) * f_val

# Using Richardson extrapolation because the functions are not defined at x=0 or x=1, 
# and we want to estimate the limits from the right (x -> 0+) and from the left (x -> 1-).
def richardson_limit(func, x0: float, side: str, h0: float = 1e-2, n_terms: int = 6):
    """
    Estimate the limit of a function as x approaches x0 using Richardson extrapolation, 
    from the right (side="right") or from the left (side="left").

    Parameters
    ----------
    func : callable
        The function for which to estimate the limit.
    x0 : float
        The point at which to estimate the limit.
    side : str
        "right" for approaching from the right (x -> x0+), "left" for approaching from the left (x -> x0-).
    h0 : float, optional
        The initial step size for the extrapolation. Default is 1e-2.
    n_terms : int, optional
        The number of terms to use in the Richardson extrapolation. Default is 6.
    """
    hs = [h0 / (2**k) for k in range(n_terms)] # step sizes halving each time
    sign = 1 if side == "right" else -1 

    raw_values = [func(x0 + sign * h) for h in hs] # evaluate the function at x0 + h (or x0 - h for left)

    for h, v in zip(hs, raw_values): # print of the raw values for each step size
        print(f"h = {h:12.3e} | valor = {v:.10e}")

    seq = [sp.Float(v) for v in raw_values]
    order = 1
    while len(seq) > 1: # continue Richardson extrapolation until we have a single estimate
        seq = [(2**order * seq[i+1] - seq[i]) / (2**order - 1) for i in range(len(seq) - 1)]
        order += 1 

    estimate = float(seq[0])
    print(f"Richardson estimate: {estimate:.10e}")
    return estimate


print("dΓ/dx complete with x² and radiative term, in the limits x -> 0+ and x -> 1-")
lim_dGamma_0 = richardson_limit(dGamma_dx_real, x0=0.0, side="right")
lim_dGamma_1 = richardson_limit(dGamma_dx_real, x0=1.0, side="left")

print("\nx² · f(x) (radiative term only), in the limits x -> 0+ and x -> 1-")
lim_f_0 = richardson_limit(x2_times_f_radiative, x0=0.0, side="right")
lim_f_1 = richardson_limit(x2_times_f_radiative, x0=1.0, side="left")

print("\nResults summary:")
print(f"lim x->0+  dΓ/dx           ≈ {lim_dGamma_0:.6e}")
print(f"lim x->1-  dΓ/dx           ≈ {lim_dGamma_1:.6e}")
print(f"lim x->0+  x²·f_radiative  ≈ {lim_f_0:.6e}")
print(f"lim x->1-  x²·f_radiative  ≈ {lim_f_1:.6e}")


# Now, let's check the behavior of dΓ/dx across the entire physical domain [0, 1], 
# focusing on regions where it might become negative, concluding the domain regions where the function is negative, if any.

def check_region(x_values: np.ndarray, label: str):
    """
    Check the differential decay rate dΓ/dx over a specified range of x values,
    and report the number of negative values found, along with their fraction and approximate region.

    Parameters
    ----------
    x_values : np.ndarray
        Array of x values to evaluate.
    label : str
        Label for the region being checked, for reporting purposes.
    """
    y = muon_mc._dGamma_dx_differential_decay_rate(x_values)

    n_total = len(y)
    n_negative = np.sum(y < 0)
    frac_negative = n_negative / n_total

    print(f"\nRegion: {label}")
    print(f"x range      : [{x_values.min():.6e}, {x_values.max():.6e}]")
    print(f"Total values : {n_total}")
    print(f"Total negative values : {n_negative}  ({frac_negative*100:.2f}%)")

    if n_negative > 0:
        idx_min = np.argmin(y)
        print(f"Minimum value     : {y[idx_min]:.6e}  at x = {x_values[idx_min]:.6e}")

        # Find the boundary of the negative region (largest x where y < 0)
        negative_mask = y < 0
        x_negative = x_values[negative_mask]
        print(f"Approximate negative region: x ∈ ({x_negative.min():.6e}, {x_negative.max():.6e})")
    else:
        print("No negative values found in this region.")

    return y, n_negative, frac_negative

# Using a uniform mesh over the full domain known to get a fair estimate of where the function is negative.
x_uniform = np.linspace(1e-8, 1 - 1e-8, 1_000_000)
y_uniform, n_neg_uniform, frac_neg_uniform = check_region(x_uniform, "Uniform mesh over full domain (fair estimate)")

positive_mask = y_uniform > 0
x_positive_region = x_uniform[positive_mask]
print(f"\n\nCONCLUSION: The function is positive for x ∈ [{x_positive_region.min():.6e}, {x_positive_region.max():.6e}]")

valor = dGamma_dx_real(4.770585e-03) # This is the value of x where the function starts behaving positevely

# This will be the value used in the main simulation as a threshold for x_min,
# to ensure that the differential decay rate is non-negative.
print(f"\n\nCheck: dΓ/dx at x=4.770585e-03: {valor:.6e}")