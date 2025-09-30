
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Normal Shock Relations Calculator (Perfect Gas)
Gamma default = 1.4.
INPUT options: M1, M2, p2/p1, rho2/rho1, T2/T1, p02/p01, p1/p02
OUTPUT: M1, M2, p02/p01, p1/p02, p2/p1, rho2/rho1, T2/T1
"""
import math
from typing import Callable

# ---------------- Core normal-shock relations ----------------

def normal_shock_from_M1(M1: float, gamma: float = 1.4) -> dict:
    if M1 <= 1.0:
        raise ValueError("Upstream Mach number M1 must be > 1 for a normal shock.")
    g = gamma
    gm1 = g - 1.0
    gp1 = g + 1.0

    # Downstream Mach number
    # M2^2 = [1 + 0.5*(g-1)*M1^2] / [g*M1^2 - 0.5*(g-1)]
    num = 1.0 + 0.5 * gm1 * M1 * M1
    den = g * M1 * M1 - 0.5 * gm1
    if den <= 0:
        raise ValueError("Invalid state (denominator <= 0).")
    M2 = math.sqrt(num / den)

    # Static pressure ratio
    p2_p1 = 1.0 + (2.0 * g / gp1) * (M1 * M1 - 1.0)

    # Density ratio
    rho2_rho1 = (gp1 * M1 * M1) / (gm1 * M1 * M1 + 2.0)

    # Temperature ratio
    T2_T1 = p2_p1 / rho2_rho1

    # Stagnation pressure ratios (isentropic to/from states)
    a = gm1 / 2.0
    # p0/p = (1 + a M^2)^{g/(g-1)}
    p01_p1 = (1.0 + a * M1 * M1) ** (g / gm1)
    p02_p2 = (1.0 + a * M2 * M2) ** (g / gm1)

    # Total pressure ratio across shock
    p02_p01 = (p2_p1 * p02_p2) / p01_p1

    # p1/p02
    p1_p02 = 1.0 / (p02_p01 * (p01_p1 / 1.0))  # since p1/p01 = 1/p01_p1

    return {
        "M1": M1,
        "M2": M2,
        "p02_p01": p02_p01,
        "p1_p02": p1_p02,
        "p2_p1": p2_p1,
        "rho2_rho1": rho2_rho1,
        "T2_T1": T2_T1,
    }

# ---------------- Numerical inversion helpers ----------------

def solve_bisection(f: Callable[[float], float], lo: float, hi: float, tol: float = 1e-12, maxiter: int = 200) -> float:
    flo = f(lo)
    fhi = f(hi)
    if math.isnan(flo) or math.isnan(fhi):
        raise ValueError("Function returned NaN at bracket endpoints.")
    # If not bracketing, try to expand hi up to a reasonable limit
    expand_count = 0
    while flo * fhi > 0 and hi < 1e6 and expand_count < 50:
        hi *= 1.5
        fhi = f(hi)
        expand_count += 1
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0:
        raise ValueError("Could not bracket a root for the provided input.")
    for _ in range(maxiter):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if fmid == 0.0 or abs(hi - lo) < tol:
            return mid
        if flo * fmid < 0:
            hi = mid
            fhi = fmid
        else:
            lo = mid
            flo = fmid
    return 0.5 * (lo + hi)

def M1_from_M2(target_M2: float, gamma: float = 1.4) -> float:
    if target_M2 <= 0 or target_M2 >= 1:
        raise ValueError("For a normal shock, M2 must be in (0,1).")
    def f(M1):
        return normal_shock_from_M1(M1, gamma)["M2"] - target_M2
    return solve_bisection(f, 1.0 + 1e-8, 10.0)

def M1_from_ratio(target: float, which: str, gamma: float = 1.4) -> float:
    which = which.lower().strip()
    if which == "p2/p1":
        if target < 1:
            raise ValueError("p2/p1 must be >= 1 across a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["p2_p1"] - target
    elif which == "rho2/rho1":
        # Upper bound tends to (g+1)/(g-1), lower bound -> 1 at M1->1+
        max_rho = (gamma + 1.0) / (gamma - 1.0)
        if target < 1 or target > max_rho:
            raise ValueError(f"rho2/rho1 must be in [1, {(max_rho):.6f}] for a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["rho2_rho1"] - target
    elif which == "t2/t1":
        if target < 1:
            raise ValueError("T2/T1 must be >= 1 across a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["T2_T1"] - target
    elif which == "p02/p01":
        if target <= 0 or target > 1:
            raise ValueError("p02/p01 is in (0,1] across a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["p02_p01"] - target
    elif which == "p1/p02":
        if target <= 0:
            raise ValueError("p1/p02 must be positive.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["p1_p02"] - target
    else:
        raise ValueError("Unknown ratio type for inversion.")
    return solve_bisection(f, 1.0 + 1e-8, 20.0)

# ---------------- Interface ----------------

def compute_from_input(input_type: str, value: float, gamma: float = 1.4) -> dict:
    """
    Compute full normal-shock outputs from any supported single input.
    input_type in {'M1','M2','p2/p1','rho2/rho1','T2/T1','p02/p01','p1/p02'}
    """
    key = input_type.strip().lower()
    if key in {"m1", "mach1", "mach 1", "upstream mach", "upstream mach number"}:
        M1 = float(value)
    elif key in {"m2", "mach2", "mach 2", "downstream mach", "downstream mach number"}:
        M1 = M1_from_M2(float(value), gamma)
    elif key in {"p2/p1", "p2_p1"}:
        M1 = M1_from_ratio(float(value), "p2/p1", gamma)
    elif key in {"rho2/rho1", "r2/r1", "rho2_rho1"}:
        M1 = M1_from_ratio(float(value), "rho2/rho1", gamma)
    elif key in {"t2/t1", "t2_t1"}:
        M1 = M1_from_ratio(float(value), "t2/t1", gamma)
    elif key in {"p02/p01", "p0_2/p0_1", "p02_p01"}:
        M1 = M1_from_ratio(float(value), "p02/p01", gamma)
    elif key in {"p1/p02", "p1_p02"}:
        M1 = M1_from_ratio(float(value), "p1/p02", gamma)
    else:
        raise ValueError("Unknown input_type.")
    return normal_shock_from_M1(M1, gamma)

def _format_output(res: dict, gamma: float) -> str:
    def fmt(x): return f"{x:.6f}"
    lines = [
        f"Normal Shock Relations — Perfect Gas (gamma = {gamma})",
        "OUTPUT:",
        f"M1 =	{fmt(res['M1'])}",
        f"M2 =	{fmt(res['M2'])}",
        f"p02/p01 =	{fmt(res['p02_p01'])}",
        f"p1/p02 =	{fmt(res['p1_p02'])}",
        f"p2/p1 =	{fmt(res['p2_p1'])}",
        f"rho2/rho1 =	{fmt(res['rho2_rho1'])}",
        f"T2/T1 =	{fmt(res['T2_T1'])}",
    ]
    return "\n".join(lines)

def main():
    print("Normal Shock Relations — Perfect Gas (gamma = 1.4)")
    print("INPUT options (enter the label exactly or choose the number):")
    options = [
        "1) M1",
        "2) M2",
        "3) p2/p1",
        "4) rho2/rho1",
        "5) T2/T1",
        "6) p02/p01",
        "7) p1/p02",
    ]
    print("\n".join(options))
    sel = input("Select input type by number or label: ").strip().lower()
    mapping = {
        "1": "M1",
        "2": "M2",
        "3": "p2/p1",
        "4": "rho2/rho1",
        "5": "T2/T1",
        "6": "p02/p01",
        "7": "p1/p02",
    }
    input_type = mapping.get(sel, sel)
    try:
        val = float(input(f"Enter value for {input_type}: ").strip())
        gamma_str = input("Gamma [default 1.4]: ").strip()
        gamma = float(gamma_str) if gamma_str else 1.4

        res = compute_from_input(input_type, val, gamma)
        print("\n" + _format_output(res, gamma))
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
