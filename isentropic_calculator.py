
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Isentropic Flow Relations Calculator (Perfect Gas)
Gamma default = 1.4, angles in degrees.
Supports inputs: M, T/T0, rho/rho0, A/A* (sub and sup branches), Mach angle mu, and Prandtl-Meyer angle nu.
Prints: Mach number, Mach angle, P-M angle, p/p0, rho/rho0, T/T0, p/p*, rho/rho*, T/T*, A/A*.
"""
import math
from typing import Callable

# --------- Core thermodynamic / gas-dynamics relations ---------

def isentropic_props(M: float, gamma: float = 1.4) -> dict:
    """
    Return a dict of isentropic relations for a given Mach number M.
    Includes p/p0, rho/rho0, T/T0, p/p*, rho/rho*, T/T*, A/A*, mu (deg), nu (deg).
    """
    if M <= 0:
        raise ValueError("Mach number must be positive.")
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    a = gm1 / 2.0

    # Stagnation relations
    T_T0 = 1.0 / (1.0 + a * M * M)
    p_p0 = T_T0 ** (gamma / gm1)
    rho_rho0 = T_T0 ** (1.0 / gm1)

    # Star (sonic) ratios (relative to same stagnation state)
    Tstar_T0 = 2.0 / gp1
    pstar_p0 = Tstar_T0 ** (gamma / gm1)
    rhostar_rho0 = Tstar_T0 ** (1.0 / gm1)

    T_Tstar = (T_T0) / Tstar_T0
    p_pstar = (p_p0) / pstar_p0
    rho_rhostar = (rho_rho0) / rhostar_rho0

    # Area-Mach relation (A/A*)
    A_Astar = area_ratio(M, gamma)

    # Mach angle (only for supersonic)
    if M >= 1.0:
        mu_deg = math.degrees(math.asin(1.0 / M))
        nu_deg = prandtl_meyer(M, gamma)  # already in degrees
    else:
        mu_deg = float("nan")
        nu_deg = 0.0  # conventionally 0 at M=1 and -> 0 as M->1-

    return {
        "M": M,
        "mu_deg": mu_deg,
        "nu_deg": nu_deg,
        "p_p0": p_p0,
        "rho_rho0": rho_rho0,
        "T_T0": T_T0,
        "p_pstar": p_pstar,
        "rho_rhostar": rho_rhostar,
        "T_Tstar": T_Tstar,
        "A_Astar": A_Astar,
    }


def area_ratio(M: float, gamma: float = 1.4) -> float:
    """A/A* as a function of Mach number (isentropic, perfect gas)."""
    if M <= 0:
        raise ValueError("Mach number must be positive.")
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    a = gm1 / 2.0
    exponent = gp1 / (2.0 * gm1)
    return (1.0 / M) * ((2.0 / gp1) * (1.0 + a * M * M)) ** exponent


def prandtl_meyer(M: float, gamma: float = 1.4) -> float:
    """
    Prandtl-Meyer function nu(M) in DEGREES for M >= 1.
    nu(M) = sqrt((g+1)/(g-1)) * atan( sqrt((g-1)/(g+1) * (M^2 - 1)) ) - atan( sqrt(M^2 - 1) )
    """
    if M < 1.0:
        return 0.0
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    root = math.sqrt((M * M) - 1.0)
    term1 = math.sqrt(gp1 / gm1) * math.atan(math.sqrt(gm1 / gp1) * root)
    term2 = math.atan(root)
    nu = term1 - term2
    return math.degrees(nu)


# --------- Inversion helpers (solve for M from various inputs) ---------

def solve_bisection(f: Callable[[float], float], a: float, b: float, tol: float = 1e-10, maxiter: int = 100) -> float:
    """Robust bisection (f(a)*f(b) must be <= 0)."""
    fa = f(a)
    fb = f(b)
    if math.isnan(fa) or math.isnan(fb):
        raise ValueError("Function returned NaN at the bracket endpoints.")
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0:
        raise ValueError("Bisection bracket does not change sign.")
    lo, hi = a, b
    flo, fhi = fa, fb
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


def M_from_T_T0(r: float, gamma: float = 1.4) -> float:
    """Invert T/T0 = 1 / (1 + (g-1)/2 M^2)."""
    if r <= 0 or r >= 1:
        raise ValueError("T/T0 must be in (0,1).")
    a = (gamma - 1.0) / 2.0
    return math.sqrt((1.0 / r - 1.0) / a)


def M_from_rho_rho0(r: float, gamma: float = 1.4) -> float:
    """Invert rho/rho0 = (1 + (g-1)/2 M^2)^(-1/(g-1))."""
    if r <= 0 or r >= 1:
        raise ValueError("rho/rho0 must be in (0,1).")
    gm1 = gamma - 1.0
    a = gm1 / 2.0
    return math.sqrt((r ** (-gm1) - 1.0) / a)


def M_from_A_Astar(area: float, branch: str, gamma: float = 1.4) -> float:
    """
    Invert A/A* to Mach on the specified branch ('sub' or 'sup').
    Valid only for area >= 1.
    """
    if area < 1.0:
        raise ValueError("A/A* must be >= 1.")
    branch = branch.lower().strip()
    def f(M): return area_ratio(M, gamma) - area
    eps = 1e-8
    if branch.startswith("sub"):
        # A/A* decreases from +inf at M->0+ down to 1 at M=1-
        return solve_bisection(f, eps, 1.0 - eps, tol=1e-12, maxiter=200)
    elif branch.startswith("sup"):
        # A/A* increases from 1 at M=1+ upward with M
        M_hi = 50.0
        while f(M_hi) < 0 and M_hi < 1e6:
            M_hi *= 2.0
        return solve_bisection(f, 1.0 + eps, M_hi, tol=1e-12, maxiter=200)
    else:
        raise ValueError("branch must be 'sub' or 'sup'.")


def M_from_mu(mu_deg: float) -> float:
    """Invert Mach angle mu (deg) to M: mu = asin(1/M). Valid for 0 < mu <= 90 deg."""
    if mu_deg <= 0 or mu_deg > 90:
        raise ValueError("Mach angle must be in (0, 90].")
    mu = math.radians(mu_deg)
    s = math.sin(mu)
    if s <= 0:
        raise ValueError("Invalid Mach angle.")
    return 1.0 / s


def M_from_nu(nu_deg: float, gamma: float = 1.4) -> float:
    """
    Invert Prandtl-Meyer angle nu (deg) to M >= 1 using bisection.
    nu increases monotonically from 0 at M=1 to a finite maximum as M->infinity.
    """
    if nu_deg < 0:
        raise ValueError("P-M angle cannot be negative.")
    target = nu_deg
    nu_max = prandtl_meyer(1e6, gamma)  # practical asymptote
    if target > nu_max + 1e-9:
        raise ValueError(f"P-M angle exceeds maximum ({nu_max:.6f} deg) for gamma={gamma}.")
    def g(M): return prandtl_meyer(M, gamma) - target
    eps = 1e-8
    lo, hi = 1.0 + eps, 50.0
    while g(hi) < 0 and hi < 1e6:
        hi *= 2.0
    return solve_bisection(g, lo, hi, tol=1e-12, maxiter=200)


# --------- Interface helpers ---------

def compute_from_input(input_type: str, value: float, gamma: float = 1.4) -> dict:
    """
    Compute full isentropic output from any supported input.
    input_type in {'M','T/T0','rho/rho0','A/A* sub','A/A* sup','Mach angle','P-M angle'}
    """
    key = input_type.strip().lower()
    if key in {"m", "mach", "mach number"}:
        M = float(value)
    elif key in {"t/t0", "t over t0", "t_t0"}:
        M = M_from_T_T0(float(value), gamma)
    elif key in {"rho/rho0", "ρ/ρ0", "rho_over_rho0", "rho_rho0"}:
        M = M_from_rho_rho0(float(value), gamma)
    elif key in {"a/a* sub", "a/astar sub", "area sub"}:
        M = M_from_A_Astar(float(value), "sub", gamma)
    elif key in {"a/a* sup", "a/astar sup", "area sup"}:
        M = M_from_A_Astar(float(value), "sup", gamma)
    elif key in {"mach angle", "mu", "μ"}:
        M = M_from_mu(float(value))
    elif key in {"p-m angle", "pm angle", "prandtl-meyer angle", "ν", "nu"}:
        M = M_from_nu(float(value), gamma)
    else:
        raise ValueError("Unknown input_type.")
    return isentropic_props(M, gamma)


def _format_output(res: dict, gamma: float) -> str:
    def fmt(x):
        if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
            return "—"
        return f"{x:.6f}"
    lines = [
        f"Isentropic Flow Relations — Perfect Gas (gamma = {gamma}), angles in degrees",
        "OUTPUT:",
        f"Mach number =	{fmt(res['M'])}",
        f"Mach angle (deg) =	{fmt(res['mu_deg'])}",
        f"P-M angle (deg) =	{fmt(res['nu_deg'])}",
        f"p/p0 =	{fmt(res['p_p0'])}",
        f"rho/rho0 =	{fmt(res['rho_rho0'])}",
        f"T/T0 =	{fmt(res['T_T0'])}",
        f"p/p* =	{fmt(res['p_pstar'])}",
        f"rho/rho* =	{fmt(res['rho_rhostar'])}",
        f"T/T* =	{fmt(res['T_Tstar'])}",
        f"A/A* =	{fmt(res['A_Astar'])}",
    ]
    return "\n".join(lines)


def main():
    print("Isentropic Flow Relations — Perfect Gas (gamma = 1.4), angles in degrees")
    print("INPUT options (enter the label exactly or choose the number):")
    options = [
        "1) Mach number",
        "2) T/T0",
        "3) rho/rho0",
        "4) A/A* sub",
        "5) A/A* sup",
        "6) Mach angle",
        "7) P-M angle",
    ]
    print("\n".join(options))
    sel = input("Select input type by number or label: ").strip().lower()
    mapping = {
        "1": "Mach number",
        "2": "T/T0",
        "3": "rho/rho0",
        "4": "A/A* sub",
        "5": "A/A* sup",
        "6": "Mach angle",
        "7": "P-M angle",
    }
    input_type = mapping.get(sel, sel)  # allow direct label input
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
