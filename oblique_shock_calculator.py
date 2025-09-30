
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Oblique Shock Relations Calculator (Perfect Gas)
Gamma default = 1.4. Angles in degrees.
Inputs supported (always provide M1, then one of the following):
  - Turn angle (Weak shock)
  - Turn angle (Strong shock)
  - Wave angle (beta)
  - M1n (normal component of M1)
Outputs:
  M2, Turn ang., Wave ang., p2/p1, rho2/rho1, T2/T1, p02/p01, M1n, M2n
"""
import math
from typing import Callable, Tuple

# ---------- Helpers ----------

def deg2rad(x): return math.radians(x)
def rad2deg(x): return math.degrees(x)

def theta_from_beta(M1: float, beta_rad: float, gamma: float) -> float:
    """Return deflection angle theta (radians) from theta-beta-M relation for given beta."""
    g = gamma
    Mn = M1 * math.sin(beta_rad)
    # Guard Mn <= 1 -> no attached oblique shock
    if Mn <= 1.0:
        return 0.0
    num = 2.0 * (Mn*Mn - 1.0)
    den = M1*M1 * (g + math.cos(2.0*beta_rad)) + 2.0
    theta = math.atan( (num / den) * (1.0 / math.tan(beta_rad)) )
    return theta

def oblique_from_M1_beta(M1: float, beta_rad: float, gamma: float = 1.4) -> dict:
    """Compute full oblique shock outputs from M1 and beta (radians)."""
    g = gamma
    gm1 = g - 1.0
    gp1 = g + 1.0
    a = gm1 / 2.0

    # Deflection
    theta = theta_from_beta(M1, beta_rad, g)

    # Normal components
    M1n = M1 * math.sin(beta_rad)
    if M1n <= 1.0:
        raise ValueError("Shock not attached (M1n <= 1). Increase beta or M1.")
    # Normal shock relations
    num = 1.0 + 0.5 * gm1 * M1n * M1n
    den = g * M1n * M1n - 0.5 * gm1
    if den <= 0:
        raise ValueError("Invalid post-shock state (denominator <= 0).")
    M2n = math.sqrt(num / den)

    p2_p1 = 1.0 + (2.0 * g / gp1) * (M1n * M1n - 1.0)
    rho2_rho1 = (gp1 * M1n * M1n) / (gm1 * M1n * M1n + 2.0)
    T2_T1 = p2_p1 / rho2_rho1

    # Tangential velocity unchanged across shock -> downstream M2 using geometry
    theta_rad = theta
    sin_term = math.sin(beta_rad - theta_rad)
    if sin_term <= 0:
        raise ValueError("Computed (beta - theta) not physically valid.")
    M2 = M2n / sin_term

    # Total pressure loss across shock
    p01_p1 = (1.0 + a * M1 * M1) ** (g / gm1)
    p02_p2 = (1.0 + a * M2 * M2) ** (g / gm1)
    p02_p01 = (p2_p1 * p02_p2) / p01_p1

    return {
        "M2": M2,
        "theta_deg": rad2deg(theta),
        "beta_deg": rad2deg(beta_rad),
        "p2_p1": p2_p1,
        "rho2_rho1": rho2_rho1,
        "T2_T1": T2_T1,
        "p02_p01": p02_p01,
        "M1n": M1n,
        "M2n": M2n,
    }

def theta_beta_brackets(M1: float, gamma: float) -> Tuple[float, float, float, float]:
    """
    Compute useful beta bounds (in radians) and approximate peak location:
    Returns (beta_min, beta_peak, beta_max, theta_max).
    beta_min ~ asin(1/M1)+eps, beta_max ~ pi/2 - eps.
    beta_peak is where theta(beta) is max (found by coarse search).
    """
    eps = 1e-6
    beta_min = math.asin(1.0 / M1) + eps
    beta_max = 0.5 * math.pi - eps

    # Coarse search for peak
    N = 500
    best_beta = beta_min
    best_theta = -1.0
    for i in range(N + 1):
        b = beta_min + (beta_max - beta_min) * i / N
        th = theta_from_beta(M1, b, gamma)
        if th > best_theta:
            best_theta = th
            best_beta = b
    return beta_min, best_beta, beta_max, best_theta

def solve_bisection(f: Callable[[float], float], lo: float, hi: float, tol: float = 1e-12, maxiter: int = 200) -> float:
    flo = f(lo); fhi = f(hi)
    if math.isnan(flo) or math.isnan(fhi):
        raise ValueError("Function returned NaN at bracket endpoints.")
    if flo == 0.0: return lo
    if fhi == 0.0: return hi
    if flo * fhi > 0:
        raise ValueError("No sign change in bracket for bisection.")
    for _ in range(maxiter):
        mid = 0.5*(lo+hi)
        fmid = f(mid)
        if fmid == 0.0 or abs(hi-lo) < tol:
            return mid
        if flo * fmid < 0:
            hi = mid; fhi = fmid
        else:
            lo = mid; flo = fmid
    return 0.5*(lo+hi)

def beta_from_theta(M1: float, theta_deg: float, branch: str, gamma: float) -> float:
    """Solve theta-beta-M for beta (radians) given M1 and theta (deg). branch in {'weak','strong'}."""
    theta = deg2rad(theta_deg)
    beta_min, beta_peak, beta_max, theta_max = theta_beta_brackets(M1, gamma)
    if theta > theta_max + 1e-12 or theta < 0:
        raise ValueError(f"Requested turn angle exceeds attached-shock limit. theta_max ≈ {rad2deg(theta_max):.6f} deg for M1={M1}.")
    def g(b): return theta_from_beta(M1, b, gamma) - theta
    eps = 1e-6
    if branch.lower().startswith("weak"):
        lo, hi = beta_min, beta_peak
    else:
        lo, hi = beta_peak, beta_max
    # Ensure sign change by adjusting small margins
    # Evaluate endpoints slightly inward if needed
    f_lo = g(lo + eps)
    f_hi = g(hi - eps)
    if f_lo == 0.0:
        return lo + eps
    if f_hi == 0.0:
        return hi - eps
    if f_lo * f_hi > 0:
        # Try scanning to find a sign change
        K = 200
        prev_b = lo + eps
        prev_f = g(prev_b)
        for i in range(1, K+1):
            b = lo + eps + (hi - lo - 2*eps) * i / K
            fb = g(b)
            if prev_f * fb <= 0:
                return solve_bisection(g, prev_b, b)
            prev_b, prev_f = b, fb
        raise ValueError("Could not bracket root for the chosen branch; check inputs.")
    return solve_bisection(g, lo + eps, hi - eps)

# ---------- Interface ----------

def compute_from_inputs(mode: str, M1: float, value: float, gamma: float = 1.4) -> dict:
    """
    mode in {
      'turn angle (weak)',
      'turn angle (strong)',
      'wave angle',
      'm1n'
    }
    value is theta_deg for 'turn angle' modes, beta_deg for 'wave angle', or M1n for 'm1n'.
    """
    mode_k = mode.strip().lower()
    if M1 <= 1.0:
        raise ValueError("M1 must be > 1 for an attached oblique shock.")
    if mode_k == "turn angle (weak)":
        beta = beta_from_theta(M1, value, "weak", gamma)
    elif mode_k == "turn angle (strong)":
        beta = beta_from_theta(M1, value, "strong", gamma)
    elif mode_k == "wave angle":
        beta = deg2rad(value)
        # Validate attached condition
        beta_min = math.asin(1.0 / M1) + 1e-9
        if not (beta_min < beta < 0.5*math.pi):
            raise ValueError("Wave angle must satisfy asin(1/M1) < beta < 90 deg for attached oblique shock.")
    elif mode_k == "m1n":
        M1n = value
        if M1n <= 1.0 or M1n >= M1:
            raise ValueError("M1n must be in (1, M1).")
        beta = math.asin(M1n / M1)
    else:
        raise ValueError("Unknown mode.")

    res = oblique_from_M1_beta(M1, beta, gamma)
    return res

def _format_output(res: dict, gamma: float) -> str:
    def fmt(x): return f"{x:.6f}"
    lines = [
        f"Oblique Shock Relations — Perfect Gas (gamma = {gamma}), angles in degrees",
        "OUTPUT:",
        f"M2 =	{fmt(res['M2'])}",
        f"Turn ang. =	{fmt(res['theta_deg'])}",
        f"Wave ang. =	{fmt(res['beta_deg'])}",
        f"p2/p1 =	{fmt(res['p2_p1'])}",
        f"rho2/rho1 =	{fmt(res['rho2_rho1'])}",
        f"T2/T1 =	{fmt(res['T2_T1'])}",
        f"p02/p01 =	{fmt(res['p02_p01'])}",
        f"M1n =	{fmt(res['M1n'])}",
        f"M2n =	{fmt(res['M2n'])}",
    ]
    return "\n".join(lines)

def main():
    print("Oblique Shock Relations — Perfect Gas (gamma = 1.4)")
    print("Always provide M1, then choose an input mode:")
    modes = [
        "1) Turn angle (Weak)",
        "2) Turn angle (Strong)",
        "3) Wave angle",
        "4) M1n",
    ]
    print("\n".join(modes))
    sel = input("Select mode by number or label: ").strip().lower()
    mapping = {
        "1": "Turn angle (Weak)",
        "2": "Turn angle (Strong)",
        "3": "Wave angle",
        "4": "M1n",
    }
    mode = mapping.get(sel, sel)
    try:
        M1 = float(input("Enter M1 (>1): ").strip())
        if mode.lower() == "turn angle (weak)" or mode.lower() == "turn angle (strong)":
            val = float(input("Enter turn angle (deg): ").strip())
        elif mode.lower() == "wave angle":
            val = float(input("Enter wave angle beta (deg): ").strip())
        elif mode.lower() == "m1n":
            val = float(input("Enter M1n (>1 and < M1): ").strip())
        else:
            raise ValueError("Unknown mode.")
        gamma_str = input("Gamma [default 1.4]: ").strip()
        gamma = float(gamma_str) if gamma_str else 1.4

        res = compute_from_inputs(mode, M1, val, gamma)
        print("\n" + _format_output(res, gamma))
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
