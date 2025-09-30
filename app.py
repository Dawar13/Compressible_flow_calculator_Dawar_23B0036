
import math
from typing import Callable, Tuple
import pandas as pd
import streamlit as st

st.set_page_config(page_title="Compressible Flow Calculator", page_icon="🛩️", layout="wide")

# ======================= Core Utilities =======================

def solve_bisection(f: Callable[[float], float], a: float, b: float, tol: float = 1e-10, maxiter: int = 200) -> float:
    fa, fb = f(a), f(b)
    if math.isnan(fa) or math.isnan(fb):
        raise ValueError("Function returned NaN at bracket endpoints.")
    if fa == 0.0: return a
    if fb == 0.0: return b
    if fa * fb > 0:
        raise ValueError("Bisection bracket does not change sign.")
    lo, hi, flo, fhi = a, b, fa, fb
    for _ in range(maxiter):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if fmid == 0.0 or abs(hi - lo) < tol:
            return mid
        if flo * fmid < 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)

# ======================= Isentropic Flow =======================

def area_ratio(M: float, gamma: float = 1.4) -> float:
    if M <= 0:
        raise ValueError("Mach number must be positive.")
    gm1, gp1 = gamma - 1.0, gamma + 1.0
    exponent = gp1 / (2.0 * gm1)
    return (1.0 / M) * ((2.0 / gp1) * (1.0 + 0.5 * gm1 * M * M)) ** exponent

def prandtl_meyer(M: float, gamma: float = 1.4) -> float:
    if M < 1.0:
        return 0.0
    gm1, gp1 = gamma - 1.0, gamma + 1.0
    root = math.sqrt(M * M - 1.0)
    term1 = math.sqrt(gp1 / gm1) * math.atan(math.sqrt(gm1 / gp1) * root)
    term2 = math.atan(root)
    return math.degrees(term1 - term2)

def isentropic_props(M: float, gamma: float = 1.4) -> dict:
    if M <= 0:
        raise ValueError("Mach number must be positive.")
    gm1, gp1 = gamma - 1.0, gamma + 1.0
    a = gm1 / 2.0
    T_T0 = 1.0 / (1.0 + a * M * M)
    p_p0 = T_T0 ** (gamma / gm1)
    rho_rho0 = T_T0 ** (1.0 / gm1)

    Tstar_T0 = 2.0 / gp1
    pstar_p0 = Tstar_T0 ** (gamma / gm1)
    rhostar_rho0 = Tstar_T0 ** (1.0 / gm1)

    return {
        "M": M,
        "mu_deg": math.degrees(math.asin(1.0 / M)) if M >= 1 else float("nan"),
        "nu_deg": prandtl_meyer(M, gamma) if M >= 1 else 0.0,
        "p_p0": p_p0,
        "rho_rho0": rho_rho0,
        "T_T0": T_T0,
        "p_pstar": p_p0 / pstar_p0,
        "rho_rhostar": rho_rho0 / rhostar_rho0,
        "T_Tstar": T_T0 / Tstar_T0,
        "A_Astar": area_ratio(M, gamma),
    }

def M_from_T_T0(r: float, gamma: float = 1.4) -> float:
    if r <= 0 or r >= 1:
        raise ValueError("T/T0 must be in (0,1).")
    a = (gamma - 1.0) / 2.0
    return math.sqrt((1.0 / r - 1.0) / a)

def M_from_rho_rho0(r: float, gamma: float = 1.4) -> float:
    if r <= 0 or r >= 1:
        raise ValueError("rho/rho0 must be in (0,1).")
    gm1 = gamma - 1.0
    a = gm1 / 2.0
    return math.sqrt((r ** (-gm1) - 1.0) / a)

def M_from_A_Astar(area: float, branch: str, gamma: float = 1.4) -> float:
    if area < 1.0:
        raise ValueError("A/A* must be >= 1.")
    def f(M): return area_ratio(M, gamma) - area
    eps = 1e-8
    if branch == "Subsonic":
        return solve_bisection(f, eps, 1.0 - eps, tol=1e-12, maxiter=200)
    else:
        M_hi = 50.0
        while f(M_hi) < 0 and M_hi < 1e6:
            M_hi *= 2.0
        return solve_bisection(f, 1.0 + eps, M_hi, tol=1e-12, maxiter=200)

def M_from_mu(mu_deg: float) -> float:
    if mu_deg <= 0 or mu_deg > 90:
        raise ValueError("Mach angle must be in (0, 90].")
    mu = math.radians(mu_deg)
    s = math.sin(mu)
    if s <= 0:
        raise ValueError("Invalid Mach angle.")
    return 1.0 / s

def M_from_nu(nu_deg: float, gamma: float = 1.4) -> float:
    if nu_deg < 0:
        raise ValueError("P–M angle cannot be negative.")
    target = nu_deg
    def g(M): return prandtl_meyer(M, gamma) - target
    eps = 1e-8
    lo, hi = 1.0 + eps, 50.0
    while g(hi) < 0 and hi < 1e6:
        hi *= 2.0
    return solve_bisection(g, lo, hi, tol=1e-12, maxiter=200)

# ======================= Normal Shock =======================

def normal_shock_from_M1(M1: float, gamma: float = 1.4) -> dict:
    if M1 <= 1.0:
        raise ValueError("Upstream Mach number M1 must be > 1 for a normal shock.")
    g = gamma
    gm1 = g - 1.0
    gp1 = g + 1.0
    num = 1.0 + 0.5 * gm1 * M1 * M1
    den = g * M1 * M1 - 0.5 * gm1
    if den <= 0:
        raise ValueError("Invalid state (denominator <= 0).")
    M2 = math.sqrt(num / den)
    p2_p1 = 1.0 + (2.0 * g / gp1) * (M1 * M1 - 1.0)
    rho2_rho1 = (gp1 * M1 * M1) / (gm1 * M1 * M1 + 2.0)
    T2_T1 = p2_p1 / rho2_rho1
    a = gm1 / 2.0
    p01_p1 = (1.0 + a * M1 * M1) ** (g / gm1)
    p02_p2 = (1.0 + a * M2 * M2) ** (g / gm1)
    p02_p01 = (p2_p1 * p02_p2) / p01_p1
    p1_p02 = 1.0 / (p02_p01 * p01_p1)
    return {
        "M1": M1, "M2": M2, "p02_p01": p02_p01, "p1_p02": p1_p02,
        "p2_p1": p2_p1, "rho2_rho1": rho2_rho1, "T2_T1": T2_T1
    }

def M1_from_M2(target_M2: float, gamma: float = 1.4) -> float:
    if target_M2 <= 0 or target_M2 >= 1:
        raise ValueError("For a normal shock, M2 must be in (0,1).")
    def f(M1): return normal_shock_from_M1(M1, gamma)["M2"] - target_M2
    return solve_bisection(f, 1.0 + 1e-8, 20.0)

def M1_from_ratio(target: float, which: str, gamma: float = 1.4) -> float:
    which = which.lower().strip()
    if which == "p2/p1":
        if target < 1: raise ValueError("p2/p1 must be >= 1 across a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["p2_p1"] - target
    elif which == "rho2/rho1":
        max_rho = (gamma + 1.0) / (gamma - 1.0)
        if target < 1 or target > max_rho:
            raise ValueError(f"rho2/rho1 must be in [1, {max_rho:.6f}] for a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["rho2_rho1"] - target
    elif which == "t2/t1":
        if target < 1: raise ValueError("T2/T1 must be >= 1 across a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["T2_T1"] - target
    elif which == "p02/p01":
        if target <= 0 or target > 1: raise ValueError("p02/p01 is in (0,1] across a normal shock.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["p02_p01"] - target
    elif which == "p1/p02":
        if target <= 0: raise ValueError("p1/p02 must be positive.")
        def f(M1): return normal_shock_from_M1(M1, gamma)["p1_p02"] - target
    else:
        raise ValueError("Unknown ratio type for inversion.")
    return solve_bisection(f, 1.0 + 1e-8, 20.0)

# ======================= Oblique Shock =======================

def theta_from_beta(M1: float, beta_rad: float, gamma: float) -> float:
    g = gamma
    Mn = M1 * math.sin(beta_rad)
    if Mn <= 1.0:
        return 0.0
    num = 2.0 * (Mn*Mn - 1.0)
    den = M1*M1 * (g + math.cos(2.0*beta_rad)) + 2.0
    theta = math.atan( (num / den) * (1.0 / math.tan(beta_rad)) )
    return theta

def oblique_from_M1_beta(M1: float, beta_rad: float, gamma: float = 1.4) -> dict:
    g = gamma
    gm1 = g - 1.0
    gp1 = g + 1.0
    a = gm1 / 2.0
    theta = theta_from_beta(M1, beta_rad, g)
    M1n = M1 * math.sin(beta_rad)
    if M1n <= 1.0:
        raise ValueError("Shock not attached (M1n <= 1). Increase beta or M1.")
    num = 1.0 + 0.5 * gm1 * M1n * M1n
    den = g * M1n * M1n - 0.5 * gm1
    if den <= 0: raise ValueError("Invalid post-shock state.")
    M2n = math.sqrt(num / den)
    p2_p1 = 1.0 + (2.0 * g / gp1) * (M1n * M1n - 1.0)
    rho2_rho1 = (gp1 * M1n * M1n) / (gm1 * M1n * M1n + 2.0)
    T2_T1 = p2_p1 / rho2_rho1
    theta_rad = theta
    sin_term = math.sin(beta_rad - theta_rad)
    if sin_term <= 0: raise ValueError("Computed (beta - theta) not physically valid.")
    M2 = M2n / sin_term
    p01_p1 = (1.0 + a * M1 * M1) ** (g / gm1)
    p02_p2 = (1.0 + a * M2 * M2) ** (g / gm1)
    p02_p01 = (p2_p1 * p02_p2) / p01_p1
    return {
        "M2": M2, "theta_deg": math.degrees(theta), "beta_deg": math.degrees(beta_rad),
        "p2_p1": p2_p1, "rho2_rho1": rho2_rho1, "T2_T1": T2_T1, "p02_p01": p02_p01,
        "M1n": M1n, "M2n": M2n
    }

def theta_beta_brackets(M1: float, gamma: float) -> Tuple[float, float, float, float]:
    eps = 1e-6
    beta_min = math.asin(1.0 / M1) + eps
    beta_max = 0.5 * math.pi - eps
    N = 400
    best_beta, best_theta = beta_min, -1.0
    for i in range(N + 1):
        b = beta_min + (beta_max - beta_min) * i / N
        th = theta_from_beta(M1, b, gamma)
        if th > best_theta:
            best_theta, best_beta = th, b
    return beta_min, best_beta, beta_max, best_theta

def beta_from_theta(M1: float, theta_deg: float, branch: str, gamma: float) -> float:
    theta = math.radians(theta_deg)
    beta_min, beta_peak, beta_max, theta_max = theta_beta_brackets(M1, gamma)
    if theta > theta_max + 1e-12 or theta < 0:
        raise ValueError(f"Requested turn angle exceeds attached-shock limit. theta_max ≈ {math.degrees(theta_max):.6f} deg for M1={M1}.")
    def g(b): return theta_from_beta(M1, b, gamma) - theta
    eps = 1e-6
    if branch == "Weak":
        lo, hi = beta_min, beta_peak
    else:
        lo, hi = beta_peak, beta_max
    f_lo, f_hi = g(lo + eps), g(hi - eps)
    if f_lo == 0.0: return lo + eps
    if f_hi == 0.0: return hi - eps
    if f_lo * f_hi > 0:
        # scan for sign change
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

# ======================= UI Helpers =======================

def df_from_dict(d: dict, order: list) -> pd.DataFrame:
    rows = [(k, d[k]) for k in order if k in d]
    return pd.DataFrame(rows, columns=["Quantity", "Value"])

def show_results(df: pd.DataFrame, title: str):
    st.subheader(title)
    st.dataframe(df, hide_index=True, use_container_width=True)

# ======================= Sidebar =======================
st.title("🛩️ Compressible Flow Calculator")
st.caption("Isentropic flows, Normal shocks, and Oblique shocks (weak/strong). Perfect gas model.")

with st.sidebar:
    st.header("Global Settings")
    gamma = st.number_input("Gamma (γ)", min_value=1.01, max_value=2.0, value=1.4, step=0.01, format="%.2f")
    module = st.radio("Module", ["Isentropic", "Normal Shock", "Oblique Shock"], horizontal=False)

# ======================= Main Panels =======================

if module == "Isentropic":
    st.markdown("**Inputs:** Mach number, T/T₀, ρ/ρ₀, A/A*, Mach angle μ (deg), P–M angle ν (deg)")
    with st.form("iso_form"):
        itype = st.selectbox("Select input type", ["Mach number", "T/T0", "rho/rho0", "A/A* (Subsonic)", "A/A* (Supersonic)", "Mach angle (deg)", "P-M angle (deg)"])
        val = st.number_input("Enter value", value=2.0, step=0.1, format="%.6f")
        submitted = st.form_submit_button("Compute")
    if submitted:
        try:
            key = itype.lower()
            if key.startswith("mach number"):
                M = float(val)
            elif key == "t/t0":
                M = M_from_T_T0(float(val), gamma)
            elif key == "rho/rho0":
                M = M_from_rho_rho0(float(val), gamma)
            elif key.startswith("a/a* (sub"):
                M = M_from_A_Astar(float(val), "Subsonic", gamma)
            elif key.startswith("a/a* (sup"):
                M = M_from_A_Astar(float(val), "Supersonic", gamma)
            elif key.startswith("mach angle"):
                M = M_from_mu(float(val))
            elif key.startswith("p-m angle"):
                M = M_from_nu(float(val), gamma)
            else:
                raise ValueError("Unknown input type.")
            res = isentropic_props(M, gamma)
            order = ["M", "mu_deg", "nu_deg", "p_p0", "rho_rho0", "T_T0", "p_pstar", "rho_rhostar", "T_Tstar", "A_Astar"]
            df = df_from_dict(res, order)
            show_results(df, "Isentropic Flow — Results")
        except Exception as e:
            st.error(f"Error: {e}")

elif module == "Normal Shock":
    st.markdown("**Inputs:** M1, M2, p2/p1, ρ2/ρ1, T2/T1, p02/p01, p1/p02")
    with st.form("ns_form"):
        itype = st.selectbox("Select input type", ["M1", "M2", "p2/p1", "rho2/rho1", "T2/T1", "p02/p01", "p1/p02"])
        val = st.number_input("Enter value", value=2.0, step=0.1, format="%.6f")
        submitted = st.form_submit_button("Compute")
    if submitted:
        try:
            key = itype.lower()
            if key == "m1":
                M1 = float(val)
            elif key == "m2":
                M1 = M1_from_M2(float(val), gamma)
            elif key in {"p2/p1", "rho2/rho1", "t2/t1", "p02/p01", "p1/p02"}:
                M1 = M1_from_ratio(float(val), key, gamma)
            else:
                raise ValueError("Unknown input type.")
            res = normal_shock_from_M1(M1, gamma)
            order = ["M1", "M2", "p02_p01", "p1_p02", "p2_p1", "rho2_rho1", "T2_T1"]
            df = df_from_dict(res, order)
            show_results(df, "Normal Shock — Results")
        except Exception as e:
            st.error(f"Error: {e}")

else:  # Oblique Shock
    st.markdown("**Inputs:** M1 and one of (Turn angle — Weak/Strong, Wave angle β, or M1n). Outputs include both normal components.")
    with st.form("obl_form"):
        M1 = st.number_input("M1 (> 1)", value=2.0, step=0.1, format="%.6f")
        mode = st.selectbox("Mode", ["Turn angle (Weak)", "Turn angle (Strong)", "Wave angle (β)", "M1n"])
        if mode.startswith("Turn angle"):
            val = st.number_input("Turn angle θ (deg)", value=10.0, step=0.5, format="%.6f")
        elif mode.startswith("Wave angle"):
            val = st.number_input("Wave angle β (deg)", value=30.0, step=0.5, format="%.6f")
        else:
            val = st.number_input("M1n (normal component)", value=1.5, step=0.1, format="%.6f")
        submitted = st.form_submit_button("Compute")
    if submitted:
        try:
            if M1 <= 1.0: raise ValueError("M1 must be > 1 for an attached oblique shock.")
            if mode == "Turn angle (Weak)":
                beta = beta_from_theta(M1, val, "Weak", gamma)
            elif mode == "Turn angle (Strong)":
                beta = beta_from_theta(M1, val, "Strong", gamma)
            elif mode.startswith("Wave angle"):
                beta = math.radians(val)
                beta_min = math.asin(1.0 / M1) + 1e-9
                if not (beta_min < beta < 0.5*math.pi):
                    raise ValueError("Wave angle must satisfy asin(1/M1) < β < 90°.")
            else:
                M1n = val
                if M1n <= 1.0 or M1n >= M1:
                    raise ValueError("M1n must be in (1, M1).")
                beta = math.asin(M1n / M1)
            res = oblique_from_M1_beta(M1, beta, gamma)
            order = ["M2", "theta_deg", "beta_deg", "p2_p1", "rho2_rho1", "T2_T1", "p02_p01", "M1n", "M2n"]
            df = df_from_dict(res, order)
            show_results(df, "Oblique Shock — Results")
        except Exception as e:
            st.error(f"Error: {e}")

# ======================= Footer =======================
st.markdown("---")
st.caption("Made for quick aerospace coursework and design sanity checks. Use responsibly; verify critical results.")    
