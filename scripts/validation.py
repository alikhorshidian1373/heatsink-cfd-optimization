"""
Independent analytical check of the CFD result.

A CFD number is only worth reporting if something outside the solver agrees
with it. Two checks are applied:

  1. Energy balance, taken from Fluent's own flux report (see README) -
     heat in 10.000001 W, net imbalance 2.3e-5 W. That verifies the solver
     conserved energy; it says nothing about whether the physics is right.

  2. A 1-D fin-array model built from a textbook correlation, compared against
     the CFD thermal resistance across the whole velocity sweep. This is the
     check that can actually catch a wrong answer.

Choice of correlation
---------------------
This heat sink sits in an open domain, not inside a duct: the fin channels are
short (L/Dh ~ 6) and most of the approaching air passes over and around the
array rather than through it. Under those conditions the fin surfaces behave
much more like flat plates growing fresh boundary layers than like a
fully-developed internal channel, so the laminar flat-plate correlation

    Nu_L = 0.664 · Re_L^(1/2) · Pr^(1/3)

is the right reference. (The fully-developed parallel-plate result, Nu = 7.54,
is also computed below - it is reported only to show how badly a ducted-flow
assumption would mis-predict the velocity trend here.)

The 1-D model is deliberately crude: it assumes one uniform surface
temperature, ignores three-dimensional flow acceleration around the array, and
ignores the leading-edge horseshoe vortex. Agreement within ~30 % with a
consistent sign is the pass criterion.

Run:  python scripts/validation.py
"""

import numpy as np
import pandas as pd
from pathlib import Path

# --- geometry (full heat sink, not the symmetric half) -----------------------
L_FIN = 0.060       # m, fin length along the flow
H_FIN = 0.025       # m, fin height
T_FIN = 0.0015      # m, fin thickness
N_FIN = 10          # fins
PITCH = 0.0062      # m, fin pitch
GAP = PITCH - T_FIN  # m, open channel width

L_BASE = 0.060      # m
W_BASE = 0.060      # m
T_BASE = 0.003      # m
K_AL = 202.0        # W/m·K, aluminium 6061

Q_TOTAL = 20.0      # W
T_INLET = 300.0     # K

# --- air: Fluent's built-in constant-property "air" --------------------------
RHO = 1.225         # kg/m³
MU = 1.7894e-5      # Pa·s
K_AIR = 0.0242      # W/m·K
CP = 1006.43        # J/kg·K
PR = MU * CP / K_AIR

ROOT = Path(__file__).resolve().parent.parent

A_FIN = N_FIN * 2 * L_FIN * H_FIN + N_FIN * L_FIN * T_FIN   # fin sides + tips
A_BASE = L_BASE * W_BASE - N_FIN * L_FIN * T_FIN            # exposed base top
A_TOTAL = A_FIN + A_BASE


def fin_efficiency(h):
    """Straight fin, uniform cross-section, adiabatic tip (corrected length)."""
    lc = H_FIN + T_FIN / 2
    m = np.sqrt(2 * h / (K_AL * T_FIN))
    return np.tanh(m * lc) / (m * lc)


def h_flat_plate(v):
    """Laminar flat plate, average over the fin length."""
    re_l = RHO * v * L_FIN / MU
    nu_l = 0.664 * re_l ** 0.5 * PR ** (1 / 3)
    return nu_l * K_AIR / L_FIN, re_l


def h_ducted(v):
    """Fully-developed parallel-plate channel - shown for contrast only."""
    dh = 2 * GAP
    nu = 7.54
    return nu * K_AIR / dh


def resistance(h):
    eta = fin_efficiency(h)
    ua = h * (eta * A_FIN + A_BASE)
    r_cond = T_BASE / (K_AL * L_BASE * W_BASE)   # 1-D spreading-free base conduction
    return 1 / ua + r_cond, eta


def main():
    print(f"Pr = {PR:.3f}   channel gap = {GAP*1000:.1f} mm   "
          f"wetted area = {A_TOTAL*1e4:.0f} cm²\n")

    cfd = pd.read_csv(ROOT / "data" / "velocity_sweep.csv").sort_values("velocity_m_s")
    cfd["R_CFD"] = (cfd.Tmax_K - T_INLET) / Q_TOTAL

    rows = []
    for v, r_cfd in zip(cfd.velocity_m_s, cfd.R_CFD):
        h, re_l = h_flat_plate(v)
        r_model, eta = resistance(h)
        r_duct, _ = resistance(h_ducted(v))
        rows.append({
            "V_m_s": v, "Re_L": re_l, "h_W_m2K": h, "eta_fin": eta,
            "R_flatplate": r_model, "R_ducted": r_duct, "R_CFD": r_cfd,
            "dev_%": 100 * (r_cfd - r_model) / r_model,
        })
    out = pd.DataFrame(rows)
    print(out.to_string(index=False, float_format=lambda v: f"{v:9.3f}"))

    dev = out["dev_%"]
    print(f"\nDeviation from the flat-plate model: {dev.min():.0f} % to {dev.max():.0f} %")
    print(f"CFD R_th scales as V^{np.polyfit(np.log(out.V_m_s), np.log(out.R_CFD), 1)[0]:.2f}; "
          f"flat-plate model V^{np.polyfit(np.log(out.V_m_s), np.log(out.R_flatplate), 1)[0]:.2f}; "
          f"ducted model V^{np.polyfit(np.log(out.V_m_s), np.log(out.R_ducted), 1)[0]:.2f}")

    ok = dev.abs().max() < 30 and (dev < 0).all()
    print("\nPASS - consistent sign, within 30 %: CFD predicts better cooling than the\n"
          "       1-D model, which is the expected direction (the model ignores flow\n"
          "       acceleration through the array and assumes an isothermal surface)."
          if ok else "\nREVIEW - deviation or sign is not as expected.")

    out.to_csv(ROOT / "data" / "validation.csv", index=False, float_format="%.5f")
    print(f"\nwrote {ROOT / 'data' / 'validation.csv'}")


if __name__ == "__main__":
    main()
