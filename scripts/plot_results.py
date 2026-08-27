"""
Post-processing and optimisation analysis for the air-cooled heat-sink CFD study.

Reads  data/velocity_sweep.csv   (one row per inlet velocity, from Fluent)
Writes figures/*.png and figures/*.svg
       data/derived_results.csv  (every quantity used in the report)

Run:   python scripts/plot_results.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from vizstyle import (apply_style, save, title, annotate, markers,
                      BLUE, ORANGE, AQUA, INK, INK_2, INK_MUTED, SURFACE)

# --- case constants ----------------------------------------------------------
Q_TOTAL = 20.0          # W, heat dissipated by the full heat sink
T_INLET = 300.0         # K, air inlet temperature
W_DOMAIN = 0.260        # m, full-model domain width  (2 x 0.130 half-model)
H_DOMAIN = 0.128        # m, domain height
A_INLET = W_DOMAIN * H_DOMAIN   # m^2, full-model inlet area

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)


def load():
    df = pd.read_csv(DATA / "velocity_sweep.csv").sort_values("velocity_m_s")

    df["dT_K"] = df["Tmax_K"] - T_INLET
    df["R_th_K_per_W"] = df["dT_K"] / Q_TOTAL
    # Volumetric flow and the fan work needed to push it through the domain.
    df["Q_flow_m3_s"] = df["velocity_m_s"] * A_INLET
    df["P_pump_W"] = df["dP_Pa"] * df["Q_flow_m3_s"]
    return df


# --- figure 1: thermal resistance --------------------------------------------
def fig_thermal(df):
    fig, ax = plt.subplots()
    markers(ax, df.velocity_m_s, df.R_th_K_per_W, BLUE)

    # Label only the endpoints and the knee - never every point.
    lo, hi = df.iloc[0], df.iloc[-1]
    annotate(ax, lo.velocity_m_s, lo.R_th_K_per_W,
             f"{lo.R_th_K_per_W:.2f}", dx=8, dy=4, color=BLUE, weight="600")
    annotate(ax, hi.velocity_m_s, hi.R_th_K_per_W,
             f"{hi.R_th_K_per_W:.2f} K/W", dx=-4, dy=-18, color=BLUE, weight="600")

    ax.set_xlabel("Inlet velocity  (m/s)")
    ax.set_ylabel("Thermal resistance  R$_{th}$  (K/W)")
    title(ax, "Cooling improves fast, then barely at all",
          f"R$_{{th}}$ = ΔT / Q,  Q = {Q_TOTAL:.0f} W,  T$_{{in}}$ = {T_INLET:.0f} K")
    ax.set_ylim(bottom=0)
    save(fig, FIGS / "01_thermal_resistance")


# --- figure 2: pressure drop -------------------------------------------------
def fig_pressure(df):
    fig, ax = plt.subplots()

    # A quadratic reference curve: fully-developed duct loss scales with V^2.
    k = np.polyfit(df.velocity_m_s ** 2, df.dP_Pa, 1)[0]
    vv = np.linspace(df.velocity_m_s.min(), df.velocity_m_s.max(), 100)
    ax.plot(vv, k * vv ** 2, color=INK_MUTED, linewidth=1.4,
            linestyle=(0, (4, 3)), zorder=1, label=f"quadratic fit  {k:.3f}·V²")

    markers(ax, df.velocity_m_s, df.dP_Pa, ORANGE, label="CFD")

    hi = df.iloc[-1]
    annotate(ax, hi.velocity_m_s, hi.dP_Pa, f"{hi.dP_Pa:.2f} Pa",
             dx=-6, dy=-20, color=ORANGE, weight="600")

    ax.set_xlabel("Inlet velocity  (m/s)")
    ax.set_ylabel("Pressure drop  Δp  (Pa)")
    title(ax, "Pressure drop follows the classic V² law",
          "Area-weighted static pressure at the inlet, outlet held at 0 Pa gauge")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left")
    save(fig, FIGS / "02_pressure_drop")


# --- figure 3: pumping power -------------------------------------------------
def fig_pumping(df):
    fig, ax = plt.subplots()

    k = np.polyfit(df.velocity_m_s ** 3, df.P_pump_W, 1)[0]
    vv = np.linspace(df.velocity_m_s.min(), df.velocity_m_s.max(), 100)
    ax.plot(vv, k * vv ** 3, color=INK_MUTED, linewidth=1.4,
            linestyle=(0, (4, 3)), zorder=1, label=f"cubic fit  {k:.4f}·V³")

    markers(ax, df.velocity_m_s, df.P_pump_W, ORANGE, label="CFD")

    hi = df.iloc[-1]
    annotate(ax, hi.velocity_m_s, hi.P_pump_W, f"{hi.P_pump_W:.2f} W",
             dx=-8, dy=-20, color=ORANGE, weight="600")

    ax.set_xlabel("Inlet velocity  (m/s)")
    ax.set_ylabel("Pumping power  P$_{pump}$  (W)")
    title(ax, "The cost of air grows with the cube of speed",
          "P$_{pump}$ = Δp · V · A$_{inlet}$  —  the reason more airflow stops paying off")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left")
    save(fig, FIGS / "03_pumping_power")


# --- figure 4: the trade-off (Pareto front) ----------------------------------
def fig_pareto(df):
    fig, ax = plt.subplots()
    markers(ax, df.P_pump_W, df.R_th_K_per_W, BLUE)

    # Every point is a design choice, so every point earns its label here.
    for _, r in df.iterrows():
        annotate(ax, r.P_pump_W, r.R_th_K_per_W, f"{r.velocity_m_s:.0f} m/s",
                 dx=9, dy=4, color=INK_2)

    ax.set_xlabel("Pumping power  (W)")
    ax.set_ylabel("Thermal resistance  R$_{th}$  (K/W)")
    title(ax, "Buying the last few degrees costs the most power",
          "Each marker is one operating point; down-and-left is better")
    ax.set_xlim(0, df.P_pump_W.max() * 1.18)   # room for the right-most label
    ax.set_ylim(0, df.R_th_K_per_W.max() * 1.12)
    save(fig, FIGS / "04_tradeoff_pareto")


# --- figure 5: combined objective -> the optimum -----------------------------
def fig_objective(df, weight=1.0):
    """Normalised objective  J = R/R_ref + w · P_pump/P_ref.

    Both terms are made dimensionless against the mid-sweep case so the weight
    is the only judgement call: w = 1 treats a 1 % thermal gain and a 1 % power
    increase as equally important.
    """
    ref = df.iloc[len(df) // 2]
    j = df.R_th_K_per_W / ref.R_th_K_per_W + weight * df.P_pump_W / ref.P_pump_W

    # Smooth the discrete sweep so the minimum can be read off a curve.
    c = np.polyfit(df.velocity_m_s, j, 3)
    vv = np.linspace(df.velocity_m_s.min(), df.velocity_m_s.max(), 400)
    jj = np.polyval(c, vv)
    v_opt, j_opt = vv[jj.argmin()], jj.min()

    fig, ax = plt.subplots()
    ax.plot(vv, jj, color=AQUA, linewidth=2.0, zorder=2)
    ax.plot(df.velocity_m_s, j, linestyle="none", marker="o", color=AQUA,
            markerfacecolor=AQUA, markeredgecolor=SURFACE, zorder=3)

    span = jj.max() - jj.min()
    ax.set_ylim(jj.min() - 0.18 * span, jj.max() + 0.10 * span)

    ax.axvline(v_opt, color=INK_MUTED, linewidth=1.2, linestyle=(0, (3, 3)), zorder=1)
    annotate(ax, v_opt, j_opt, f"optimum ≈ {v_opt:.1f} m/s",
             dx=10, dy=10, color=INK, weight="600", size=10)

    ax.set_xlabel("Inlet velocity  (m/s)")
    ax.set_ylabel("Combined objective  J  (dimensionless)")
    title(ax, f"The best trade-off sits near {v_opt:.1f} m/s",
          f"J = R/R$_{{ref}}$ + {weight:g}·P/P$_{{ref}}$,  normalised at V = {ref.velocity_m_s:.0f} m/s")
    save(fig, FIGS / "05_combined_objective")
    return v_opt, j


def main():
    apply_style()
    df = load()

    print("Sweep summary")
    print(df.to_string(index=False,
                       float_format=lambda v: f"{v:8.4f}"))
    print()

    fig_thermal(df)
    fig_pressure(df)
    fig_pumping(df)
    fig_pareto(df)
    v_opt, j = fig_objective(df)

    df["objective_J"] = j
    df.to_csv(DATA / "derived_results.csv", index=False,
              float_format="%.5f")
    print(f"\n  wrote {DATA / 'derived_results.csv'}")
    print(f"\nOptimum inlet velocity (w = 1): {v_opt:.2f} m/s")


if __name__ == "__main__":
    main()
