# =============================================================================
#  Heat-sink CFD study - figure generation
#  Run top-to-bottom in Google Colab. No uploads needed: all data is inline.
#  Outputs: PNG (200 dpi) + SVG into ./figures/, then zips them for download.
# =============================================================================

import os, json, zipfile
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

OUT = "figures"
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------
# 1. DATA  (fine mesh, 204 575 cells, 8 prism layers, k-omega SST)
# ----------------------------------------------------------------------------
V      = np.array([1, 2, 3, 4, 5, 6], float)              # inlet velocity [m/s]
TMAX   = np.array([340.9997, 329.1974, 323.5351,
                   320.0085, 317.5530, 315.7402])          # max solid temp [K]
DP     = np.array([0.0247929, 0.0805411, 0.1668621,
                   0.2848628, 0.4348673, 0.6179630])       # inlet gauge pressure [Pa]
YPAVG  = np.array([0.195669, 0.317226, 0.431602,
                   0.542471, 0.651344, 0.758252])
YPMAX  = np.array([0.710580, 1.184328, 1.595034,
                   1.989020, 2.342731, 2.668583])

T_AMB   = 300.0          # inlet air temperature [K]
Q_TOTAL = 20.0           # heat load of the full device [W] (10 W in the half model)
A_DUCT  = 0.03328        # full duct cross-section 0.260 x 0.128 [m^2]

DT    = TMAX - T_AMB                 # temperature rise [K]
R_TH  = DT / Q_TOTAL                 # thermal resistance [K/W]
QFLOW = V * A_DUCT                   # volumetric flow [m^3/s]
PPUMP = DP * QFLOW                   # ideal pumping power [W]

# numerical uncertainty, from the three-mesh GCI study at 6 m/s
GCI_R  = 0.0175          # 1.75 %  (observed order p = 1.05, monotonic)
GCI_DP = 0.0423          # 4.23 %  (non-asymptotic -> formal order p = 2 assumed)

# grid convergence triplet (6 m/s)
GCI_N   = np.array([76186, 122656, 204575], float)       # coarse -> fine
GCI_DT  = np.array([15.8309, 15.78310, 15.73992])        # temperature rise [K]
GCI_H   = (0.960 * 0.130 * 0.128 / GCI_N) ** (1 / 3) * 1e3   # mean cell size [mm]
DT_EXT  = 15.5197                                        # Richardson extrapolation
P_OBS   = 1.051                                          # observed order

# turbulence-model sensitivity: k-omega SST vs Transition SST (4 eqn)
MODEL = [
    # label,                 SST,        Transition SST
    (r"$R_{th}$  @ 1 m/s",   2.04999,    2.01249),
    (r"$R_{th}$  @ 6 m/s",   0.786996,   0.779484),
    (r"$\Delta p$  @ 1 m/s", 0.0247929,  0.0253315),
    (r"$\Delta p$  @ 6 m/s", 0.6176570,  0.6121786),
]

LAMBDA_REF = 2.2803      # reference weighting in J = R_th + lambda * P_pump

# ----------------------------------------------------------------------------
# 2. POWER-LAW FITS
# ----------------------------------------------------------------------------
def powerfit(x, y):
    n, lnA = np.polyfit(np.log(x), np.log(y), 1)
    A = np.exp(lnA)
    err = np.max(np.abs(A * x ** n / y - 1)) * 100
    return float(A), float(n), float(err)

A_R, n_R, e_R = powerfit(V, R_TH)
A_P, n_P, e_P = powerfit(V, DP)
A_W, n_W, e_W = powerfit(V, PPUMP)

def v_opt(lam):
    return (-A_R * n_R / (lam * A_W * n_W)) ** (1.0 / (n_W - n_R))

V_OPT = v_opt(LAMBDA_REF)

print(f"R_th   = {A_R:.5f} V^{n_R:+.4f}   max deviation {e_R:.2f} %")
print(f"dp     = {A_P:.5f} V^{n_P:+.4f}   max deviation {e_P:.2f} %")
print(f"P_pump = {A_W:.6f} V^{n_W:+.4f}   max deviation {e_W:.2f} %")
print(f"V_opt(lambda={LAMBDA_REF}) = {V_OPT:.3f} m/s")

# ----------------------------------------------------------------------------
# 3. STYLE
# ----------------------------------------------------------------------------
C = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100",
         magenta="#e87ba4", green="#008300", violet="#4a3aa7", red="#e34948")
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURFACE, GRIDC  = "#ffffff", "#e8e7e3"

mpl.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 200,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "savefig.bbox": "tight", "savefig.pad_inches": 0.28,
    "font.family": "DejaVu Sans", "font.size": 10.5,
    "text.color": INK, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": INK3, "ytick.color": INK3,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.edgecolor": GRIDC, "axes.linewidth": 1.0,
    "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.9,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 9.5,
    "axes.titlesize": 12.5, "axes.titleweight": "semibold",
    "axes.titlepad": 14, "axes.labelpad": 8,
    "lines.linewidth": 2.0, "lines.markersize": 6.5,
    "xtick.major.size": 0, "ytick.major.size": 0,
})

def finish(ax, title, sub=None, xlabel=None, ylabel=None, xgrid=False):
    ax.set_title(title, loc="left", pad=32 if sub else 14)
    if sub:
        ax.annotate(sub, xy=(0, 1), xycoords="axes fraction",
                    xytext=(0, 9), textcoords="offset points",
                    fontsize=9.5, color=INK2, va="bottom", ha="left")
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.grid(axis="y", zorder=0)
    ax.grid(axis="x", visible=xgrid, zorder=0)
    ax.set_axisbelow(True)

def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT}/{name}.{ext}")
    plt.close(fig)
    print("  wrote", name)

VSMOOTH = np.linspace(0.85, 6.4, 300)

# ----------------------------------------------------------------------------
# fig00 - domain and geometry
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(11.2, 6.4))
gs  = fig.add_gridspec(2, 2, height_ratios=[1, 1.9], width_ratios=[1, 1.25],
                       hspace=0.18, wspace=0.06, top=0.94, bottom=0.02)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])

DOM_BLUE = "#f3f7fd"
L_UP, L_HS, L_DN, H_DOM, W_HALF = 300., 60., 600., 128., 130.
L_TOT = L_UP + L_HS + L_DN

# --- (a) streamwise plane ----------------------------------------------------
axA.add_patch(Rectangle((0, 0), L_TOT, H_DOM, fc=DOM_BLUE, ec=GRIDC, lw=1.2, zorder=1))
axA.add_patch(Rectangle((L_UP, 0), L_HS, 3, fc=C["violet"], ec="none", zorder=3))
for i in range(6):
    axA.add_patch(Rectangle((L_UP + 3 + i * 10.5, 3), 5, 25,
                            fc=C["violet"], ec="none", zorder=3))
axA.add_patch(Rectangle((L_UP + 10, -3.4), 40, 3.4, fc=C["red"], ec="none", zorder=4))
axA.plot([0, 0], [0, H_DOM], color=C["blue"], lw=4, solid_capstyle="butt", zorder=5)
axA.plot([L_TOT, L_TOT], [0, H_DOM], color=C["orange"], lw=4,
         solid_capstyle="butt", zorder=5)

for x0, x1, lab in [(0, L_UP, "300 mm   (10.7 H)"),
                    (L_UP + L_HS, L_TOT, "600 mm   (21.4 H)")]:
    axA.annotate("", xy=(x0, 104), xytext=(x1, 104),
                 arrowprops=dict(arrowstyle="<->", color=INK3, lw=1.1))
    axA.text((x0 + x1) / 2, 110, lab, ha="center", fontsize=9.5, color=INK2)

axA.text(0, -16, "velocity inlet\n1 - 6 m/s,  300 K", ha="left", va="top",
         fontsize=9.5, color=C["blue"])
axA.text(L_TOT, -16, "pressure outlet\n0 Pa gauge", ha="right", va="top",
         fontsize=9.5, color=C["orange"])
axA.annotate("heat sink", xy=(L_UP + 30, 28), xytext=(L_UP + 165, 74),
             fontsize=9.5, color=C["violet"], va="center",
             arrowprops=dict(arrowstyle="-", color=C["violet"], lw=1, alpha=0.6))
axA.annotate("40 x 40 mm heated patch,  20 W", xy=(L_UP + 30, -3.4),
             xytext=(L_UP + 175, -30), fontsize=9.5, color=C["red"], va="center",
             arrowprops=dict(arrowstyle="-", color=C["red"], lw=1, alpha=0.6))
axA.set_xlim(-30, L_TOT + 30); axA.set_ylim(-54, 142)
axA.set_aspect("equal"); axA.axis("off")
axA.set_title("(a)   streamwise plane,  960 mm long", loc="left", pad=8, fontsize=11)

# --- (b) front view of the computed half -------------------------------------
axB.add_patch(Rectangle((0, 0), W_HALF, H_DOM, fc=DOM_BLUE, ec=GRIDC, lw=1.2, zorder=1))
axB.add_patch(Rectangle((0, 0), 30, 3, fc=C["violet"], ec="none", zorder=3))
for i in range(5):                                  # fin centres 3.1 + i*6.2 mm
    axB.add_patch(Rectangle((3.1 + i * 6.2 - 0.75, 3), 1.5, 25,
                            fc=C["violet"], ec="none", zorder=3))
axB.plot([0, 0], [0, H_DOM], color=C["aqua"], lw=4, solid_capstyle="butt", zorder=5)
axB.text(-6, H_DOM / 2, "symmetry plane", color=C["aqua"], fontsize=9.5,
         rotation=90, va="center", ha="center")

axB.annotate("", xy=(30, 104), xytext=(W_HALF, 104),
             arrowprops=dict(arrowstyle="<->", color=INK3, lw=1.1))
axB.text((30 + W_HALF) / 2, 109, "100 mm  (3.6 H)", ha="center", fontsize=9.5, color=INK2)
axB.annotate("", xy=(19, 28), xytext=(19, H_DOM),
             arrowprops=dict(arrowstyle="<->", color=INK3, lw=1.1))
axB.text(23, 76, "100 mm  (3.6 H)", fontsize=9.5, color=INK2, va="center")
axB.annotate("5 fins modelled\n(10 on the full sink)", xy=(16, 16),
             xytext=(52, 28), fontsize=9.5, color=C["violet"], va="center",
             arrowprops=dict(arrowstyle="-", color=C["violet"], lw=1, alpha=0.6))
axB.set_xlim(-18, W_HALF + 8); axB.set_ylim(-14, 142)
axB.set_aspect("equal"); axB.axis("off")
axB.set_title("(b)   half-domain cross-section", loc="left", pad=8, fontsize=11)

# --- (c) setup summary -------------------------------------------------------
axC.axis("off")
ROWS = [
    ("Heat sink",  "60 x 60 x 3 mm aluminium base  (k = 202.4 W/m-K)"),
    ("",           "10 fins  -  25 mm tall, 1.5 mm thick, 6.2 mm pitch"),
    ("Heat load",  "20 W over a central 40 x 40 mm patch  (12.5 kW/m²)"),
    ("Fluid",      "air, constant properties, 300 K at the inlet"),
    ("Domain",     "960 x 260 x 128 mm, solved as the symmetric half"),
    ("Blockage",   "5.05 %   -   sink frontal height H = 28 mm"),
    ("Physics",    "steady RANS, conjugate heat transfer, k-omega SST"),
    ("Mesh",       "204 575 cells, poly-hexcore, 8 prism layers"),
    ("",           "min orthogonal quality 0.195,  y+ = 0.20 - 2.67"),
]
y = 0.97
for k, v in ROWS:
    if k:
        axC.text(0.0, y, k, fontsize=9.5, color=INK2, va="top", fontweight="semibold")
    axC.text(0.30, y, v, fontsize=9.5, color=INK, va="top")
    y -= 0.108
axC.set_xlim(0, 1); axC.set_ylim(0, 1)
axC.set_title("(c)   case setup", loc="left", pad=8, fontsize=11)

fig.suptitle("Computational domain and case setup",
             x=0.005, y=1.045, ha="left", fontsize=13.5, fontweight="semibold")
fig.text(0.005, 1.002, "the sink is a 60 mm plate-fin extrusion in a 960 mm duct; "
                       "only half is computed, split on the centre plane",
         ha="left", fontsize=10, color=INK2)
save(fig, "fig00_domain_and_geometry")

# ----------------------------------------------------------------------------
# fig01 - thermal resistance
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
fit = A_R * VSMOOTH ** n_R
ax.fill_between(VSMOOTH, fit * (1 - GCI_R), fit * (1 + GCI_R),
                color=C["blue"], alpha=0.13, lw=0, zorder=1)
ax.plot(VSMOOTH, fit, color=C["blue"], lw=1.6, alpha=0.55, zorder=2)
ax.plot(V, R_TH, "o", color=C["blue"], mec=SURFACE, mew=1.6, zorder=3)
for v, r in [(V[0], R_TH[0]), (V[-1], R_TH[-1])]:
    ax.annotate(f"{r:.3f}", (v, r), textcoords="offset points",
                xytext=(0, 13), ha="center", fontsize=9.5, color=INK)
ax.text(4.3, 1.62, f"$R_{{th}} = {A_R:.3f}\\,V^{{{n_R:.3f}}}$\n"
                   f"fit within {e_R:.1f} % of every point",
        fontsize=10, color=INK2, va="top")
ax.text(2.35, 1.34, f"shaded: ±{GCI_R*100:.2f} % discretisation\nuncertainty (GCI)",
        fontsize=9, color=C["blue"], va="top")
finish(ax, "Thermal resistance falls off as a power law",
       "junction-to-ambient, 20 W total load, fine mesh (204 575 cells)",
       "inlet velocity  [m/s]", "$R_{th}$  [K/W]")
ax.set_xlim(0.7, 6.5); ax.set_ylim(0.6, 2.25)
save(fig, "fig01_thermal_resistance")

# ----------------------------------------------------------------------------
# fig02 - pressure drop
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
fit = A_P * VSMOOTH ** n_P
ax.fill_between(VSMOOTH, fit * (1 - GCI_DP), fit * (1 + GCI_DP),
                color=C["orange"], alpha=0.14, lw=0, zorder=1)
ax.plot(VSMOOTH, fit, color=C["orange"], lw=1.6, alpha=0.55, zorder=2)
ax.plot(V, DP, "o", color=C["orange"], mec=SURFACE, mew=1.6, zorder=3)
ax.annotate(f"{DP[-1]:.3f} Pa", (V[-1], DP[-1]), textcoords="offset points",
            xytext=(-6, 10), ha="right", fontsize=9.5, color=INK)
ax.text(1.05, 0.50, f"$\\Delta p = {A_P:.4f}\\,V^{{{n_P:.3f}}}$\n"
                    f"shaded: ±{GCI_DP*100:.2f} % (formal 2nd order;\n"
                    f"the triplet was not asymptotic for $\\Delta p$)",
        fontsize=10, color=INK2, va="top")
finish(ax, "Pressure drop rises as $V^{1.8}$, not $V^2$",
       "area-weighted static pressure at the inlet, outlet held at 0 Pa gauge",
       "inlet velocity  [m/s]", "$\\Delta p$  [Pa]")
ax.set_xlim(0.7, 6.5); ax.set_ylim(0, 0.72)
save(fig, "fig02_pressure_drop")

# ----------------------------------------------------------------------------
# fig03 - pumping power
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
fit = A_W * VSMOOTH ** n_W
ax.fill_between(VSMOOTH, fit * (1 - GCI_DP), fit * (1 + GCI_DP),
                color=C["aqua"], alpha=0.16, lw=0, zorder=1)
ax.plot(VSMOOTH, fit, color=C["aqua"], lw=1.6, alpha=0.6, zorder=2)
ax.plot(V, PPUMP, "o", color=C["aqua"], mec=SURFACE, mew=1.6, zorder=3)
for v, p in [(V[2], PPUMP[2]), (V[-1], PPUMP[-1])]:
    ax.annotate(f"{p*1000:.1f} mW", (v, p), textcoords="offset points",
                xytext=(-8, 6), ha="right", fontsize=9.5, color=INK)
ax.text(1.1, 0.104, f"$P_{{pump}} = \\Delta p\\,\\cdot\\,Q \\propto V^{{{n_W:.2f}}}$\n"
                    "doubling the flow costs ~7x the fan power,\n"
                    "and buys only a 31 % drop in $R_{th}$",
        fontsize=10, color=INK2, va="top")
finish(ax, "Fan power is the reason there is an optimum at all",
       "ideal pumping power $\\Delta p \\cdot Q$, half-domain flow scaled to the full duct",
       "inlet velocity  [m/s]", "$P_{pump}$  [W]")
ax.set_xlim(0.7, 6.5); ax.set_ylim(0, 0.135)
save(fig, "fig03_pumping_power")

# ----------------------------------------------------------------------------
# fig04 - trade-off front
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 4.6))
pw = A_W * VSMOOTH ** n_W
rt = A_R * VSMOOTH ** n_R
ax.plot(pw, rt, color=INK3, lw=1.6, alpha=0.5, zorder=2)
ax.plot(PPUMP, R_TH, "o", color=C["blue"], mec=SURFACE, mew=1.6, zorder=4)
for v, p, r in zip(V, PPUMP, R_TH):
    ax.annotate(f"{v:.0f} m/s", (p, r), textcoords="offset points",
                xytext=(9, 7), fontsize=9, color=INK2)
p_opt, r_opt = A_W * V_OPT ** n_W, A_R * V_OPT ** n_R
ax.plot([p_opt], [r_opt], "o", ms=11, color=C["red"], mec=SURFACE, mew=2, zorder=5)
ax.annotate(f"optimum  {V_OPT:.2f} m/s\n$\\lambda = {LAMBDA_REF}$ K/W²",
            (p_opt, r_opt), textcoords="offset points", xytext=(16, 24),
            ha="left", fontsize=9.5, color=C["red"])
finish(ax, "Every extra watt of cooling is bought with fan power",
       "each point is one converged solve; the optimum minimises "
       "$J = R_{th} + \\lambda\\,P_{pump}$",
       "pumping power  [W]", "$R_{th}$  [K/W]", xgrid=True)
ax.set_xlim(-0.006, 0.142); ax.set_ylim(0.6, 2.25)
save(fig, "fig04_tradeoff_front")

# ----------------------------------------------------------------------------
# fig05 - where the optimum sits, as a function of the weighting
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
lam = np.logspace(np.log10(0.15), np.log10(40), 400)
ax.plot(lam, v_opt(lam), color=C["violet"], zorder=3)
ax.axvline(LAMBDA_REF, color=GRIDC, lw=1.4, zorder=1)
ax.plot([LAMBDA_REF], [V_OPT], "o", ms=10, color=C["red"], mec=SURFACE, mew=2, zorder=4)
ax.annotate(f"reference weighting\n$\\lambda$ = {LAMBDA_REF}  ->  {V_OPT:.2f} m/s",
            (LAMBDA_REF, V_OPT), textcoords="offset points", xytext=(16, 14),
            fontsize=9.5, color=C["red"])
ax.axhspan(1, 6, color=C["blue"], alpha=0.07, lw=0, zorder=0)
ax.text(0.17, 6.15, "velocities actually simulated", fontsize=9, color=C["blue"])
ax.set_xscale("log")
ax.set_xticks([0.2, 0.5, 1, 2, 5, 10, 20, 40])
ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
finish(ax, "The optimum is a soft one - it moves with how you value fan power",
       "a tenfold change in $\\lambda$ shifts the best velocity by about a factor of two",
       "weighting  $\\lambda$  [K/W²]", "optimum velocity  [m/s]", xgrid=True)
ax.set_ylim(1.5, 11)
save(fig, "fig05_optimum_vs_lambda")

# ----------------------------------------------------------------------------
# fig06 - y+ verification
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
ax.axhspan(0, 5, color=C["aqua"], alpha=0.10, lw=0, zorder=0)
ax.axhline(5, color=C["aqua"], lw=1.4, ls=(0, (5, 4)), zorder=1)
ax.text(6.42, 5.28, "edge of the viscous sublayer,  $y^+ = 5$",
        fontsize=9, color=C["aqua"], ha="right")
ax.plot(V, YPMAX, "o-", color=C["orange"], mec=SURFACE, mew=1.6,
        label="maximum over the sink surface", zorder=3)
ax.plot(V, YPAVG, "o-", color=C["blue"], mec=SURFACE, mew=1.6,
        label="area-weighted average", zorder=3)
ax.annotate(f"{YPMAX[-1]:.2f}", (V[-1], YPMAX[-1]), textcoords="offset points",
            xytext=(0, 12), ha="center", fontsize=9.5, color=INK)
ax.annotate(f"{YPAVG[-1]:.2f}", (V[-1], YPAVG[-1]), textcoords="offset points",
            xytext=(0, -20), ha="center", fontsize=9.5, color=INK)
ax.legend(loc="upper left", bbox_to_anchor=(0.005, 0.97))
finish(ax, "The boundary layer is resolved at every velocity in the sweep",
       "8 prism layers on every wetted wall - no wall function is doing the work",
       "inlet velocity  [m/s]", "wall  $y^+$")
ax.set_xlim(0.7, 6.5); ax.set_ylim(0, 6.4)
save(fig, "fig06_yplus_verification")

# ----------------------------------------------------------------------------
# fig07 - grid convergence
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
hh = np.linspace(0, GCI_H[0] * 1.06, 200)
ax.plot(hh, DT_EXT + (GCI_DT[-1] - DT_EXT) * (hh / GCI_H[-1]) ** P_OBS,
        color=INK3, lw=1.5, ls=(0, (5, 4)), alpha=0.7, zorder=2,
        label=f"Richardson fit, observed order $p$ = {P_OBS:.2f}")
ax.plot(GCI_H, GCI_DT, "o", color=C["blue"], mec=SURFACE, mew=1.6, zorder=4,
        label="three systematically refined meshes")
ax.errorbar([0], [DT_EXT], yerr=[GCI_DT[-1] * GCI_R], fmt="o", ms=9,
            color=C["red"], mec=SURFACE, mew=1.8, ecolor=C["red"], elinewidth=2,
            capsize=5, zorder=5, label="extrapolated to zero cell size")
for h, d, n in zip(GCI_H, GCI_DT, GCI_N):
    ax.annotate(f"{n/1000:.0f}k", (h, d), textcoords="offset points",
                xytext=(0, 12), ha="center", fontsize=9, color=INK2)
ax.annotate(f"GCI$_{{fine}}$ = {GCI_R*100:.2f} %",
            (0, DT_EXT), textcoords="offset points", xytext=(14, -22),
            fontsize=9.5, color=C["red"])
ax.legend(loc="lower right")
finish(ax, "Grid convergence is monotonic and close to first order",
       "temperature rise at 6 m/s against mean cell size; Celik et al. (2008) procedure",
       "mean cell size  [mm]", "$\\Delta T$  [K]", xgrid=True)
ax.set_xlim(-0.35, 6.4); ax.set_ylim(15.10, 15.95)
save(fig, "fig07_grid_convergence")

# ----------------------------------------------------------------------------
# fig08 - turbulence-model sensitivity
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 4.0))
labels = [m[0] for m in MODEL]
delta = np.array([(m[2] / m[1] - 1) * 100 for m in MODEL])
ypos = np.arange(len(labels))[::-1]
ax.axvspan(-2, 2, color=GRIDC, alpha=0.55, lw=0, zorder=0)
for y, d in zip(ypos, delta):
    col = C["blue"] if d < 0 else C["red"]
    ax.barh(y, d, height=0.46, color=col, zorder=3)
    ax.annotate(f"{d:+.2f} %", (d, y), textcoords="offset points",
                xytext=(7 if d > 0 else -7, 0), ha="left" if d > 0 else "right",
                va="center", fontsize=10, color=INK)
ax.axvline(0, color=INK3, lw=1.2, zorder=4)
ax.set_yticks(ypos); ax.set_yticklabels(labels)
ax.text(2.55, -0.45, "±2 % band", fontsize=9, color=INK2, ha="center", va="center")
finish(ax, "Treating the flow as fully turbulent barely moves the answer",
       "Transition SST (4 eqn) relative to k-ω SST; mean intermittency on the "
       "sink surface was 0.03-0.04, i.e. the boundary layers stay laminar",
       "change when the transition model is used  [%]", None)
ax.grid(axis="y", visible=False); ax.grid(axis="x", zorder=0)
ax.set_xlim(-3.1, 3.1); ax.set_ylim(-0.6, 3.5)
save(fig, "fig08_model_sensitivity")

# ----------------------------------------------------------------------------
# fig09 - summary panel
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.6))
(a1, a2), (a3, a4) = axes

a1.fill_between(VSMOOTH, A_R * VSMOOTH ** n_R * (1 - GCI_R),
                A_R * VSMOOTH ** n_R * (1 + GCI_R), color=C["blue"], alpha=0.13, lw=0)
a1.plot(V, R_TH, "o-", color=C["blue"], mec=SURFACE, mew=1.5)
finish(a1, "Thermal resistance", None, "inlet velocity  [m/s]", "$R_{th}$  [K/W]")

a2.fill_between(VSMOOTH, A_P * VSMOOTH ** n_P * (1 - GCI_DP),
                A_P * VSMOOTH ** n_P * (1 + GCI_DP), color=C["orange"], alpha=0.14, lw=0)
a2.plot(V, DP, "o-", color=C["orange"], mec=SURFACE, mew=1.5)
finish(a2, "Pressure drop", None, "inlet velocity  [m/s]", "$\\Delta p$  [Pa]")

a3.plot(V, YPMAX, "o-", color=C["orange"], mec=SURFACE, mew=1.5, label="max")
a3.plot(V, YPAVG, "o-", color=C["blue"], mec=SURFACE, mew=1.5, label="average")
a3.axhline(5, color=C["aqua"], lw=1.3, ls=(0, (5, 4)))
a3.text(6.4, 5.25, "$y^+ = 5$", fontsize=9, color=C["aqua"], ha="right")
a3.legend(loc="upper left")
finish(a3, "Near-wall resolution", None, "inlet velocity  [m/s]", "wall  $y^+$")
a3.set_ylim(0, 6.4)

a4.plot(A_W * VSMOOTH ** n_W, A_R * VSMOOTH ** n_R, color=INK3, lw=1.5, alpha=0.5)
a4.plot(PPUMP, R_TH, "o", color=C["blue"], mec=SURFACE, mew=1.5)
a4.plot([p_opt], [r_opt], "o", ms=10, color=C["red"], mec=SURFACE, mew=2)
a4.annotate(f"{V_OPT:.2f} m/s", (p_opt, r_opt), textcoords="offset points",
            xytext=(-12, 16), ha="right", fontsize=9.5, color=C["red"])
finish(a4, "Where the two meet", None, "pumping power  [W]", "$R_{th}$  [K/W]", xgrid=True)

fig.suptitle("Air-cooled heat sink   -   verified velocity sweep, 1 to 6 m/s",
             x=0.005, y=1.035, ha="left", fontsize=13.5, fontweight="semibold")
fig.text(0.005, 0.995,
         f"k-ω SST with the boundary layer resolved  |  discretisation uncertainty "
         f"{GCI_R*100:.2f} % on $R_{{th}}$  |  turbulence-model sensitivity ≤ 2.2 %",
         ha="left", fontsize=10, color=INK2)
fig.tight_layout(h_pad=3.2, w_pad=3.6)
save(fig, "fig09_summary_panel")

# ----------------------------------------------------------------------------
# fig10 - analytical cross-check
# ----------------------------------------------------------------------------
# 1-D fin-array model: laminar flat plate Nu = 0.664 Re^0.5 Pr^(1/3),
# straight-fin efficiency, plus base conduction.  See scripts/validation.py.
R_FLAT = np.array([1.95715, 1.39662, 1.14826, 1.00019, 0.89913, 0.82452])
R_DUCT = 1.57687   # fully-developed parallel-plate channel, Nu = 7.54

fig, ax = plt.subplots(figsize=(7.2, 4.5))
ax.fill_between(VSMOOTH, A_R * VSMOOTH ** n_R * (1 - GCI_R),
                A_R * VSMOOTH ** n_R * (1 + GCI_R), color=C["blue"], alpha=0.13, lw=0)
ax.plot(V, R_TH, "o-", color=C["blue"], mec=SURFACE, mew=1.6, zorder=4,
        label="CFD, boundary layer resolved")
ax.plot(V, R_FLAT, "s--", color=C["orange"], mec=SURFACE, mew=1.6, ms=6, zorder=3,
        label="1-D model, laminar flat plate")
ax.axhline(R_DUCT, color=C["aqua"], lw=2, ls=(0, (5, 4)), zorder=2,
           label="1-D model, fully-developed duct")
ax.annotate("a ducted-flow assumption predicts\nno velocity dependence at all",
            (4.1, R_DUCT), textcoords="offset points", xytext=(0, 10),
            fontsize=9, color=C["aqua"])
ax.annotate("within 5 % of the CFD at every point,\nwith a crossover near 4 m/s", (1.15, 0.78),
            fontsize=9.5, color=INK2, ha="left")
ax.legend(loc="upper right")
finish(ax, "An independent 1-D model agrees to within 5 %",
       "the same comparison was 8-28 % off before the boundary layer was resolved",
       "inlet velocity  [m/s]", "$R_{th}$  [K/W]")
ax.set_xlim(0.7, 6.5); ax.set_ylim(0.6, 2.35)
save(fig, "fig10_analytical_validation")

# ----------------------------------------------------------------------------
# 4. ZIP FOR DOWNLOAD
# ----------------------------------------------------------------------------
zpath = "heatsink_figures.zip"
with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(os.listdir(OUT)):
        z.write(os.path.join(OUT, f), f)
print("\nzipped ->", zpath)

try:
    from google.colab import files
    files.download(zpath)
except Exception:
    print("(not on Colab - the figures are in ./figures/)")
