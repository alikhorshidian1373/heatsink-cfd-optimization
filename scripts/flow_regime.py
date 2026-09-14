"""
Flow-regime, domain-sizing and near-wall (y+) analysis for the heat-sink study.

Three questions a reviewer will ask about any external-flow CFD model, answered
with numbers rather than assertion:

  1. Is the computational domain big enough that its outer boundaries are not
     influencing the result?  ->  blockage ratio and clearance-to-height ratios
     against the usual external-aerodynamics guidelines.

  2. What flow regime is this, and does the turbulence model suit it?
     ->  channel Reynolds number across the whole velocity sweep.

  3. Is the near-wall mesh adequate for the chosen turbulence model?
     ->  the first-cell height required to hit a target y+, and the y+ that a
     given first-cell height actually produces.

Question 3 is history now: the published mesh carries 8 `last-ratio` prism
layers and a measured y+ of 0.20-2.67 (see README section 4).  The near-wall
block below is kept because it explains *why* a default inflation stack could
not be used here - a 15-layer stack at growth 1.2 is 2.3x thicker than the
half-gap it has to fit into - and because it shows what the wall spacing would
have had to be had the layers been left out.  It is sizing rationale, not a
description of the final mesh.

Run:  python scripts/flow_regime.py
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# --- air: Fluent's built-in constant-property "air" --------------------------
RHO = 1.225          # kg/m3
MU = 1.7894e-5       # Pa.s
NU = MU / RHO        # m2/s

# --- heat sink ---------------------------------------------------------------
L_FIN = 0.060        # m, fin length along the flow
H_FIN = 0.025        # m, fin height
T_FIN = 0.0015       # m, fin thickness
PITCH = 0.0062       # m, fin pitch
N_FIN = 10
GAP = PITCH - T_FIN  # m, open channel width between two fins
T_BASE = 0.003       # m
W_BASE = 0.060       # m
H_SINK = T_BASE + H_FIN   # m, total frontal height of the sink

# --- air domain (full width; the solved model is a symmetric half) -----------
L_UP, L_DOWN = 0.300, 0.600      # m, upstream and downstream of the sink
CLEAR_SIDE, CLEAR_TOP = 0.100, 0.100   # m
DOM_X = L_UP + L_FIN + L_DOWN
DOM_Y = W_BASE + 2 * CLEAR_SIDE
DOM_Z = H_SINK + CLEAR_TOP

VELOCITIES = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])


# ---------------------------------------------------------------- domain -----
def domain_metrics():
    """Blockage ratio and clearance ratios against external-flow guidelines.

    The reference length for the clearance ratios is the frontal height of the
    obstacle, H, which is the convention used in the bluff-body and building-
    aerodynamics literature (e.g. the COST 732 / AIJ guidance: >= 5H upstream,
    >= 10-15H downstream, blockage below a few per cent).
    """
    frontal = W_BASE * H_SINK
    cross = DOM_Y * DOM_Z
    return {
        "domain_x_mm": DOM_X * 1e3,
        "domain_y_mm": DOM_Y * 1e3,
        "domain_z_mm": DOM_Z * 1e3,
        "frontal_area_mm2": frontal * 1e6,
        "cross_section_mm2": cross * 1e6,
        "blockage_ratio_pct": 100 * frontal / cross,
        "upstream_over_H": L_UP / H_SINK,
        "downstream_over_H": L_DOWN / H_SINK,
        "side_clearance_over_H": CLEAR_SIDE / H_SINK,
        "top_clearance_over_H": CLEAR_TOP / H_SINK,
    }


# ---------------------------------------------------------------- regime -----
def hydraulic_diameter():
    """D_h of one fin channel, wetted on three sides.

    The channel is open at the top: two fin side walls and the base plate are
    wetted, the top is not. Counting only the wetted perimeter is what makes
    this an unducted heat sink rather than a closed duct, and it is the same
    assumption that makes the flat-plate Nusselt correlation the right choice
    in validation.py.
    """
    area = GAP * H_FIN
    wetted = 2 * H_FIN + GAP
    return 4 * area / wetted


def regime_table():
    d_h = hydraulic_diameter()
    rows = []
    for v in VELOCITIES:
        re_dh = v * d_h / NU
        re_l = v * L_FIN / NU
        if re_dh < 2300:
            regime = "laminar"
        elif re_dh < 4000:
            regime = "transitional"
        else:
            regime = "turbulent"
        rows.append({"V_m_s": v, "Re_Dh": re_dh, "Re_L": re_l,
                     "channel_regime": regime})
    return pd.DataFrame(rows), d_h


# --------------------------------------------------------------- near wall ---
def cf_laminar(re_l):
    """Average skin-friction coefficient, laminar flat plate (Blasius)."""
    return 1.328 / np.sqrt(re_l)


def cf_turbulent(re_l):
    """Average skin-friction coefficient, turbulent flat plate (1/5 power law)."""
    return 0.074 * re_l ** -0.2


def u_tau(v, cf):
    """Friction velocity from the skin-friction coefficient."""
    tau_w = 0.5 * RHO * v ** 2 * cf
    return np.sqrt(tau_w / RHO)


def first_cell_for_target(v, y_plus_target, turbulent=False):
    """Wall-normal distance to the first cell centroid that gives a target y+."""
    re_l = v * L_FIN / NU
    cf = cf_turbulent(re_l) if turbulent else cf_laminar(re_l)
    return y_plus_target * NU / u_tau(v, cf)


def y_plus_for_spacing(v, y_centroid, turbulent=False):
    """The y+ that a given first-cell-centroid distance actually produces."""
    re_l = v * L_FIN / NU
    cf = cf_turbulent(re_l) if turbulent else cf_laminar(re_l)
    return u_tau(v, cf) * y_centroid / NU


def near_wall_table():
    """First-layer sizing targets, and the y+ implied by plausible tet spacings.

    The right-hand block answers a counterfactual: if no prism layers were
    used, the wall-adjacent cell would be a tetrahedron sized by the local
    body/proximity sizing, and its centroid would sit roughly one third of the
    cell height from the wall. The columns bracket the plausible range for that
    case. The mesh actually used does carry prism layers, and its y+ was
    measured rather than estimated - 0.20 to 2.67 across the sweep.
    """
    rows = []
    for v in VELOCITIES:
        row = {
            "V_m_s": v,
            "y_for_yplus1_mm": first_cell_for_target(v, 1.0) * 1e3,
            "y_for_yplus30_mm": first_cell_for_target(v, 30.0) * 1e3,
        }
        # y+ produced by a wall-adjacent tet of edge h, centroid at ~h/3
        for h_mm in (0.2, 0.5, 1.0):
            y_c = (h_mm / 3.0) * 1e-3
            row[f"yplus_at_h{h_mm:g}mm"] = y_plus_for_spacing(v, y_c)
        rows.append(row)
    return pd.DataFrame(rows)


def boundary_layer_table():
    """Laminar boundary-layer thickness at the fin trailing edge, vs the channel.

    Blasius, delta = 5 L / sqrt(Re_L). Comparing it with half the channel gap
    says whether the boundary layers growing off the two facing fin walls meet
    before the air leaves the array. Where they meet, the channel is behaving
    like a developing duct; where they do not, each fin wall is still a free
    flat plate. This is the physical test behind the correlation choice in
    validation.py, and it is not a single answer for the whole sweep.
    """
    rows = []
    for v in VELOCITIES:
        re_l = v * L_FIN / NU
        delta = 5 * L_FIN / np.sqrt(re_l)
        rows.append({
            "V_m_s": v,
            "delta_te_mm": delta * 1e3,
            "half_gap_mm": (GAP / 2) * 1e3,
            "delta_over_half_gap": delta / (GAP / 2),
            "boundary_layers_merge": bool(delta > GAP / 2),
        })
    return pd.DataFrame(rows)


def inflation_spec(y_plus_target=1.0, n_layers=15, growth=1.2):
    """Total inflation-layer thickness for the worst case (highest velocity).

    Sized at the top of the sweep because the boundary layer is thinnest and
    the required first layer smallest there; a layer stack that resolves 6 m/s
    also resolves every slower case.
    """
    v = VELOCITIES.max()
    y_centroid = first_cell_for_target(v, y_plus_target)
    first_height = 2 * y_centroid          # centroid sits mid-height of layer 1
    total = first_height * (growth ** n_layers - 1) / (growth - 1)
    return {
        "sized_at_V_m_s": v,
        "y_plus_target": y_plus_target,
        "first_layer_height_mm": first_height * 1e3,
        "n_layers": n_layers,
        "growth_rate": growth,
        "total_thickness_mm": total * 1e3,
        "channel_gap_mm": GAP * 1e3,
        "pct_of_half_gap": 100 * total / (GAP / 2),
    }


def feasible_inflation(y_plus_target=1.0, growth=1.2, max_fill=0.40):
    """Largest layer count whose stack still fits the fin channel.

    Two inflation stacks grow towards each other from the two facing fin walls,
    so each may occupy at most half the gap, and in practice should leave a
    core of tetrahedra between them: max_fill caps the stack at a fraction of
    the half-gap. This is the constraint that a default 15-layer setting
    violates in a 4.7 mm channel, which is why the mesher failed on it.
    """
    v = VELOCITIES.max()
    first_height = 2 * first_cell_for_target(v, y_plus_target)
    budget = max_fill * (GAP / 2)
    n = 0
    while True:
        total = first_height * (growth ** (n + 1) - 1) / (growth - 1)
        if total > budget:
            break
        n += 1
    total = first_height * (growth ** n - 1) / (growth - 1)
    return {
        "y_plus_target": y_plus_target,
        "growth_rate": growth,
        "max_fill_of_half_gap": max_fill,
        "n_layers_that_fit": n,
        "first_layer_height_mm": first_height * 1e3,
        "total_thickness_mm": total * 1e3,
        "pct_of_half_gap": 100 * total / (GAP / 2),
    }


# --------------------------------------------------------------------- main --
def main():
    dom = domain_metrics()
    print("Domain sizing")
    print(f"  air domain           {dom['domain_x_mm']:.0f} x {dom['domain_y_mm']:.0f}"
          f" x {dom['domain_z_mm']:.0f} mm (full width)")
    print(f"  obstacle height H    {H_SINK*1e3:.0f} mm")
    print(f"  blockage ratio       {dom['blockage_ratio_pct']:.2f} %"
          "   (guideline: a few per cent)")
    print(f"  upstream             {dom['upstream_over_H']:.1f} H"
          "   (guideline: >= 5 H)")
    print(f"  downstream           {dom['downstream_over_H']:.1f} H"
          "   (guideline: >= 10-15 H)")
    print(f"  side / top clearance {dom['side_clearance_over_H']:.1f} H"
          "   (guideline: >= 5 H)  <-- below guideline, see README")

    reg, d_h = regime_table()
    print(f"\nChannel Reynolds number   (D_h = {d_h*1e3:.2f} mm, "
          f"gap = {GAP*1e3:.1f} mm, three sides wetted)")
    print(reg.to_string(index=False, float_format=lambda v: f"{v:9.1f}"))

    nw = near_wall_table()
    print("\nNear-wall sizing")
    print("  left block: wall distance needed to hit a target y+")
    print("  right block: counterfactual y+ if the wall cell were a tetrahedron\n"
          "               of edge h -- the mesh used has 8 prism layers, y+ 0.20-2.67")
    print(nw.to_string(index=False, float_format=lambda v: f"{v:9.3f}"))

    bl = boundary_layer_table()
    print("\nBoundary-layer development along the fin (Blasius, laminar)")
    print(bl.to_string(index=False, float_format=lambda v: f"{v:9.2f}"))
    merge_v = bl.loc[bl.boundary_layers_merge, "V_m_s"]
    if len(merge_v):
        print(f"  boundary layers meet inside the channel up to "
              f"V = {merge_v.max():g} m/s; above that each fin wall is still a "
              "free flat plate at the trailing edge")

    infl = inflation_spec()
    print(f"\nInflation layer, default 15-layer stack, y+ = "
          f"{infl['y_plus_target']:g} at V = {infl['sized_at_V_m_s']:g} m/s")
    print(f"  first layer height   {infl['first_layer_height_mm']:.4f} mm")
    print(f"  {infl['n_layers']} layers at growth {infl['growth_rate']}"
          f"  ->  total {infl['total_thickness_mm']:.2f} mm")
    print(f"  channel half-gap     {infl['channel_gap_mm']/2:.2f} mm"
          f"   ({infl['pct_of_half_gap']:.0f} % of it) "
          "-> does not fit; the final mesh uses 8 last-ratio layers instead")

    fit = feasible_inflation()
    print(f"\n  a stack that does fit: {fit['n_layers_that_fit']} layers at "
          f"growth {fit['growth_rate']}, first layer "
          f"{fit['first_layer_height_mm']:.4f} mm")
    print(f"  -> total {fit['total_thickness_mm']:.2f} mm = "
          f"{fit['pct_of_half_gap']:.0f} % of the half-gap, leaving a "
          "tetrahedral core between the two stacks")

    DATA.mkdir(exist_ok=True)
    pd.DataFrame([dom]).to_csv(DATA / "domain_metrics.csv", index=False,
                               float_format="%.4f")
    reg.to_csv(DATA / "flow_regime.csv", index=False, float_format="%.2f")
    nw.to_csv(DATA / "near_wall_sizing.csv", index=False, float_format="%.5f")
    bl.to_csv(DATA / "boundary_layer.csv", index=False, float_format="%.4f")
    pd.DataFrame([infl, ]).to_csv(DATA / "inflation_spec.csv", index=False,
                                  float_format="%.5f")
    pd.DataFrame([fit]).to_csv(DATA / "inflation_feasible.csv", index=False,
                               float_format="%.5f")
    print(f"\nwrote 6 CSVs to {DATA}")


if __name__ == "__main__":
    main()
