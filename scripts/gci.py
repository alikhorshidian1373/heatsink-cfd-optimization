"""
Grid Convergence Index for the three-mesh study, following

    Celik, Ghia, Roache, Freitas, Coleman & Raad (2008),
    "Procedure for Estimation and Reporting of Uncertainty Due to
     Discretization in CFD Applications", J. Fluids Eng. 130(7).

Run:  python scripts/gci.py
      python scripts/gci.py --selftest      reproduces the paper's Table 1

The self-test matters: a GCI number that nobody has checked against a known
answer is just arithmetic with a citation attached.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def apparent_order(e21, e32, r21, r32, tol=1e-12, itmax=200):
    """Solve Celik eq. (3) for p by fixed-point iteration."""
    s = np.sign(e32 / e21)
    p = abs(np.log(abs(e32 / e21))) / np.log(r21)
    for _ in range(itmax):
        q = np.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
        if abs(p_new - p) < tol:
            return p_new, s
        p = p_new
    return p, s


def gci(f1, f2, f3, n1, n2, n3, dim=3):
    """f1 = finest. Returns a dict of the Celik quantities."""
    r21 = (n1 / n2) ** (1 / dim)
    r32 = (n2 / n3) ** (1 / dim)
    e21, e32 = f2 - f1, f3 - f2
    ratio = e32 / e21
    out = {"r21": r21, "r32": r32, "eps21": e21, "eps32": e32, "ratio": ratio,
           "e_a21": abs((f1 - f2) / f1)}
    if ratio > 1:                                   # monotonic convergence
        p, _ = apparent_order(e21, e32, r21, r32)
        f_ext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
        out.update(regime="monotonic convergence", p=p, f_ext=f_ext,
                   e_ext21=abs((f_ext - f1) / f_ext),
                   gci_fine=1.25 * out["e_a21"] / (r21 ** p - 1))
    else:                                           # not in the asymptotic range
        out.update(regime="non-asymptotic", p=np.nan, f_ext=np.nan,
                   e_ext21=np.nan,
                   gci_fine_p1=1.25 * out["e_a21"] / (r21 - 1),
                   gci_fine_p2=1.25 * out["e_a21"] / (r21 ** 2 - 1))
    return out


def selftest():
    """Celik et al. (2008), Table 1 - 2-D turbulent flow over a backward-facing step."""
    r = gci(6.063, 5.972, 5.863, 18000, 8000, 4500, dim=2)
    print("self-test against Celik et al. (2008) Table 1")
    print(f"  apparent order p   = {r['p']:.3f}   (paper: 1.53)")
    print(f"  extrapolated value = {r['f_ext']:.4f}  (paper: 6.1685)")
    print(f"  GCI_fine           = {r['gci_fine']*100:.3f} %  (paper: 2.17 %)")
    ok = abs(r["p"] - 1.53) < 0.02 and abs(r["gci_fine"] * 100 - 2.17) < 0.05
    print("  PASS" if ok else "  FAIL")
    return ok


def report(name, f1, f2, f3, n, unit):
    r = gci(f1, f2, f3, *n)
    print(f"\n=== {name} ===")
    print(f"  fine {f1:.6g}   medium {f2:.6g}   coarse {f3:.6g}  {unit}")
    print(f"  r21 = {r['r21']:.4f}   r32 = {r['r32']:.4f}")
    print(f"  eps32/eps21 = {r['ratio']:.4f}  ->  {r['regime']}")
    print(f"  e_a21 (fine vs medium) = {r['e_a21']*100:.3f} %")
    if r["regime"] == "monotonic convergence":
        print(f"  apparent order p = {r['p']:.3f}")
        print(f"  Richardson extrapolation = {r['f_ext']:.6g} {unit}"
              f"   (e_ext21 = {r['e_ext21']*100:.3f} %)")
        print(f"  GCI_fine = {r['gci_fine']*100:.2f} %"
              f"  ->  {f1:.5g} +/- {r['gci_fine']*abs(f1):.5g} {unit}")
    else:
        print("  successive differences GROW with refinement, so the three-grid")
        print("  observed order cannot be established. Reported conservatively:")
        print(f"    GCI_fine (assumed p = 1) = {r['gci_fine_p1']*100:.2f} %")
        print(f"    GCI_fine (assumed p = 2) = {r['gci_fine_p2']*100:.2f} %")
    return r


def main():
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)

    df = pd.read_csv(ROOT / "data" / "grid_convergence.csv")
    df = df.set_index("mesh").loc[["fine", "medium", "coarse"]]
    n = tuple(df.cells)
    print(f"meshes: fine {n[0]:,}  medium {n[1]:,}  coarse {n[2]:,} cells"
          f"   (all at 6 m/s)")
    report("temperature rise  dT  (= 20 x R_th)", *df.dT_K, n=n, unit="K")
    report("pressure drop  dp", *df.dP_Pa, n=n, unit="Pa")
    print("\nnote: the prism first-layer height was set per mesh to keep y+ sensible")
    print("      rather than scaled by r, so near-wall resolution does not refine")
    print("      systematically (y+max 2.67 / 3.26 / 4.06). dp is dominated by wall")
    print("      shear and inherits that; dT is not, and converges cleanly.")


if __name__ == "__main__":
    main()
