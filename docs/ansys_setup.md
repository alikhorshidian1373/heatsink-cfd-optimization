# ANSYS setup — everything needed to reproduce the runs

ANSYS Fluent 2026 R1, Student licence (512 000 cell cap), 4 parallel processes,
double precision.

---

## 1 · Geometry

Single solid body: a 60 × 60 × 3 mm base with **10 fins pulled from the base
body itself** (25 mm tall, 1.5 mm thick, 6.2 mm pitch). The fins must not be
separate solids assembled onto the base — see
[`change_log.md`](change_log.md), Fault 1.

An enclosure is built around the sink, then cut on the fin-array centre plane.
The computed half spans:

| | from | to |
|---|---:|---:|
| X (streamwise) | −0.300 m (inlet) | +0.660 m (outlet) |
| Y (spanwise) | −0.070 m (far side) | +0.060 m (symmetry plane) |
| Z (vertical) | 0.000 m (floor) | +0.128 m (top) |

The sink occupies X 0–0.060, Y 0.030–0.060, Z 0–0.028. The heated patch is
X 0.010–0.050, Y 0.040–0.060 on the base underside.

Share topology must report **Not-shared 0** before export.

---

## 2 · Fluent Meshing — watertight geometry workflow

| Task | Setting |
|---|---|
| Import Geometry | Use Body Labels = Yes |
| Add Local Sizings | none |
| Generate the Surface Mesh | Cells Per Gap 1.78 · growth 1.2 · min/max size per mesh, see below |
| Describe Geometry | geometry has both fluid and solid regions, share topology = Yes |
| Update Boundaries | see below |
| Create / Update Regions | 1 fluid region, 1 solid region, 0 voids |
| Add Boundary Layers | `last-ratio`, **8 layers**, transition ratio 0.272, first height per mesh, see below |
| Generate the Volume Mesh | poly-hexcore · Buffer 1 · Peel 1 · min/max cell length per mesh, see below |

A correct geometry reports `2 fluid/solid regions, 0 voids` and `0 marked
faces` at Describe Geometry, and completes the surface mesh in a few seconds.

### Boundary zones

The enclosure arrives as one wall zone. Split it before assigning types:

```
Boundary → Zone → Separate ,  Angle = 40
```

which gives six zones. The two smallest (equal face counts) are the inlet and
outlet; verify by area-weighted average of X-Coordinate — the inlet reads
−0.300.

| Zone | Type |
|---|---|
| `enclosure-enclosure:1:23` | velocity-inlet |
| `enclosure-enclosure:1:24` | pressure-outlet |
| `enclosure-enclosure:1`, `:1:21`, `:1:22`, `:1:25` | symmetry |
| `base-hs:1`, `base-hs-enclosure-enclosure` | wall |

### Mesh statistics

| Mesh | Cells | Min orthogonal quality | y+ max at 6 m/s |
|---|---:|---:|---:|
| fine | 204 575 | 0.195 | 2.67 |
| medium | 122 656 | 0.197 | 3.26 |
| coarse | 76 186 | 0.195 | 4.06 |

The three meshes come from the same workflow with the sizing scaled between
them. The coarse settings are recorded exactly:

| | coarse mesh |
|---|---|
| surface min / max size | 0.002028 m / 0.04056 m |
| volume min / max cell length | 0.002028 m / 0.032448 m |
| prism first height | 1.2675 × 10⁻⁴ m |

The fine and medium meshes use the same parameters scaled down by the
refinement ratio; **their exact sizing values were not transcribed at build
time** and should be read back from `heatsink-phc-8bl.msh.h5` and
`phc-medium.msh.h5` if the meshes need to be rebuilt from scratch rather than
re-used. This is a gap in the record, and it is noted rather than filled in
from memory.

---

## 3 · Solver setup

| | |
|---|---|
| Solver | pressure-based, steady, double precision |
| Energy | on |
| Viscous | **SST k-ω**; Low-Re Corrections **off**; near-wall treatment `correlation`; Production Limiter **on**; Curvature and Corner Flow Correction **off**; model constants at defaults |
| Materials | fluid `air` (constant properties); solid: Fluent's built-in `aluminum`, k = 202.4 W/m·K (pure aluminium, not a 6061 alloy) |
| Cell zones | `enclosure-enclosure` → air, `base-hs` → aluminium |
| Pressure-velocity coupling | SIMPLE |
| Gradient | least-squares cell based, **warped-face gradient correction on** |
| Discretisation | pressure second order; momentum, k, ω, energy second-order upwind |

### Boundary conditions

| Zone | Setting |
|---|---|
| velocity-inlet | magnitude normal to boundary, 1–6 m/s · 300 K · turbulence: **intensity 5 %, hydraulic diameter 0.129 m** |
| pressure-outlet | 0 Pa gauge, backflow 300 K |
| symmetry zones | symmetry |
| `base-hs-enclosure-enclosure` | coupled wall (created automatically with its shadow) |
| `base-hs:1.1` | symmetry (the solid's cut face on the centre plane) |
| `base-hs:1.2` | wall, **heat flux** (value below) |
| `base-hs:1` | wall, adiabatic (the rest of the base underside) |

### Creating the heated patch in the solver

The patch is carved out of the base underside after the mesh is read, so it does
not depend on the CAD:

1. `Domain → Zones → Separate → Faces`, option **Angle**, angle 40, zone
   `base-hs:1`. This splits the solid's non-wetted boundary into the base
   underside (~0.0018 m²) and the symmetry-plane cut face (~0.00018 m²).
   Set the small one to type `symmetry`.
2. Create a hexahedral cell register named `patch` covering
   X 0.010–0.050, Y 0.040–0.060, Z below the base.
3. `Domain → Zones → Separate → Faces`, option **Mark**, register `patch`, zone
   `base-hs:1`. The new zone `base-hs:1.2` is the heated patch.
4. Measure its area (`Reports → Surface Integrals → Area`) and set the flux to
   **10 W ÷ area**, so exactly 10 W enters the half-model whatever the mesh
   resolves the patch boundary to:

| Mesh | Patch area (m²) | Applied flux (W/m²) | Power in (W) |
|---|---:|---:|---:|
| fine | 0.00080238 | 12 463 | 10.000044 |
| medium | 0.00079721 | 12 543.8 | 10.000043 |
| coarse | 0.00081152 | 12 322.54 | 9.999988 |

The fine-mesh figure is Fluent's own flux report; the other two are the product
of the columns beside them. The fine-mesh area is quoted to the precision it was
read at, so its product differs from the reported power in the seventh digit.

---

## 4 · Running

Hybrid initialisation, then 600 iterations at 1–5 m/s and 1 500 at 6 m/s. The
residual convergence criteria stop the run far too early on this case — a
solution that "converges" at ~450 iterations is still drifting in T_max — so
**uncheck every box in the Check Convergence column of the Residual Monitors**
and let the reports flatten instead. Convergence was judged on the standard deviation of
T_max and Δp over the last 30 iterations (recorded in
`data/velocity_sweep.csv`).

The four-equation Transition SST runs need **2 000–3 000** iterations.

### Report definitions used

| Name | Definition |
|---|---|
| `tmax-solid` | volume maximum of static temperature over cell zone `base-hs` |
| `p-inlet` | area-weighted average of static pressure on the velocity inlet |
| `yplus-avg`, `yplus-max` | area-weighted average / facet maximum of wall y+ on `base-hs-enclosure-enclosure` |

### Transition SST comparison

Select **Transition SST (4 eqn)** in the Viscous panel, then **uncheck
Production Kato-Launder**, which Fluent switches on automatically for that model
— leaving it on would change two things at once and make the comparison
uninterpretable. Set Intermittency and Momentum Thickness Re to second-order
upwind in `Solution → Methods`.

---

## 5 · Traps

- **Replace Mesh deletes solver-side zone separations.** The heated patch and
  its flux disappear and the solid equilibrates at 300 K. A run returning
  T_max = 300.000 K exactly has lost its heat source, not converged.
- **Clicking the parent "Add Boundary Layers" task** creates a second scoped
  prism object that collides with the first
  (`Scoped Prism last-ratio_2 is trying to reset settings applied on zones`).
  Edit `last-ratio_1` directly.
- **Update Boundaries remembers the previous mesh's zone names.** Run the
  angle separation first, then refresh the task list by clicking away and back.
- **Camera TUI order matters**: set `target` and `position` before `up-vector`,
  or Fluent rejects the up-vector as collinear with the line of sight.

---

## 6 · Files produced

```
phc-V1..V6.cas.h5 / .dat.h5     fine-mesh sweep, k-omega SST
phc-medium-V6.cas.h5 / .dat.h5  medium mesh at 6 m/s
phc-coarse-V6.cas.h5 / .dat.h5  coarse mesh at 6 m/s
phc-V1-transSST.cas.h5          Transition SST at 1 m/s
phc-V6-transSST.cas.h5          Transition SST at 6 m/s
```

These are not in the repository — they total several gigabytes — but every
number extracted from them is in `data/`.
