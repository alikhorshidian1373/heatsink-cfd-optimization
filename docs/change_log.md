# What was wrong, how it was found, and what it cost

This file is the audit trail behind §7 of the README. It is kept because the
corrections were large enough that a reader is entitled to see how they were
arrived at, and because two of the three faults were found by questioning an
assumption rather than by running anything.

---

## Fault 1 — a 0.05° geometry defect that had blocked meshing for weeks

### Symptom

Every attempt to put a proper near-wall mesh on the assembly failed, in ways
that looked unrelated to each other:

- ANSYS Meshing reported *"the mesh generation did not complete"*, with the
  enclosure body producing 0 elements.
- SpaceClaim / DesignModeler reported **Failed Share Topology**.
- Fluent Meshing's watertight workflow found **7 659 self-intersecting faces**,
  all clustered at one repeating coordinate, `(0, 0.033100, 0.003000)`.
- Fluent reported **6 regions** where there should have been 2.
- `Combine` on the fin and base bodies did nothing at all when clicked.

Each of these was chased as its own problem. Three separate diagnoses — an
octree patch-independent method, lost scoping attachments, a MultiZone method
scoped to the wrong body — were tried and each was refuted by the next run.

### Diagnosis

The defect coordinate `z = 0.003000` is exactly the top of the 3 mm base, i.e.
the fin root. That pointed at the fin-to-base contact rather than at any mesher
setting. Measuring between the bodies in Discovery gave:

```
Minimum distance between objects   28.0035 mm
Angle between objects              0.05°
```

The five fin solids were tilted by **0.05°–0.07°** relative to the base. Over a
60 mm fin that opens a wedge-shaped gap running from 0 at one end to about
**73 µm** at the other. Every symptom follows from it:

- two faces 73 µm apart are not coincident, so share topology has nothing to
  share;
- the bodies touch along a line rather than a face, so `Combine` finds zero
  intersection volume and does nothing;
- the sliver volume between them is a separate region, which is why Fluent
  counted 6;
- the surface mesh bridges the sliver inconsistently, producing
  self-intersections at the fin roots.

`Force Share` made it worse — 13 380 self-intersections and
`Mesh topology corrupted (v-m-v-w)` — and was undone.

### Fix

The fins were deleted and rebuilt as **pulled features on the base body** rather
than as separate solids, so no contact exists to be misaligned. The enclosure
was regenerated around the single body. Fluent then reported:

```
2 fluid/solid regions,  0 voids
0 marked faces
maximum skewness 0.46
surface mesh complete in 5 s
```

### What it cost

Weeks of meshing attempts, and — until it was fixed — a model that could not
carry an inflation layer at all, which is Fault 2.

---

## Fault 2 — no prism layers, and therefore an unresolved boundary layer

### Symptom

The first version's mesh was tetrahedral with **no inflation layers**. y+ was
never measured; the README of that version estimated it at 0.5–9 from
flat-plate skin friction, which spans the viscous sublayer, the buffer layer and
the log layer. k-ω SST had been chosen specifically for its ability to integrate
to the wall, and nothing verified that it was being allowed to.

There was also a second, independent reason the layers had failed: a default
15-layer stack at growth ratio 1.2 sized for y+ = 1 at 6 m/s is **5.39 mm**
thick, and two such stacks have to fit inside a **4.7 mm** channel. The stack is
2.3× thicker than the half-gap available to it. `python scripts/flow_regime.py`
reproduces that calculation.

### Fix

With the geometry repaired, a poly-hexcore mesh was built with **8 prism layers**
using the `last-ratio` offset method (transition ratio 0.272), which sizes the
stack to the local gap instead of growing blindly:

| | |
|---|---|
| cells | 204 575 |
| min orthogonal quality | 0.195 |
| measured y+ | 0.20–2.67 across the whole sweep |

The full 1–6 m/s sweep was re-run on it.

### What it cost

| V (m/s) | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| R_th, no prism layers (K/W) | 1.794 | 1.155 | 0.880 | 0.732 | 0.650 | 0.596 |
| R_th, 8 prism layers (K/W) | 2.050 | 1.460 | 1.177 | 1.000 | 0.878 | 0.787 |
| change in R_th | +14 % | +26 % | +34 % | +37 % | +35 % | +32 % |
| change in Δp | −14 % | −18 % | −20 % | −21 % | −20 % | −18 % |

Both signs are the ones an under-resolved boundary layer predicts: too much
heat leaving the wall, too little shear at it. The velocity exponent moved from
**−0.62 to −0.53** and the optimum velocity from **4.5 to 5.0 m/s** at the same
weighting.

The independent confirmation is the analytical check. Against the same 1-D
fin-array model, the unresolved sweep was **8–28 % off with a consistent
one-sided bias**; the resolved sweep is **within ±5 % with a crossover near
4 m/s** that the boundary-layer merging criterion predicts. A correction that
moves a result *towards* an independent model is worth more than one that
merely moves it.

---

## Fault 3 — nothing carried an uncertainty

The first version reported six numbers from one mesh with one turbulence model.
Three things were added:

| Added | Result |
|---|---|
| Three-mesh grid convergence (Celik et al. 2008) | R_th uncertainty **1.75 %**; Δp not in the asymptotic range, bounded at **4.2 %** |
| Measured y+ across the sweep | 0.20–2.67, inside the viscous sublayer everywhere |
| Transition SST vs k-ω SST at 1 and 6 m/s | model-form sensitivity **≤ 2.2 %** |

The GCI implementation is self-tested against the worked example in the Celik
paper (`python scripts/gci.py --selftest`), because a GCI number nobody has
checked against a known answer is arithmetic with a citation attached.

A reproducibility check was also added: the 6 m/s case was re-solved from the
saved case file and returned **315.73992 K** against **315.7402 K**, and
**0.617657 Pa** against **0.617963 Pa** — 0.0003 K and 0.05 %.

---

## Corrections to the first version's text

| First version said | Correct |
|---|---|
| "Domain 960 × 260 × 128 mm, solved as a symmetric half" — ambiguous | The full extent is 960 × 260 × 128 mm; the domain actually solved is the 960 × **130** × 128 mm half |
| "10 fins" without qualification | 10 on the full sink, **5 in the model**; the symmetry plane falls in the central gap and cuts no fin |
| Optimum "2.74 m/s", later revised to "≈ 3 m/s" | 5.05 m/s at λ = 2.28 K/W². The earlier figures came from a different objective function *and* from the unresolved mesh, so neither is a like-for-like comparison |
| "y+ probably between 0.5 and 9" | Measured: 0.20–2.67 |

---

## What was tried and did not work

Recorded so nobody repeats it.

| Attempt | Outcome |
|---|---|
| Patch Independent method on the enclosure | enclosure went from 0 to 12 388 elements but never completed |
| Fluent Meshing `Apply Share Topology` | 0 pairs shared |
| DesignModeler Join-and-Intersect | no change |
| Discovery `Share` | 31 faces / 76 edges shared, counter still 0 |
| Discovery `Force Share` | **made it worse**: 13 380 self-intersections, corrupted topology, enclosure sharing with itself |
| `Combine` on fin and base bodies | no effect — the bodies barely touch |
| Selecting bodies from the tree in Meshing Details fields | not possible; Named Selection scoping is required |

Two workflow traps were also hit in the solver and are worth flagging:

1. Clicking the parent **Add Boundary Layers** task in the watertight workflow
   creates a *second* prism object that collides with the existing one
   (`Scoped Prism last-ratio_2 is trying to reset settings applied on zones`).
   Edit `last-ratio_1` directly instead.
2. **Replace Mesh** silently deletes solver-side zone separations, so the heated
   patch and its heat flux vanish and the solid equilibrates at the inlet
   temperature. A run that returns T_max = 300.000 K exactly is that, not a
   converged result.
