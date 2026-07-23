# Changelog

## v1.1

A maintenance release containing two correctness fixes and a substantial
speed-up of the sphere generation and the ODB post-processing.

**The published algorithm is unchanged.** Everything described in

> Y. Lim, M. Song, S. Ha, *Spatium: An open-source Abaqus/CAE plug-in for
> computational homogenization of particle-reinforced composites*,
> SoftwareX **31** (2025) 102219, doi:10.1016/j.softx.2025.102219

still applies to this version: random sequential placement of spheres with
radii drawn from `N(r_avg, r_std)`, the two-stage voxel + octree overlap check
of Fig. 2(b), the octree search radius `r + min_dist + L/8`, the `Tol_ = 0.4`
minimum-overlap-portion rule at the cell boundary, the dummy SFM3D4 surface
approach for periodic boundary conditions, and the `IVOL`-weighted volume
average of Eq. (11).

Only the implementation changed. The GUI, the input parameters, the generated
Abaqus model and the output file format are all identical to v1.0, so no
change is required in existing workflows — replace `Spatium_Kernel.py` and
restart Abaqus/CAE.

---

### Fixed

**Stage-1 of the overlap check was inactive at run time.**
The module-level `voxelGrid` object was declared as `None` and never
instantiated, and `is_overlapping` called `VoxelGrid.mark_occupied()`, a method
that does not exist on the class. As a result the voxel filter shown as
Stage-1 in Fig. 2(b) of the paper was never executed and every candidate sphere
went straight to the octree search. The voxel grid is now created in
`GeneratePBCell` and updated whenever a sphere is accepted.

This affects run time only. Sphere positions produced by v1.0 were valid;
they were simply obtained without the intended acceleration.

**Volume averaging was mis-weighted for elements with several integration
points.**
`generate_ss_curve` built its volume lookup as `{v.elementLabel: v.data}` from
the `IVOL` field. `IVOL` is an *integration-point* volume, so for a C3D10,
C3D10M or C3D8 element the dictionary retained only the last integration
point's volume and applied it to the stress at every integration point of that
element. Eq. (11) was therefore evaluated with incorrect weights. The lookup is
now keyed on `(elementLabel, integrationPoint)`.

The error is small but systematic. On a synthetic 4-integration-point mesh the
revised implementation reproduces Eq. (11) to machine precision, while the v1.0
implementation deviates by roughly 0.08 %. Homogenized curves generated with
v1.0 on a tetrahedral mesh with a single integration point (C3D4) are
unaffected.

---

### Performance

#### Sphere generation

*Voxel filter (Stage-1).* Rewritten as an exact rejection test rather than a
conservative one. A sphere marks every voxel fully contained in a ball of
radius `r2 + min_dist + r_pad - half_diagonal`, where
`r_pad = max(0, r_avg - 3*r_std)` is a lower bound on the radius a candidate
can have. A candidate with `r >= r_pad` whose centre falls in a marked voxel is
then provably overlapping and can be discarded without entering the octree;
candidates smaller than `r_pad` skip Stage-1 and are handled by the exact test.
The filter therefore never rejects a valid candidate. The grid is capped at
128³ voxels (about 2 MB) so that a very small `min_distance` cannot exhaust
memory. Marking was also moved out of `is_overlapping`, so a query no longer
mutates state.

*Octree (Stage-2).* Only the original spheres are inserted. Periodicity is
handled by `OctreeNode.query_periodic`, which prunes a node using the
minimum-image distance from the query centre to that node's bounding box.
Previously the periodic images were inserted as separate entries to make the
plain bounding-box query reachable across the cell faces; storing them is no
longer necessary, which shrinks the tree by roughly 1.5–2× at high volume
fraction, and the periodic distance test is a tighter prune than the
axis-aligned overlap it replaces.

For the same reason the mirrored spheres are no longer re-tested for overlap:
under the minimum-image convention already used in `is_overlapping`, testing an
image at `c + nL` gives exactly the same result as testing the original at `c`.

*Bookkeeping.* `Sphere_array` is accumulated in a Python list and converted to
an array once, instead of `np.vstack` inside the loop, which copied the whole
array on every acceptance (O(N²)). Removal of duplicate periodic images uses a
hash set instead of an O(k²) `np.allclose` scan.

*Vectorization.* Candidate spheres are drawn in batches. The radius positivity
check, the Stage-1 voxel test and the boundary overlap-portion rule are
evaluated with NumPy over a whole batch, and the distance test inside
`is_overlapping` is vectorized over the candidates returned by the octree. The
batch size adapts between 64 and 4096 as the acceptance rate collapses near the
target volume fraction.

`iterations` still counts the number of candidate spheres drawn, so
`max_iterations` has exactly the same meaning as in v1.0.

#### ODB post-processing

`generate_ss_curve` reads `S` and `IVOL` through `FieldOutput.bulkDataBlocks`,
which returns NumPy arrays, instead of `FieldOutput.values`, which constructs
one Python object per integration point. The permutation that aligns the two
fields is computed on the first frame and reused for the remaining frames, and
the accumulation of Eq. (11) is a single `numpy.dot`.

---

### Benchmarks

Sphere generation, `L = 0.2`, `min_dist = 0.001`, `max_iter = 100000`,
identical random seed. The number of real overlaps in the resulting array was
verified to be zero in every case (minimum-image distance, after de-duplicating
periodic images).

| `r_avg` | `r_std` | target VF | v1.0 | v1.1 | speed-up |
|---|---|---|---|---|---|
| 0.040 | 0.005  | 0.2 | 0.00 s | 0.06 s | — |
| 0.040 | 0.005  | 0.3 | 0.17 s | 0.09 s | 1.9× |
| 0.040 | 0.005  | 0.4 | 4.54 s | 0.16 s | 28× |
| 0.015 | 0.002  | 0.2 | 0.16 s | 0.11 s | 1.5× |
| 0.015 | 0.002  | 0.3 | 4.95 s | 0.44 s | 11× |
| 0.015 | 0.002  | 0.4 | 21.8 s | 1.15 s | 19× |
| 0.010 | 0.0015 | 0.2 | 0.79 s | 0.27 s | 2.9× |
| 0.010 | 0.0015 | 0.3 | 28.9 s | 2.72 s | 11× |
| 0.010 | 0.0015 | 0.4 | 41.4 s | 4.72 s | 8.8× |

The first row of the table corresponds to the input parameters of Table 5 in
the paper. Achieved volume fractions are statistically identical to v1.0, and
the approximately 40 % ceiling of the placement algorithm reported in the paper
is unchanged — it is a property of random sequential placement, not of the
implementation.

ODB post-processing, 500 000 integration points (about 125 000 C3D10 elements):
**3.21 s → 0.017 s per frame**.

---

### Known limitations

The octree search radius `r + min_dist + L/8` uses `L/8` in place of the
largest radius present in the cell. When `r_avg` is large relative to `L` this
under-estimates that radius — for the parameters of Table 5 in the paper,
`L/8 = 0.025` against `r_avg + 3*r_std = 0.055` — so in principle a pair of
large spheres can fall outside the near-by search. No overlap was observed in
any of the benchmark cases above, and the expression is kept as published. It
is a candidate for revision in a future version, together with a placement
algorithm able to exceed the current volume-fraction ceiling (for example a
force-biased or Lubachevsky–Stillinger scheme).

---

### Verification

To confirm the installation, regenerate the RVEs of Table 5 (volume fractions
of 10, 20, 30 and 40 %) and reproduce Fig. 6 and Fig. 7 of the paper. The
homogenized stress–strain curves and the predicted elastic moduli should
overlay the published results.
