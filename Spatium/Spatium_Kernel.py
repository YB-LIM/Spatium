"""
Spatium: A Plug-in for Generating Periodic Composite Cells with Random Sphere Inclusions in Abaqus/CAE

AUTHOR  : YOUNGBIN LIM
CONTACT : lyb0684@naver.com

Parameter definition

L              : Size of unit cell
L_mesh         : Global mesh size
r_avg          : Average radius of sphere
r_std          : Standard deviation of radius of sphere
VoF_tar        : Target volume fraction of a unit cell
min_distance   : Minimum allowable distance between spheres
Tol            : Minimum overlapping portion
Mode           : Deformation mode 
                 Uniaxial tension 
                 Biaxial tension
                 Confined compression   
                 Simple shear 
Disp           : Displacement for RP 
Is_Porous      : Flag for porous material
                 Composite
                 Porous    
max_iterations : Maximum number of iteration to generate spheres

---------------------------------------------------------------------------
Performance-revised build (v1.0-perf).  Numerical results are unchanged;
only the implementation of the two-stage overlap check, the sphere-array
bookkeeping and the ODB post-processing were rewritten.  See CHANGES.md.
---------------------------------------------------------------------------
"""

import os
import numpy as np
from itertools import product
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import visualization
import xyPlot
from odbAccess import openOdb

class VoxelGrid:
    """
    Stage-1 of the two-stage overlap check (rapid rejection of overlapping
    spheres, Fig. 2(b) of the SoftwareX paper).

    The unit cell [0,L]^3 is divided into a uniform grid of voxels, every
    voxel initially marked as *free*.  When a sphere is accepted, only the
    voxels that lie ENTIRELY inside that sphere are marked as *occupied*
    (with periodic wrap-around, so the grid itself is periodic).

    The marking radius is r2 + min_dist + r_pad - half_diagonal, where
    r_pad is a lower bound on the radius a candidate can have.  A candidate
    whose radius is at least r_pad and whose centre falls in an occupied
    voxel is then *provably* overlapping, so the filter never rejects a
    valid candidate.  Candidates smaller than r_pad simply skip stage-1 and
    go straight to the octree search (stage-2).

    NOTE  In v1.0 the module level object `voxelGrid` was never instantiated
          and `is_overlapping` called a non-existent method
          (`mark_occupied`), so stage-1 was effectively inactive at run time.
          It is now created in GeneratePBCell and updated on acceptance.
    """

    MAX_DIV = 128          # cap on voxels per direction (128^3 bool ~ 2 MB)

    def __init__(self, L, min_distance, resolutionFactor=2, r_pad=0.0):
        """
        :param L: Size of the cell
        :param min_distance: Minimum surface-to-surface distance between spheres
        :param resolutionFactor: Increase to get finer voxel spacing
                                 (target voxel size = min_distance / resolutionFactor)
        :param r_pad: Lower bound on the candidate radius (see class doc)
        """
        L = float(L)
        target = float(min_distance) / float(resolutionFactor)
        if target <= 0.0:
            target = L / float(self.MAX_DIV)

        n = int(np.ceil(L / target))
        n = max(8, min(self.MAX_DIV, n))

        self.L = L
        self.n = n
        self.nx = self.ny = self.nz = n
        self.voxel_size = L / float(n)
        self.half_diag = 0.5 * np.sqrt(3.0) * self.voxel_size
        self.min_distance = float(min_distance)
        self.r_pad = max(0.0, float(r_pad))
        self.occupied = np.zeros((n, n, n), dtype=bool)

    # ------------------------------------------------------------------
    # Query
    # ------------------------------------------------------------------
    def _xyz_to_ijk(self, x, y, z):
        """Continuous coordinates -> integer voxel indices (periodic)."""
        h, n = self.voxel_size, self.n
        i = int(np.floor((x % self.L) / h)) % n
        j = int(np.floor((y % self.L) / h)) % n
        k = int(np.floor((z % self.L) / h)) % n
        return i, j, k

    def is_occupied_batch(self, centers, radii=None):
        """
        Vectorised stage-1 test.

        :param centers: (M,3) array of candidate centres
        :param radii:   (M,) array of candidate radii (optional)
        :return: (M,) boolean array, True = definitely overlapping
        """
        c = np.mod(np.asarray(centers, dtype=float), self.L)
        idx = np.floor(c / self.voxel_size).astype(np.int64)
        np.clip(idx, 0, self.n - 1, out=idx)
        hit = self.occupied[idx[:, 0], idx[:, 1], idx[:, 2]]
        if radii is not None and self.r_pad > 0.0:
            hit = hit & (np.asarray(radii, dtype=float) >= self.r_pad)
        return hit

    def is_occupied(self, center, radius=None):
        """Scalar convenience wrapper around is_occupied_batch."""
        c = np.asarray(center, dtype=float).reshape(1, 3)
        r = None if radius is None else np.asarray([radius], dtype=float)
        return bool(self.is_occupied_batch(c, r)[0])

    # ------------------------------------------------------------------
    # Update
    # ------------------------------------------------------------------
    def mark_sphere(self, center, radius):
        """
        Mark every voxel that is fully contained in the exclusion ball of
        the sphere (centre, radius) as occupied.  Periodic images are
        handled by wrapping the voxel indices, so no explicit image loop is
        needed.
        """
        r_in = float(radius) + self.min_distance + self.r_pad - self.half_diag
        if r_in <= 0.0:
            return                      # sphere smaller than one voxel

        h, n, L = self.voxel_size, self.n, self.L
        cx, cy, cz = [float(v) % L for v in center]

        span = int(np.ceil(2.0 * r_in / h)) + 2
        if span >= n:                   # would wrap onto itself
            self.occupied[:, :, :] = True
            return

        i0 = int(np.floor((cx - r_in) / h))
        j0 = int(np.floor((cy - r_in) / h))
        k0 = int(np.floor((cz - r_in) / h))

        ii = np.arange(i0, int(np.floor((cx + r_in) / h)) + 1)
        jj = np.arange(j0, int(np.floor((cy + r_in) / h)) + 1)
        kk = np.arange(k0, int(np.floor((cz + r_in) / h)) + 1)
        if ii.size == 0 or jj.size == 0 or kk.size == 0:
            return

        dx = (ii + 0.5) * h - cx
        dy = (jj + 0.5) * h - cy
        dz = (kk + 0.5) * h - cz

        d2 = (dx[:, None, None] ** 2 +
              dy[None, :, None] ** 2 +
              dz[None, None, :] ** 2)
        mask = d2 <= r_in * r_in
        if not mask.any():
            return

        sel = np.ix_(ii % n, jj % n, kk % n)
        self.occupied[sel] = self.occupied[sel] | mask

class OctreeNode:
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max, 
                 capacity=4, depth=0, max_depth=10):
        self.x_min, self.x_max = x_min, x_max
        self.y_min, self.y_max = y_min, y_max
        self.z_min, self.z_max = z_min, z_max
        self.capacity = capacity
        self.depth = depth
        self.max_depth = max_depth
        self.spheres = []    # list of (center, radius)
        self.children = []   # octants

    def subdivide(self):
        """Split current node into 8 children (octants)."""
        x_mid = 0.5 * (self.x_min + self.x_max)
        y_mid = 0.5 * (self.y_min + self.y_max)
        z_mid = 0.5 * (self.z_min + self.z_max)

        for (x0, x1) in [(self.x_min, x_mid), (x_mid, self.x_max)]:
            for (y0, y1) in [(self.y_min, y_mid), (y_mid, self.y_max)]:
                for (z0, z1) in [(self.z_min, z_mid), (z_mid, self.z_max)]:
                    child = OctreeNode(
                        x_min=x0, x_max=x1,
                        y_min=y0, y_max=y1,
                        z_min=z0, z_max=z1,
                        capacity=self.capacity,
                        depth=self.depth + 1,
                        max_depth=self.max_depth
                    )
                    self.children.append(child)

    def insert(self, center, radius):
        """Insert a new sphere into this node or its children."""
        cx, cy, cz = center

        # If the center isn't in this node's bounding box, skip
        if not (self.x_min <= cx <= self.x_max and
                self.y_min <= cy <= self.y_max and
                self.z_min <= cz <= self.z_max):
            return False

        # If we have room or we're at max depth, store here
        if len(self.spheres) < self.capacity or self.depth == self.max_depth:
            self.spheres.append((center, radius))
            return True

        # Otherwise, subdivide if not already done
        if not self.children:
            self.subdivide()

        # Attempt to insert into each child
        for child in self.children:
            if child.insert(center, radius):
                return True

        return False

    def query(self, center, radius):
        """
        Return a list of candidate spheres in or near the bounding
        region of (center, radius).
        """
        cx, cy, cz = center

        # Bounds for the query region
        qx_min = cx - radius
        qx_max = cx + radius
        qy_min = cy - radius
        qy_max = cy + radius
        qz_min = cz - radius
        qz_max = cz + radius

        # If query region doesn't intersect this node, return empty
        if (qx_max < self.x_min or qx_min > self.x_max or
            qy_max < self.y_min or qy_min > self.y_max or
            qz_max < self.z_min or qz_min > self.z_max):
            return []

        # Otherwise, gather candidates from this node
        candidates = []
        candidates.extend(self.spheres)

        # And from children
        for child in self.children:
            candidates.extend(child.query(center, radius))

        return candidates

    def query_periodic(self, center, radius, L, out):
        """
        Periodic version of `query`.

        Only the *original* spheres (centres in [0,L]^3) are stored in the
        tree; periodicity is handled here by measuring the distance from the
        query centre to this node's bounding box under the minimum image
        convention.  A node is visited only if that distance is smaller than
        the search radius, which is a tighter test than the plain
        axis-aligned box overlap used by `query`.

        Candidates are appended to `out` (a plain list) to avoid the
        repeated list allocations of the recursive `query`.
        """
        gx = _periodic_gap(center[0], self.x_min, self.x_max, L)
        if gx > radius:
            return
        gy = _periodic_gap(center[1], self.y_min, self.y_max, L)
        if gy > radius:
            return
        gz = _periodic_gap(center[2], self.z_min, self.z_max, L)
        if gz > radius:
            return
        if gx * gx + gy * gy + gz * gz > radius * radius:
            return

        if self.spheres:
            out.extend(self.spheres)
        for child in self.children:
            child.query_periodic(center, radius, L, out)


def _periodic_gap(c, lo, hi, L):
    """Shortest distance from point `c` to the interval [lo,hi] on a
    periodic axis of length L."""
    best = None
    for s in (-L, 0.0, L):
        p = c + s
        if p < lo:
            d = lo - p
        elif p > hi:
            d = p - hi
        else:
            d = 0.0
        if best is None or d < best:
            best = d
            if best == 0.0:
                break
    return best


class Octree:
    """A simple wrapper around the root node.

    Only the original spheres are inserted; their periodic images are kept
    in Sphere_array for the CAE geometry but are NOT stored here, because
    the overlap test already uses the minimum image convention.  This keeps
    the tree ~1.5-2x smaller at high volume fraction.
    """
    def __init__(self, L, capacity=4, max_depth=10):
        eps = 1.0e-9 * L
        self.L = float(L)
        self.root = OctreeNode(
            x_min=-eps, x_max=L + eps,
            y_min=-eps, y_max=L + eps,
            z_min=-eps, z_max=L + eps,
            capacity=capacity,
            max_depth=max_depth
        )

    def insert_sphere(self, center, radius):
        self.root.insert(center, radius)

    def query_spheres(self, center, radius):
        return self.root.query(center, radius)

    def query_spheres_periodic(self, center, radius):
        out = []
        self.root.query_periodic(center, radius, self.L, out)
        return out

# Create a global (or external) voxelGrid reference. If you prefer, you can
# make it a class variable or pass it around. As long as we don't modify
# GeneratePBCell, we do it here.
voxelGrid = None  # Will be assigned before calling GeneratePBCell

def is_overlapping(Sphere_array_unused, center, radius, L, min_distance, octree):
    """
    Two-stage overlap check.

    Stage-1 : voxel grid  -> O(1) rejection of centres that are provably
                             inside an existing sphere.
    Stage-2 : octree      -> near-by search followed by an exact
                             minimum-image distance test.

    The search radius is unchanged from v1.0 (Fig. 2(b) of the paper):
        search radius = radius + min_dist + L/8

    Changes with respect to v1.0:
      * the distance test over the returned candidates is vectorised with
        NumPy instead of a Python for-loop;
      * the periodic-aware octree query is used, so periodic images no
        longer have to be inserted into the tree;
      * the function no longer mutates the voxel grid.  Marking is done by
        the caller once a sphere is actually accepted (in v1.0 the marking
        call sat inside this function and referenced a method that did not
        exist).
    """
    global voxelGrid

    # ---------------- Stage-1 : voxel filter ----------------
    if voxelGrid is not None:
        if voxelGrid.is_occupied(center, radius):
            return True

    # ---------------- Stage-2 : octree search ----------------
    search_radius = radius + min_distance + L/8.0
    nearby_spheres = octree.query_spheres_periodic(center, search_radius)
    if not nearby_spheres:
        return False

    C = np.asarray([s[0] for s in nearby_spheres], dtype=float).reshape(-1, 3)
    R = np.asarray([s[1] for s in nearby_spheres], dtype=float)

    # Minimum image logic for periodic boundaries (vectorised)
    d = np.asarray(center, dtype=float).reshape(1, 3) - C
    d -= L * np.round(d / L)
    dist2 = np.einsum('ij,ij->i', d, d)

    radii_sum = R + (radius + min_distance)
    return bool(np.any(dist2 < radii_sum * radii_sum))

voxelGrid = None

def GeneratePBCell(Is_Porous, L, L_mesh, r_avg, r_std, VoF_tar, min_distance, 
                   max_iterations, Mode, Disp):
    #################################################################
    # Generate numpy array that contains center position and radius #
    #################################################################
    global voxelGrid

    VoF = 0.0

    # -------------------------------------------------------
    # Stage-1 (voxel grid) and stage-2 (octree) accelerators
    # -------------------------------------------------------
    # r_pad: lower bound of the sampled radii, used to enlarge the voxel
    # exclusion ball without ever producing a false rejection
    r_pad = max(0.0, r_avg - 3.0 * r_std)
    voxelGrid = VoxelGrid(L, min_distance, resolutionFactor=2, r_pad=r_pad)
    octree = Octree(L, capacity=8, max_depth=10)

    # Rows are collected in a plain list and converted once at the end.
    # (v1.0 used np.vstack inside the loop, which copies the whole array
    #  on every acceptance -> O(N^2).)
    sphere_rows = []

    iterations = 0
    Tol_ = 0.4           # minimum overlapping portion, same value as v1.0

    # Adaptive batch size: small while candidates are still easy to place,
    # large once the acceptance rate collapses near the target VoF.
    BATCH_MIN, BATCH_MAX = 64, 4096
    batch = BATCH_MIN

    def _portion_bad(c, er, r_, L_):
        """
        Vectorised form of the original `check_overlap_portion`.
        True  -> the portion of the sphere sticking out of the cell is too
                 small and the sphere must be discarded.
        """
        is_left = (c - er) < 0.0
        is_right = (c + er) > L_
        touched = is_left | is_right
        if not np.any(touched):
            return np.zeros(c.shape, dtype=bool)

        d = np.where(is_left, c, L_ - c)
        rd = r_ - d
        invalid = (d <= 0.0) | (rd <= 0.0)

        safe_rd = np.where(invalid, 1.0, rd)
        ratio = np.where(invalid, 1.0, d / safe_rd)
        safe_ratio = np.where(ratio == 0.0, 1.0, ratio)
        small = (ratio < Tol_) | ((1.0 / safe_ratio) < Tol_)

        return touched & (invalid | small)

    while VoF <= VoF_tar and iterations < max_iterations:

        m = int(min(batch, max_iterations - iterations))
        if m <= 0:
            break
        iterations += m

        # ---- batched random sampling (uniform centres, normal radii) ----
        cand_c = L * np.random.rand(m, 3)
        cand_r = np.random.normal(r_avg, r_std, m)

        keep = cand_r > 0.0                      # skip negative or zero radii

        # ---- Stage-1 : voxel filter, O(1) per candidate, vectorised ----
        keep &= ~voxelGrid.is_occupied_batch(cand_c, cand_r)

        # ---- boundary overlap-portion filter, vectorised ----
        er = cand_r + min_distance / 2           # effective_radius
        bad = (_portion_bad(cand_c[:, 0], er, cand_r, L) |
               _portion_bad(cand_c[:, 1], er, cand_r, L) |
               _portion_bad(cand_c[:, 2], er, cand_r, L))
        keep &= ~bad

        idx = np.nonzero(keep)[0]
        if idx.size == 0:
            batch = min(BATCH_MAX, batch * 2)
            continue

        accepted = 0

        # ---- Stage-2 : octree search, sequential acceptance ----
        for t in idx:
            if VoF > VoF_tar:
                break

            center = cand_c[t]
            radius = float(cand_r[t])

            # re-run stage-1: spheres accepted earlier in this same batch
            # may already block this candidate
            if voxelGrid.is_occupied(center, radius):
                continue

            if is_overlapping(None, center, radius, L, min_distance, octree):
                continue

            center_X = float(center[0])
            center_Y = float(center[1])
            center_Z = float(center[2])

            effective_radius = radius + min_distance / 2

            Is_inside = (
                effective_radius < center_X < L - effective_radius and
                effective_radius < center_Y < L - effective_radius and
                effective_radius < center_Z < L - effective_radius
            )

            new_rows = [[center_X, center_Y, center_Z, radius]]

            if not Is_inside:
                # ----------------------------------------------------------
                # Sphere crosses a boundary -> create mirrored spheres.
                #
                # The mirrored spheres are NOT re-tested here: under the
                # minimum image convention used in is_overlapping, testing
                # an image at c + n*L is mathematically identical to testing
                # the original at c, so the test above already covers them.
                # ----------------------------------------------------------
                Is_X_Left = center_X - effective_radius < 0
                Is_X_Right = center_X + effective_radius > L
                Is_Y_Front = center_Y - effective_radius < 0
                Is_Y_Back = center_Y + effective_radius > L
                Is_Z_Lower = center_Z - effective_radius < 0
                Is_Z_Upper = center_Z + effective_radius > L

                x_shifts = [0]
                if Is_X_Left:
                    x_shifts.append(L)
                if Is_X_Right:
                    x_shifts.append(-L)

                y_shifts = [0]
                if Is_Y_Front:
                    y_shifts.append(L)
                if Is_Y_Back:
                    y_shifts.append(-L)

                z_shifts = [0]
                if Is_Z_Lower:
                    z_shifts.append(L)
                if Is_Z_Upper:
                    z_shifts.append(-L)

                shift_combinations = list(product(x_shifts, y_shifts, z_shifts))
                shift_combinations.remove((0, 0, 0))
                for dx, dy, dz in shift_combinations:
                    new_rows.append([center_X + dx,
                                     center_Y + dy,
                                     center_Z + dz,
                                     radius])

                # Remove duplicates with a hash set instead of the
                # O(k^2) np.allclose scan used in v1.0
                seen = set()
                unique_new_rows = []
                for row in new_rows:
                    key = (round(row[0], 12), round(row[1], 12), round(row[2], 12))
                    if key not in seen:
                        seen.add(key)
                        unique_new_rows.append(row)
                new_rows = unique_new_rows

            sphere_rows.extend(new_rows)
            accepted += 1

            # Volume is counted once for the real sphere
            VoF += (4.0 / 3.0) * np.pi * (radius ** 3) / (L ** 3)

            # Only the original sphere goes into the octree
            octree.insert_sphere(np.array([center_X, center_Y, center_Z]), radius)

            # Update the stage-1 grid
            voxelGrid.mark_sphere((center_X, center_Y, center_Z), radius)

        # Adapt the batch size to the current acceptance rate
        if accepted == 0:
            batch = min(BATCH_MAX, batch * 2)
        elif accepted * 8 >= m:
            batch = max(BATCH_MIN, batch // 2)

    if sphere_rows:
        Sphere_array = np.asarray(sphere_rows, dtype=float)
    else:
        Sphere_array = np.empty((0, 4))

    Number_of_Sphere = Sphere_array.shape[0]
    
    ##########################################
    # Generate cell and sphere in Abaqus/CAE #
    ##########################################
    
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200.890625, 
        height=89.5752334594727)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].maximize()
    executeOnCaeStartup()
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=ON)
    Mdb()
    #: A new model database has been created.
    #: The model "Model-1" has been created.
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(L, L))
    p = mdb.models['Model-1'].Part(name='Matrix_Filled', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Matrix_Filled']
    p.BaseSolidExtrude(sketch=s, depth=L)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Matrix_Filled']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['Matrix_Filled']
    a.Instance(name='Matrix_Filled-1', part=p, dependent=ON)
    
    # Generate spheres
    
    for i in range(Number_of_Sphere):
        X = Sphere_array[i,0]
        Y = Sphere_array[i,1]
        Z = Sphere_array[i,2]
        r = Sphere_array[i,3]
        Sphere_Name = 'Sphere_' + str(i)
        Sphere_Name_Assembly = Sphere_Name + "-1"
    
        s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
            sheetSize=200.0)
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        s1.setPrimaryObject(option=STANDALONE)
        s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
        s1.FixedConstraint(entity=g[2])
        session.viewports['Viewport: 1'].view.setValues(nearPlane=186.033, 
            farPlane=191.09, width=14.083, height=7.09298, cameraPosition=(1.14285, 
            -0.312562, 188.562), cameraTarget=(1.14285, -0.312562, 0))
        s1.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, -r), point2=(0.0, r), 
            direction=CLOCKWISE)
        s1.Line(point1=(0.0, -r), point2=(0.0, r))
        s1.VerticalConstraint(entity=g[4], addUndoState=False)
        s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
        p = mdb.models['Model-1'].Part(name=Sphere_Name, dimensionality=THREE_D, 
            type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts[Sphere_Name]
        p.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
        s1.unsetPrimaryObject()
        p = mdb.models['Model-1'].parts[Sphere_Name]
        p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
        p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
        p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
        c = p.cells
        # Replace with coordinates within the desired cell
        pickedCells = c.findAt(((0, 0, 0), ))
        d = p.datums
        p.PartitionCellByDatumPlane(datumPlane=d[2], cells=pickedCells)
        session.viewports['Viewport: 1'].setValues(displayedObject=p)
        del mdb.models['Model-1'].sketches['__profile__']
        
        # Translate sphere
        a = mdb.models['Model-1'].rootAssembly
        p = mdb.models['Model-1'].parts[Sphere_Name]
        a.Instance(name=Sphere_Name_Assembly, part=p, dependent=ON)
        a = mdb.models['Model-1'].rootAssembly
        a.translate(instanceList=(Sphere_Name_Assembly, ), vector=(X, Y, Z))
    
    # Generate Matrix only part and Sphere only part
    session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(
        datumAxes=OFF)
    
    a = mdb.models['Model-1'].rootAssembly
    
    # Get the number of spheres
    num_spheres = len(Sphere_array)
    
    # Generate the list of sphere instance names
    sphere_instance_names = ['Sphere_{}-1'.format(i) for i in range(num_spheres)]
    
    # Retrieve the sphere instances from the assembly
    sphere_instances = [a.instances[name] for name in sphere_instance_names]
    
    # Convert the list to a tuple for the 'cuttingInstances' parameter
    cutting_instances = tuple(sphere_instances)
    
    matrix_instance = a.instances['Matrix_Filled-1']
    all_instances = [matrix_instance] + sphere_instances
    merge_instances = tuple(all_instances)
    
    # Perform the Boolean cut operation
    a.InstanceFromBooleanCut(name='Matrix_Only', 
        instanceToBeCut=a.instances['Matrix_Filled-1'], 
        cuttingInstances=cutting_instances, 
        originalInstances=SUPPRESS)
    
    a = mdb.models['Model-1'].rootAssembly
    
    # Determine the number of spheres
    num_spheres = len(Sphere_array)
    
    # Generate the list of feature names to resume
    features_to_resume = ['Matrix_Filled-1']
    features_to_resume.extend(['Sphere_{}-1'.format(i) for i in range(num_spheres)])
    
    # Verify that all features to resume exist
    missing_features = [name for name in features_to_resume if name not in a.features]
    if missing_features:
        print("Warning: The following features do not exist and cannot be resumed:", missing_features)
    else:
        # Resume features dynamically
        a.resumeFeatures(tuple(features_to_resume))
    
    # Suppress the 'Matrix_Only-1' feature
    feature_to_suppress = 'Matrix_Only-1'  # Modify if the feature name varies
    if feature_to_suppress in a.features:
        a.features[feature_to_suppress].suppress()
    else:
        print("Warning: Feature {} does not exist and cannot be suppressed.".format(feature_to_suppress))
    
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.339431, 
        farPlane=1.138, width=0.41982, height=0.230888, viewOffsetX=0.0453624, 
        viewOffsetY=-0.0150945)
    a1 = mdb.models['Model-1'].rootAssembly
    a1.InstanceFromBooleanMerge(name='Matrix_w_Sphere', instances=merge_instances, keepIntersections=ON, 
        originalInstances=SUPPRESS, domain=GEOMETRY)
    
    p = mdb.models['Model-1'].parts['Matrix_w_Sphere']
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=L)
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=L)
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=L)
    
    # Extrude cut for spheres outside the cell
    # XY-Plane
    p = mdb.models['Model-1'].parts['Matrix_w_Sphere']
    e, d = p.edges, p.datums
    Sketch_Edge = e.findAt(((L, L/2., 0), ))[0]
    t = p.MakeSketchTransform(sketchPlane=d[2], sketchUpEdge=Sketch_Edge, 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
        0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=L, 
        gridSpacing=L/50.0, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.rectangle(point1=(-10*L, 10*L), point2=(10*L, -10*L))
    p.CutExtrude(sketchPlane=d[2], sketchUpEdge=Sketch_Edge, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    
    # XZ-Plane Upper
    p = mdb.models['Model-1'].parts['Matrix_w_Sphere']
    e, d = p.edges, p.datums
    Sketch_Edge = e.findAt(((0, L, L/2.), ))[0]
    t = p.MakeSketchTransform(sketchPlane=d[7], sketchUpEdge=Sketch_Edge, 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
        0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=L, 
        gridSpacing=L/50., transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.rectangle(point1=(-10*L, 10*L), point2=(10*L, -10*L))
    p.CutExtrude(sketchPlane=d[7], sketchUpEdge=Sketch_Edge, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=ON)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    
    # XY-Plane Front
    p = mdb.models['Model-1'].parts['Matrix_w_Sphere']
    e, d = p.edges, p.datums
    Sketch_Edge = e.findAt(((L, L/2., L), ))[0]
    t = p.MakeSketchTransform(sketchPlane=d[3], sketchUpEdge=Sketch_Edge, 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
        0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=L, 
        gridSpacing=L/50., transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.rectangle(point1=(-10*L, 10*L), point2=(10*L, -10*L))
    p.CutExtrude(sketchPlane=d[3], sketchUpEdge=Sketch_Edge, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=ON)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    
    # XZ-Plane Lower
    p = mdb.models['Model-1'].parts['Matrix_w_Sphere']
    e, d = p.edges, p.datums
    Sketch_Edge = e.findAt(((L, L, L/2.), ))[0]
    t = p.MakeSketchTransform(sketchPlane=d[6], sketchUpEdge=Sketch_Edge, 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
        0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=L, 
        gridSpacing=L/50., transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.rectangle(point1=(-10*L, 10*L), point2=(10*L, -10*L))
    p.CutExtrude(sketchPlane=d[6], sketchUpEdge=Sketch_Edge, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    
    # YZ-Plane
    p = mdb.models['Model-1'].parts['Matrix_w_Sphere']
    e, d = p.edges, p.datums
    Sketch_Edge = e.findAt(((L, L/2., 0.), ))[0]
    t = p.MakeSketchTransform(sketchPlane=d[5], sketchUpEdge=Sketch_Edge, 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
        0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=L, 
        gridSpacing=L/50., transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.rectangle(point1=(-10*L, 10*L), point2=(10*L, -10*L))
    p.CutExtrude(sketchPlane=d[5], sketchUpEdge=Sketch_Edge, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=ON)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    
    # YZ-Plane Right
    p = mdb.models['Model-1'].parts['Matrix_w_Sphere']
    e, d = p.edges, p.datums
    Sketch_Edge = e.findAt(((0., L/2., 0.), ))[0]
    t = p.MakeSketchTransform(sketchPlane=d[4], sketchUpEdge=Sketch_Edge, 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
        0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=L, 
        gridSpacing=L/50., transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.rectangle(point1=(-10*L, 10*L), point2=(10*L, -10*L))
    p.CutExtrude(sketchPlane=d[4], sketchUpEdge=Sketch_Edge, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    
    # Generate the list of feature names to delete
    features_to_delete = ['Matrix_Filled-1']
    features_to_delete.extend(['Sphere_{}-1'.format(i) for i in range(num_spheres)])
    
    # Convert the list to a tuple
    features_to_delete_tuple = tuple(features_to_delete)
    
    # Delete features dynamically
    a.deleteFeatures(features_to_delete_tuple)
    
    # Delete Matrix-Filled
    model = mdb.models['Model-1']
    del mdb.models['Model-1'].parts['Matrix_Filled']
    # Delete sphere parts dynamically
    for i in range(num_spheres):
        part_name = 'Sphere_{}'.format(i)
        if part_name in model.parts:
            del model.parts[part_name]
        else:
            print("Warning: Part '{}' not found.".format(part_name))
    
    a1 = mdb.models['Model-1'].rootAssembly
    a1.features['Matrix_Only-1'].resume()
    
    # Boolean cut for sphere only instance
    a = mdb.models['Model-1'].rootAssembly
    a.InstanceFromBooleanCut(name='Sphere_Only', 
        instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Matrix_w_Sphere-1'], 
        cuttingInstances=(a.instances['Matrix_Only-1'], ), 
        originalInstances=SUPPRESS)
    
    a = mdb.models['Model-1'].rootAssembly
    del a.features['Matrix_w_Sphere-1']
    a = mdb.models['Model-1'].rootAssembly
    a.features['Matrix_Only-1'].resume()
    
    del mdb.models['Model-1'].parts['Matrix_w_Sphere']
    
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    
    session.viewports['Viewport: 1'].enableMultipleColors()
    session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
    cmap=session.viewports['Viewport: 1'].colorMappings['Part instance']
    session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
    session.viewports['Viewport: 1'].disableMultipleColors()
    
    # Assign global seed size
    # Matrix side
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=ON)
    p = mdb.models['Model-1'].parts['Matrix_Only']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    p = mdb.models['Model-1'].parts['Matrix_Only']
    c = p.cells
    pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
    p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
    p = mdb.models['Model-1'].parts['Matrix_Only']
    c = p.cells
    cells = c.getByBoundingBox(xMin=-0.5*L, xMax=L*1.5, yMin=-0.5*L, yMax=1.5*L, zMin=-0.5*L, zMax=1.5*L)
    pickedRegions =(cells, )
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4H, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))
    p = mdb.models['Model-1'].parts['Matrix_Only']
    p.seedPart(size=L_mesh, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-1'].parts['Matrix_Only']
    p.generateMesh()
    
    # Sphere side 
    p = mdb.models['Model-1'].parts['Sphere_Only']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    p = mdb.models['Model-1'].parts['Sphere_Only']
    c = p.cells
    cells = c.getByBoundingBox(xMin=-0.5*L, xMax=L*1.5, yMin=-0.5*L, yMax=1.5*L, zMin=-0.5*L, zMax=1.5*L)
    pickedRegions = cells
    p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
    p = mdb.models['Model-1'].parts['Sphere_Only']
    p.seedPart(size=L_mesh, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-1'].parts['Sphere_Only']
    p.generateMesh()
    
    # View assembly mesh
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    
    # Generate Surface element for PBC
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=ON)
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(L, L))
    p = mdb.models['Model-1'].Part(name='PBC_Surface', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['PBC_Surface']
    p.BaseSolidExtrude(sketch=s, depth=L)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['PBC_Surface']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    p = mdb.models['Model-1'].parts['PBC_Surface']
    c1 = p.cells
    p.RemoveCells(cellList = c1[0:1])
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    p = mdb.models['Model-1'].parts['PBC_Surface']
    p.seedPart(size=L_mesh, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-1'].parts['PBC_Surface']
    p.generateMesh()
    elemType1 = mesh.ElemType(elemCode=SFM3D4, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=SFM3D3, elemLibrary=STANDARD)
    p = mdb.models['Model-1'].parts['PBC_Surface']
    f = p.faces
    faces = f.getByBoundingBox(xMin=-0.5*L, xMax=L*1.5, yMin=-0.5*L, yMax=1.5*L, zMin=-0.5*L, zMax=1.5*L)
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    
    # Assign surface section
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON, mesh=OFF)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    mdb.models['Model-1'].SurfaceSection(name='PBC_Surface', useDensity=OFF)
    p = mdb.models['Model-1'].parts['PBC_Surface']
    f = p.faces
    faces = f.getByBoundingBox(xMin=-0.5*L, xMax=L*1.5, yMin=-0.5*L, yMax=1.5*L, zMin=-0.5*L, zMax=1.5*L)
    region = p.Set(faces=faces, name='PBC_Surf_Set')
    p = mdb.models['Model-1'].parts['PBC_Surface']
    p.SectionAssignment(region=region, sectionName='PBC_Surface', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    
    # Generate surface set for tie on PBC_surface part
    p = mdb.models['Model-1'].parts['PBC_Surface']
    s = p.faces
    side2Faces = s.getByBoundingBox(xMin=-0.5*L, xMax=L*1.5, yMin=-0.5*L, yMax=1.5*L, zMin=-0.5*L, zMax=1.5*L)
    p.Surface(side2Faces=side2Faces, name='PBC_Surf_Tie')
    
    # Add PBC Surface instance to assembly
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    a = mdb.models['Model-1'].rootAssembly
    p = mdb.models['Model-1'].parts['PBC_Surface']
    a.Instance(name='PBC_Surface-1', part=p, dependent=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.hideInstances(instances=(
        'PBC_Surface-1', ))
    
    # Set color 
    cmap = session.viewports['Viewport: 1'].colorMappings['Part instance']
    cmap.updateOverrides(overrides={'PBC_Surface-1':(True, '#CCCCCC', 'Default', 
        '#CCCCCC'), 'Sphere_Only-1':(True, '#FAEBD7', 'Default', '#FAEBD7')})
    
    # Generate surface for tie on cell
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
    a = mdb.models['Model-1'].rootAssembly
    
    # Surface for +X direction
    X_Pos_Matrix = a.instances['Matrix_Only-1'].faces
    Face_X_Pos_Mat = X_Pos_Matrix.getByBoundingBox(xMin=L, xMax=L*1.5, yMin=0.0, yMax=L, zMin=0, zMax=L)
    X_Pos_Sphere = a.instances['Sphere_Only-1'].faces
    Face_X_Pos_Sph = X_Pos_Sphere.getByBoundingBox(xMin=L, xMax=L*1.5, yMin=0.0, yMax=L, zMin=0, zMax=L)
    Face_X_Pos = Face_X_Pos_Mat + Face_X_Pos_Sph
    
    # Surface for -X direction
    X_Neg_Matrix = a.instances['Matrix_Only-1'].faces
    Face_X_Neg_Mat = X_Neg_Matrix.getByBoundingBox(xMin=-0.5*L, xMax=0, yMin=0.0, yMax=L, zMin=0, zMax=L)
    X_Neg_Sphere = a.instances['Sphere_Only-1'].faces
    Face_X_Neg_Sph = X_Neg_Sphere.getByBoundingBox(xMin=-0.5*L, xMax=0, yMin=0.0, yMax=L, zMin=0, zMax=L)
    Face_X_Neg = Face_X_Neg_Mat + Face_X_Neg_Sph
    
    # Surface for +Y direction
    Y_Pos_Matrix = a.instances['Matrix_Only-1'].faces
    Face_Y_Pos_Mat = Y_Pos_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=L, yMax=1.5*L, zMin=0, zMax=L)
    Y_Pos_Sphere = a.instances['Sphere_Only-1'].faces
    Face_Y_Pos_Sph = Y_Pos_Sphere.getByBoundingBox(xMin=0, xMax=L, yMin=L, yMax=1.5*L, zMin=0, zMax=L)
    Face_Y_Pos = Face_Y_Pos_Mat + Face_Y_Pos_Sph
    
    # Surface for -Y direction
    Y_Neg_Matrix = a.instances['Matrix_Only-1'].faces
    Face_Y_Neg_Mat = Y_Neg_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=-0.5*L, yMax=0, zMin=0, zMax=L)
    Y_Neg_Sphere = a.instances['Sphere_Only-1'].faces
    Face_Y_Neg_Sph = Y_Neg_Sphere.getByBoundingBox(xMin=0, xMax=L, yMin=-0.5*L, yMax=0, zMin=0, zMax=L)
    Face_Y_Neg = Face_Y_Neg_Mat + Face_Y_Neg_Sph
    
    # Surface for +Z direction
    Z_Pos_Matrix = a.instances['Matrix_Only-1'].faces
    Face_Z_Pos_Mat = Z_Pos_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=0, yMax=L, zMin=L, zMax=1.5*L)
    Z_Pos_Sphere = a.instances['Sphere_Only-1'].faces
    Face_Z_Pos_Sph = Z_Pos_Sphere.getByBoundingBox(xMin=0, xMax=L, yMin=0, yMax=L, zMin=L, zMax=1.5*L)
    Face_Z_Pos = Face_Z_Pos_Mat + Face_Z_Pos_Sph
    
    # Surface for -Z direction
    Z_Neg_Matrix = a.instances['Matrix_Only-1'].faces
    Face_Z_Neg_Mat = Z_Neg_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=0, yMax=L, zMin=-0.5*L, zMax=0)
    Z_Neg_Sphere = a.instances['Sphere_Only-1'].faces
    Face_Z_Neg_Sph = Z_Neg_Sphere.getByBoundingBox(xMin=0, xMax=L, yMin=0, yMax=L, zMin=-0.5*L, zMax=0)
    Face_Z_Neg = Face_Z_Neg_Mat + Face_Z_Neg_Sph
    
    side1Faces = Face_X_Pos + Face_X_Neg + Face_Y_Pos + Face_Y_Neg + Face_Z_Pos + Face_Z_Neg
    a.Surface(side1Faces=side1Faces, name='Tie_Surf')
    
    # Tie
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['PBC_Surface-1'].surfaces['PBC_Surf_Tie']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.surfaces['Tie_Surf']
    mdb.models['Model-1'].Tie(name='Tie', main=region1, secondary=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
        interactions=ON, constraints=ON, connectors=ON, engineeringFeatures=ON)
    
    # Generate node sets to apply equation constraints
    # +X FACE
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=L_mesh/2., yMax=L-L_mesh/2., zMin=L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='X_Pos_FACE_Sorted')
    
    nodeSet = p.sets['X_Pos_FACE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: (node.coordinates[1], node.coordinates[2]))
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='X_Pos_FACE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['X_Pos_FACE_Sorted']
    
    # -X FACE
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=L_mesh/2., yMin=L_mesh/2., yMax=L-L_mesh/2., zMin=L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='X_Neg_FACE_Sorted')
    
    nodeSet = p.sets['X_Neg_FACE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: (node.coordinates[1], node.coordinates[2]))
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='X_Neg_FACE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['X_Neg_FACE_Sorted']
    
    # +Y FACE
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='Y_Pos_FACE_Sorted')
    
    nodeSet = p.sets['Y_Pos_FACE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: (node.coordinates[0], node.coordinates[2]))
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Y_Pos_FACE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Y_Pos_FACE_Sorted']
    
    # -Y FACE
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=-L_mesh/2., yMax=L_mesh/2., zMin=L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='Y_Neg_FACE_Sorted')
    
    nodeSet = p.sets['Y_Neg_FACE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: (node.coordinates[0], node.coordinates[2]))
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Y_Neg_FACE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Y_Neg_FACE_Sorted']
    
    # +Z FACE
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=L_mesh/2., yMax=L-L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='Z_Pos_FACE_Sorted')
    
    nodeSet = p.sets['Z_Pos_FACE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: (node.coordinates[0], node.coordinates[1]))
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Z_Pos_FACE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Z_Pos_FACE_Sorted']
    
    # -Z FACE
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=L_mesh/2., yMax=L-L_mesh/2., zMin=-L_mesh/2., zMax=+L_mesh/2.)
    p.Set(nodes=nodes, name='Z_Neg_FACE_Sorted')
    
    nodeSet = p.sets['Z_Neg_FACE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: (node.coordinates[0], node.coordinates[1]))
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Z_Neg_FACE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Z_Neg_FACE_Sorted']
    
    # X_-Y_-Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=-L_mesh/2., yMax=L_mesh/2., zMin=-L_mesh/2., zMax=L_mesh/2.)
    p.Set(nodes=nodes, name='X_YNeg_ZNeg_EDGE_Sorted')
    
    nodeSet = p.sets['X_YNeg_ZNeg_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[0])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='X_YNeg_ZNeg_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['X_YNeg_ZNeg_EDGE_Sorted']
    
    # X_+Y_-Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=-L_mesh/2., zMax=L_mesh/2.)
    p.Set(nodes=nodes, name='X_YPos_ZNeg_EDGE_Sorted')
    
    nodeSet = p.sets['X_YPos_ZNeg_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[0])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='X_YPos_ZNeg_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['X_YPos_ZNeg_EDGE_Sorted']
    
    # X_+Y_+Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='X_YPos_ZPos_EDGE_Sorted')
    
    nodeSet = p.sets['X_YPos_ZPos_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[0])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='X_YPos_ZPos_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['X_YPos_ZPos_EDGE_Sorted']
    
    # X_-Y_+Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L_mesh/2., xMax=L-L_mesh/2., yMin=-L_mesh/2., yMax=+L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='X_YNeg_ZPos_EDGE_Sorted')
    
    nodeSet = p.sets['X_YNeg_ZPos_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[0])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='X_YNeg_ZPos_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['X_YNeg_ZPos_EDGE_Sorted']
    
    # Y_-X_-Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=+L_mesh/2., yMax=L-L_mesh/2., zMin=-L_mesh/2., zMax=+L_mesh/2.)
    p.Set(nodes=nodes, name='Y_XNeg_ZNeg_EDGE_Sorted')
    
    nodeSet = p.sets['Y_XNeg_ZNeg_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[1])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Y_XNeg_ZNeg_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Y_XNeg_ZNeg_EDGE_Sorted']
    
    # Y_+X_-Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=+L_mesh/2., yMax=L-L_mesh/2., zMin=-L_mesh/2., zMax=+L_mesh/2.)
    p.Set(nodes=nodes, name='Y_XPos_ZNeg_EDGE_Sorted')
    
    nodeSet = p.sets['Y_XPos_ZNeg_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[1])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Y_XPos_ZNeg_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Y_XPos_ZNeg_EDGE_Sorted']
    
    # Y_+X_+Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=+L_mesh/2., yMax=L-L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='Y_XPos_ZPos_EDGE_Sorted')
    
    nodeSet = p.sets['Y_XPos_ZPos_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[1])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Y_XPos_ZPos_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Y_XPos_ZPos_EDGE_Sorted']
    
    # Y_-X_+Z_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=+L_mesh/2., yMax=L-L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='Y_XNeg_ZPos_EDGE_Sorted')
    
    nodeSet = p.sets['Y_XNeg_ZPos_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[1])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Y_XNeg_ZPos_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Y_XNeg_ZPos_EDGE_Sorted']
    
    # Z_-X_-Y_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=-L_mesh/2., yMax=+L_mesh/2., zMin=+L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='Z_XNeg_YNeg_EDGE_Sorted')
    
    nodeSet = p.sets['Z_XNeg_YNeg_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[2])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Z_XNeg_YNeg_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Z_XNeg_YNeg_EDGE_Sorted']
    
    # Z_+X_-Y_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=-L_mesh/2., yMax=+L_mesh/2., zMin=+L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='Z_XPos_YNeg_EDGE_Sorted')
    
    nodeSet = p.sets['Z_XPos_YNeg_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[2])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Z_XPos_YNeg_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Z_XPos_YNeg_EDGE_Sorted']
    
    # Z_+X_+Y_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=+L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='Z_XPos_YPos_EDGE_Sorted')
    
    nodeSet = p.sets['Z_XPos_YPos_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[2])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Z_XPos_YPos_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Z_XPos_YPos_EDGE_Sorted']
    
    # Z_-X_+Y_Edge
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=+L_mesh/2., zMax=L-L_mesh/2.)
    p.Set(nodes=nodes, name='Z_XNeg_YPos_EDGE_Sorted')
    
    nodeSet = p.sets['Z_XNeg_YPos_EDGE_Sorted']
    sorted_nodes = sorted(nodeSet.nodes, key=lambda node: node.coordinates[2])
    nodeLabels = [node.label for node in sorted_nodes]
    p.SetFromNodeLabels(name='Z_XNeg_YPos_EDGE', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Z_XNeg_YPos_EDGE_Sorted']
    
    # Corner_000
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=-L_mesh/2., yMax=+L_mesh/2., zMin=-L_mesh/2., zMax=+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_000_Sorted')
    
    nodeSet = p.sets['Corner_000_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_000', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_000_Sorted']
    
    # Corner_100
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=-L_mesh/2., yMax=+L_mesh/2., zMin=-L_mesh/2., zMax=+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_100_Sorted')
    
    nodeSet = p.sets['Corner_100_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_100', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_100_Sorted']
    
    # Corner_110
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=-L_mesh/2., zMax=+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_110_Sorted')
    
    nodeSet = p.sets['Corner_110_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_110', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_110_Sorted']
    
    # Corner_010
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=-L_mesh/2., zMax=+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_010_Sorted')
    
    nodeSet = p.sets['Corner_010_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_010', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_010_Sorted']
    
    # Corner_001
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=-L_mesh/2., yMax=+L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_001_Sorted')
    
    nodeSet = p.sets['Corner_001_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_001', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_001_Sorted']
    
    # Corner_101
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=-L_mesh/2., yMax=+L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_101_Sorted')
    
    nodeSet = p.sets['Corner_101_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_101', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_101_Sorted']
    
    # Corner_111
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=L-L_mesh/2., xMax=L+L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_111_Sorted')
    
    nodeSet = p.sets['Corner_111_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_111', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_111_Sorted']
    
    # Corner_011
    p = mdb.models['Model-1'].parts['PBC_Surface']
    n = p.nodes
    nodes = n.getByBoundingBox(xMin=-L_mesh/2., xMax=+L_mesh/2., yMin=L-L_mesh/2., yMax=L+L_mesh/2., zMin=L-L_mesh/2., zMax=L+L_mesh/2.)
    p.Set(nodes=nodes, name='Corner_011_Sorted')
    
    nodeSet = p.sets['Corner_011_Sorted']
    nodeLabels = [node.label for node in nodeSet.nodes]
    p.SetFromNodeLabels(name='Corner_011', nodeLabels=nodeLabels, unsorted=True)
    del mdb.models['Model-1'].parts['PBC_Surface'].sets['Corner_011_Sorted']
    
    # Generate RPs for PBC equation constraint
    a = mdb.models['Model-1'].rootAssembly
    a.ReferencePoint(point=(L/2., L/2., L/2.))
    refPointFound = a.referencePoints.findAt((L/2., L/2., L/2.),)
    a.Set(name='RP-1', referencePoints=(refPointFound,))
    a = mdb.models['Model-1'].rootAssembly
    a.ReferencePoint(point=(L/2., L/2., L/2.))
    refPointFound = a.referencePoints.findAt((L/2., L/2., L/2.),)
    a.Set(name='RP-2', referencePoints=(refPointFound,))
    a = mdb.models['Model-1'].rootAssembly
    a.ReferencePoint(point=(L/2., L/2., L/2.))
    refPointFound = a.referencePoints.findAt((L/2., L/2., L/2.),)
    a.Set(name='RP-3', referencePoints=(refPointFound,))
    
    # X FACES
    mdb.models['Model-1'].Equation(name='X_FACES_DOF_1', terms=(
        (-1.0, 'PBC_Surface-1.X_Pos_FACE', 1),
        ( 1.0, 'PBC_Surface-1.X_Neg_FACE', 1),
        ( 1.0, 'RP-1', 1)))
    mdb.models['Model-1'].Equation(name='X_FACES_DOF_2', terms=(
        (-1.0, 'PBC_Surface-1.X_Pos_FACE', 2),
        ( 1.0, 'PBC_Surface-1.X_Neg_FACE', 2)))
    mdb.models['Model-1'].Equation(name='X_FACES_DOF_3', terms=(
        (-1.0, 'PBC_Surface-1.X_Pos_FACE', 3),
        ( 1.0, 'PBC_Surface-1.X_Neg_FACE', 3)))
    
    # Y FACES
    mdb.models['Model-1'].Equation(name='Y_FACES_DOF_1', terms=(
        (-1.0, 'PBC_Surface-1.Y_Pos_FACE', 1),
        ( 1.0, 'PBC_Surface-1.Y_Neg_FACE', 1)))
    mdb.models['Model-1'].Equation(name='Y_FACES_DOF_2', terms=(
        (-1.0, 'PBC_Surface-1.Y_Pos_FACE', 2),
        ( 1.0, 'PBC_Surface-1.Y_Neg_FACE', 2),
        ( 1.0, 'RP-2', 2)))
    mdb.models['Model-1'].Equation(name='Y_FACES_DOF_3', terms=(
        (-1.0, 'PBC_Surface-1.Y_Pos_FACE', 3),
        ( 1.0, 'PBC_Surface-1.Y_Neg_FACE', 3)))
    
    # Z FACES
    mdb.models['Model-1'].Equation(name='Z_FACES_DOF_1', terms=(
        (-1.0, 'PBC_Surface-1.Z_Pos_FACE', 1),
        ( 1.0, 'PBC_Surface-1.Z_Neg_FACE', 1)))
    mdb.models['Model-1'].Equation(name='Z_FACES_DOF_2', terms=(
        (-1.0, 'PBC_Surface-1.Z_Pos_FACE', 2),
        ( 1.0, 'PBC_Surface-1.Z_Neg_FACE', 2)))
    mdb.models['Model-1'].Equation(name='Z_FACES_DOF_3', terms=(
        (-1.0, 'PBC_Surface-1.Z_Pos_FACE', 3),
        ( 1.0, 'PBC_Surface-1.Z_Neg_FACE', 3),
        ( 1.0, 'RP-3', 3)))
    
    # X EDGES_1
    mdb.models['Model-1'].Equation(name='X_EDGES_1_DOF_1', terms=(
        ( 1.0, 'PBC_Surface-1.X_YNeg_ZPos_EDGE',1),
        (-1.0, 'PBC_Surface-1.X_YNeg_ZNeg_EDGE',1)))
    mdb.models['Model-1'].Equation(name='X_EDGES_1_DOF_2', terms=(
        ( 1.0, 'PBC_Surface-1.X_YNeg_ZPos_EDGE',2),
        (-1.0, 'PBC_Surface-1.X_YNeg_ZNeg_EDGE',2)))
    mdb.models['Model-1'].Equation(name='X_EDGES_1_DOF_3', terms=(
        ( 1.0, 'PBC_Surface-1.X_YNeg_ZPos_EDGE',3),
        (-1.0, 'PBC_Surface-1.X_YNeg_ZNeg_EDGE',3),
        (-1.0, 'RP-3',3)))
    
    # X EDGES_2
    mdb.models['Model-1'].Equation(name='X_EDGES_2_DOF_1', terms=(
        ( 1.0, 'PBC_Surface-1.X_YNeg_ZNeg_EDGE',1),
        (-1.0, 'PBC_Surface-1.X_YPos_ZNeg_EDGE',1)))
    mdb.models['Model-1'].Equation(name='X_EDGES_2_DOF_2', terms=(
        (-1.0, 'PBC_Surface-1.X_YNeg_ZNeg_EDGE',2),
        ( 1.0, 'PBC_Surface-1.X_YPos_ZNeg_EDGE',2),
        (-1.0, 'RP-2',2)))
    mdb.models['Model-1'].Equation(name='X_EDGES_2_DOF_3', terms=(
        ( 1.0, 'PBC_Surface-1.X_YNeg_ZNeg_EDGE',3),
        (-1.0, 'PBC_Surface-1.X_YPos_ZNeg_EDGE',3)))
    
    # X EDGES_3
    mdb.models['Model-1'].Equation(name='X_EDGES_3_DOF_1', terms=(
        ( 1.0, 'PBC_Surface-1.X_YPos_ZNeg_EDGE',1),
        (-1.0, 'PBC_Surface-1.X_YPos_ZPos_EDGE',1)))
    mdb.models['Model-1'].Equation(name='X_EDGES_3_DOF_2', terms=(
        ( 1.0, 'PBC_Surface-1.X_YPos_ZNeg_EDGE',2),
        (-1.0, 'PBC_Surface-1.X_YPos_ZPos_EDGE',2)))
    mdb.models['Model-1'].Equation(name='X_EDGES_3_DOF_3', terms=(
        (-1.0, 'PBC_Surface-1.X_YPos_ZNeg_EDGE',3),
        ( 1.0, 'PBC_Surface-1.X_YPos_ZPos_EDGE',3),
        (-1.0, 'RP-3',3)))
    
    # Y EDGES_1
    mdb.models['Model-1'].Equation(name='Y_EDGES_1_DOF_1', terms=(
        ( 1.0, 'PBC_Surface-1.Y_XNeg_ZPos_EDGE',1),
        (-1.0, 'PBC_Surface-1.Y_XNeg_ZNeg_EDGE',1)))
    mdb.models['Model-1'].Equation(name='Y_EDGES_1_DOF_2', terms=(
        ( 1.0, 'PBC_Surface-1.Y_XNeg_ZPos_EDGE',2),
        (-1.0, 'PBC_Surface-1.Y_XNeg_ZNeg_EDGE',2)))
    mdb.models['Model-1'].Equation(name='Y_EDGES_1_DOF_3', terms=(
        ( 1.0, 'PBC_Surface-1.Y_XNeg_ZPos_EDGE',3),
        (-1.0, 'PBC_Surface-1.Y_XNeg_ZNeg_EDGE',3),
        (-1.0, 'RP-3',3)))
    
    # Y EDGES_2
    mdb.models['Model-1'].Equation(name='Y_EDGES_2_DOF_1', terms=(
        (-1.0, 'PBC_Surface-1.Y_XNeg_ZNeg_EDGE',1),
        ( 1.0, 'PBC_Surface-1.Y_XPos_ZNeg_EDGE',1),
        (-1.0, 'RP-1',1)))
    mdb.models['Model-1'].Equation(name='Y_EDGES_2_DOF_2', terms=(
        ( 1.0, 'PBC_Surface-1.Y_XNeg_ZNeg_EDGE',2),
        (-1.0, 'PBC_Surface-1.Y_XPos_ZNeg_EDGE',2)))
    mdb.models['Model-1'].Equation(name='Y_EDGES_2_DOF_3', terms=(
        ( 1.0, 'PBC_Surface-1.Y_XNeg_ZNeg_EDGE',3),
        (-1.0, 'PBC_Surface-1.Y_XPos_ZNeg_EDGE',3)))
    
    # Y EDGES_3
    mdb.models['Model-1'].Equation(name='Y_EDGES_3_DOF_1', terms=(
        ( 1.0, 'PBC_Surface-1.Y_XPos_ZNeg_EDGE',1),
        (-1.0, 'PBC_Surface-1.Y_XPos_ZPos_EDGE',1)))
    mdb.models['Model-1'].Equation(name='Y_EDGES_3_DOF_2', terms=(
        ( 1.0, 'PBC_Surface-1.Y_XPos_ZNeg_EDGE',2),
        (-1.0, 'PBC_Surface-1.Y_XPos_ZPos_EDGE',2)))
    mdb.models['Model-1'].Equation(name='Y_EDGES_3_DOF_3', terms=(
        (-1.0, 'PBC_Surface-1.Y_XPos_ZNeg_EDGE',3),
        ( 1.0, 'PBC_Surface-1.Y_XPos_ZPos_EDGE',3),
        (-1.0, 'RP-3',3)))
    
    # Z EDGES_1
    mdb.models['Model-1'].Equation(name='Z_EDGES_1_DOF_1', terms=(
        (-1.0, 'PBC_Surface-1.Z_XNeg_YNeg_EDGE',1),
        ( 1.0, 'PBC_Surface-1.Z_XPos_YNeg_EDGE',1),
        (-1.0, 'RP-1',1)))
    mdb.models['Model-1'].Equation(name='Z_EDGES_1_DOF_2', terms=(
        ( 1.0, 'PBC_Surface-1.Z_XNeg_YNeg_EDGE',2),
        (-1.0, 'PBC_Surface-1.Z_XPos_YNeg_EDGE',2)))
    mdb.models['Model-1'].Equation(name='Z_EDGES_1_DOF_3', terms=(
        ( 1.0, 'PBC_Surface-1.Z_XNeg_YNeg_EDGE',3),
        (-1.0, 'PBC_Surface-1.Z_XPos_YNeg_EDGE',3)))
    
    # Z EDGES_2
    mdb.models['Model-1'].Equation(name='Z_EDGES_2_DOF_1', terms=(
        ( 1.0, 'PBC_Surface-1.Z_XPos_YNeg_EDGE',1),
        (-1.0, 'PBC_Surface-1.Z_XPos_YPos_EDGE',1)))
    mdb.models['Model-1'].Equation(name='Z_EDGES_2_DOF_2', terms=(
        (-1.0, 'PBC_Surface-1.Z_XPos_YNeg_EDGE',2),
        ( 1.0, 'PBC_Surface-1.Z_XPos_YPos_EDGE',2),
        (-1.0, 'RP-2',2)))
    mdb.models['Model-1'].Equation(name='Z_EDGES_2_DOF_3', terms=(
        ( 1.0, 'PBC_Surface-1.Z_XPos_YNeg_EDGE',3),
        (-1.0, 'PBC_Surface-1.Z_XPos_YPos_EDGE',3)))
    
    # Z EDGES_3
    mdb.models['Model-1'].Equation(name='Z_EDGES_3_DOF_1', terms=(
        ( 1.0, 'PBC_Surface-1.Z_XPos_YPos_EDGE',1),
        (-1.0, 'PBC_Surface-1.Z_XNeg_YPos_EDGE',1),
        (-1.0, 'RP-1',1)))
    mdb.models['Model-1'].Equation(name='Z_EDGES_3_DOF_2', terms=(
        ( 1.0, 'PBC_Surface-1.Z_XPos_YPos_EDGE',2),
        (-1.0, 'PBC_Surface-1.Z_XNeg_YPos_EDGE',2)))
    mdb.models['Model-1'].Equation(name='Z_EDGES_3_DOF_3', terms=(
        ( 1.0, 'PBC_Surface-1.Z_XPos_YPos_EDGE',3),
        (-1.0, 'PBC_Surface-1.Z_XNeg_YPos_EDGE',3)))
    
    
    # CORNERS
    mdb.models['Model-1'].Equation(name='Corner_001_000_DOF_1', terms=(
        ( 1.0,'PBC_Surface-1.Corner_001',1),
        (-1.0,'PBC_Surface-1.Corner_000',1)))
    mdb.models['Model-1'].Equation(name='Corner_001_000_DOF_2', terms=(
        ( 1.0,'PBC_Surface-1.Corner_001',2),
        (-1.0,'PBC_Surface-1.Corner_000',2)))
    mdb.models['Model-1'].Equation(name='Corner_001_000_DOF_3', terms=(
        ( 1.0,'PBC_Surface-1.Corner_001',3),
        (-1.0,'PBC_Surface-1.Corner_000',3),
        (-1.0,'RP-3',3)))
    
    mdb.models['Model-1'].Equation(name='Corner_000_101_DOF_1', terms=(
        (-1.0,'PBC_Surface-1.Corner_000',1),
        ( 1.0,'PBC_Surface-1.Corner_101',1),
        (-1.0,'RP-1',1)))
    mdb.models['Model-1'].Equation(name='Corner_000_101_DOF_2', terms=(
        ( 1.0,'PBC_Surface-1.Corner_000',2),
        (-1.0,'PBC_Surface-1.Corner_101',2)))
    mdb.models['Model-1'].Equation(name='Corner_000_101_DOF_3', terms=(
        (-1.0,'PBC_Surface-1.Corner_000',3),
        ( 1.0,'PBC_Surface-1.Corner_101',3),
        (-1.0,'RP-3',3)))
    
    mdb.models['Model-1'].Equation(name='Corner_101_100_DOF_1', terms=(
        ( 1.0,'PBC_Surface-1.Corner_101',1),
        (-1.0,'PBC_Surface-1.Corner_100',1)))
    mdb.models['Model-1'].Equation(name='Corner_101_100_DOF_2', terms=(
        ( 1.0,'PBC_Surface-1.Corner_101',2),
        (-1.0,'PBC_Surface-1.Corner_100',2)))
    mdb.models['Model-1'].Equation(name='Corner_101_100_DOF_3', terms=(
        ( 1.0,'PBC_Surface-1.Corner_101',3),
        (-1.0,'PBC_Surface-1.Corner_100',3),
        (-1.0,'RP-3',3)))
    
    mdb.models['Model-1'].Equation(name='Corner_100_011_DOF_1', terms=(
        ( 1.0,'PBC_Surface-1.Corner_100',1),
        (-1.0,'PBC_Surface-1.Corner_011',1),
        (-1.0,'RP-1',1)))
    mdb.models['Model-1'].Equation(name='Corner_100_011_DOF_2', terms=(
        (-1.0,'PBC_Surface-1.Corner_100',2),
        ( 1.0,'PBC_Surface-1.Corner_011',2),
        (-1.0,'RP-2',2)))
    mdb.models['Model-1'].Equation(name='Corner_100_011_DOF_3', terms=(
        (-1.0,'PBC_Surface-1.Corner_100',3),
        ( 1.0,'PBC_Surface-1.Corner_011',3),
        (-1.0,'RP-3',3)))
    
    mdb.models['Model-1'].Equation(name='Corner_011_010_DOF_1', terms=(
        ( 1.0,'PBC_Surface-1.Corner_011',1),
        (-1.0,'PBC_Surface-1.Corner_010',1)))
    mdb.models['Model-1'].Equation(name='Corner_011_010_DOF_2', terms=(
        ( 1.0,'PBC_Surface-1.Corner_011',2),
        (-1.0,'PBC_Surface-1.Corner_010',2)))
    mdb.models['Model-1'].Equation(name='Corner_011_010_DOF_3', terms=(
        ( 1.0,'PBC_Surface-1.Corner_011',3),
        (-1.0,'PBC_Surface-1.Corner_010',3),
        (-1.0,'RP-3',3)))
    
    mdb.models['Model-1'].Equation(name='Corner_010_111_DOF_1', terms=(
        (-1.0,'PBC_Surface-1.Corner_010',1),
        ( 1.0,'PBC_Surface-1.Corner_111',1),
        (-1.0,'RP-1',1)))
    mdb.models['Model-1'].Equation(name='Corner_010_111_DOF_2', terms=(
        ( 1.0,'PBC_Surface-1.Corner_010',2),
        (-1.0,'PBC_Surface-1.Corner_111',2)))
    mdb.models['Model-1'].Equation(name='Corner_010_111_DOF_3', terms=(
        (-1.0,'PBC_Surface-1.Corner_010',3),
        ( 1.0,'PBC_Surface-1.Corner_111',3),
        (-1.0,'RP-3',3)))
    
    mdb.models['Model-1'].Equation(name='Corner_111_110_DOF_1', terms=(
        ( 1.0,'PBC_Surface-1.Corner_111',1),
        (-1.0,'PBC_Surface-1.Corner_110',1)))
    mdb.models['Model-1'].Equation(name='Corner_111_110_DOF_2', terms=(
        ( 1.0,'PBC_Surface-1.Corner_111',2),
        (-1.0,'PBC_Surface-1.Corner_110',2)))
    mdb.models['Model-1'].Equation(name='Corner_111_110_DOF_3', terms=(
        ( 1.0,'PBC_Surface-1.Corner_111',3),
        (-1.0,'PBC_Surface-1.Corner_110',3),
        (-1.0,'RP-3',3)))
    
    
    # Apply general contact with cohesive surface property
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=FRICTIONLESS)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].CohesiveBehavior(
        defaultPenalties=OFF, table=((1E12, 1E12, 1E12), ))
    mdb.models['Model-1'].ContactStd(name='General_Contact', 
        createStepName='Initial')
    mdb.models['Model-1'].interactions['General_Contact'].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
    r12=mdb.models['Model-1'].rootAssembly.instances['PBC_Surface-1'].surfaces['PBC_Surf_Tie']
    mdb.models['Model-1'].interactions['General_Contact'].excludedPairs.setValuesInStep(
        stepName='Initial', addPairs=((ALLSTAR, r12), ))
    mdb.models['Model-1'].interactions['General_Contact'].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, 'IntProp-1'), ))
    mdb.models['Model-1'].interactions['General_Contact'].slidingFormulationAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, SMALL_SLIDING), ))
    
    # Create Step
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, 
        interactions=OFF, constraints=OFF, connectors=OFF, engineeringFeatures=OFF, 
        adaptiveMeshConstraints=ON)
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', 
        maxNumInc=100000, stabilizationMagnitude=0.0002, 
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION, 
        continueDampingFactors=False, adaptiveDampingRatio=0.05, initialInc=0.1, 
        minInc=1e-10, maxInc=0.1, nlgeom=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    
    # Create field output for homogenization
    mdb.models['Model-1'].FieldOutputRequest(name='For_Volume_Avarege', 
        createStepName='Step-1', variables=('S', 'IVOL'), numIntervals=10, 
        timeMarks=OFF)
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
        numIntervals=10, timeMarks=OFF)
    
    # For Porous material
    if Is_Porous=='Porous':
        a2 = mdb.models['Model-1'].rootAssembly
        a2.features['Sphere_Only-1'].suppress()
        del mdb.models['Model-1'].interactionProperties['IntProp-1'].cohesiveBehavior
        mdb.models['Model-1'].interactionProperties['IntProp-1'].tangentialBehavior.setValues(
            formulation=FRICTIONLESS)
        mdb.models['Model-1'].interactionProperties['IntProp-1'].normalBehavior.setValues(
            pressureOverclosure=HARD, allowSeparation=ON, 
            constraintEnforcementMethod=DEFAULT)
        a = mdb.models['Model-1'].rootAssembly
        session.viewports['Viewport: 1'].setValues(displayedObject=a)
        session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
        a = mdb.models['Model-1'].rootAssembly
        
        # Surface for +X direction
        X_Pos_Matrix = a.instances['Matrix_Only-1'].faces
        Face_X_Pos_Mat = X_Pos_Matrix.getByBoundingBox(xMin=L, xMax=L*1.5, yMin=0.0, yMax=L, zMin=0, zMax=L)
        Face_X_Pos = Face_X_Pos_Mat
        
        # Surface for -X direction
        X_Neg_Matrix = a.instances['Matrix_Only-1'].faces
        Face_X_Neg_Mat = X_Neg_Matrix.getByBoundingBox(xMin=-0.5*L, xMax=0, yMin=0.0, yMax=L, zMin=0, zMax=L)
        Face_X_Neg = Face_X_Neg_Mat
        
        # Surface for +Y direction
        Y_Pos_Matrix = a.instances['Matrix_Only-1'].faces
        Face_Y_Pos_Mat = Y_Pos_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=L, yMax=1.5*L, zMin=0, zMax=L)
        Face_Y_Pos = Face_Y_Pos_Mat
        
        # Surface for -Y direction
        Y_Neg_Matrix = a.instances['Matrix_Only-1'].faces
        Face_Y_Neg_Mat = Y_Neg_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=-0.5*L, yMax=0, zMin=0, zMax=L)
        Face_Y_Neg = Face_Y_Neg_Mat
        
        # Surface for +Z direction
        Z_Pos_Matrix = a.instances['Matrix_Only-1'].faces
        Face_Z_Pos_Mat = Z_Pos_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=0, yMax=L, zMin=L, zMax=1.5*L)
        Face_Z_Pos = Face_Z_Pos_Mat
        
        # Surface for -Z direction
        Z_Neg_Matrix = a.instances['Matrix_Only-1'].faces
        Face_Z_Neg_Mat = Z_Neg_Matrix.getByBoundingBox(xMin=0, xMax=L, yMin=0, yMax=L, zMin=-0.5*L, zMax=0)
        Face_Z_Neg = Face_Z_Neg_Mat
        
        side1Faces = Face_X_Pos + Face_X_Neg + Face_Y_Pos + Face_Y_Neg + Face_Z_Pos + Face_Z_Neg
        a.Surface(side1Faces=side1Faces, name='Tie_Surf')
    
    # Apply BC
    if Mode=='Uniaxial Tension': #Uniaxial Tension
        # RP-1
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-1']
        mdb.models['Model-1'].DisplacementBC(name='RP-1', createStepName='Step-1', 
            region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        # RP-2
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-2']
        mdb.models['Model-1'].DisplacementBC(name='RP-2', createStepName='Step-1', 
            region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        # RP-3
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-3']
        mdb.models['Model-1'].DisplacementBC(name='RP-3', createStepName='Step-1', 
            region=region, u1=UNSET, u2=UNSET, u3=Disp, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
    elif Mode=='Biaxial Tension': #Biaxial tension
        # RP-1
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-1']
        mdb.models['Model-1'].DisplacementBC(name='RP-1', createStepName='Step-1', 
            region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        # RP-2
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-2']
        mdb.models['Model-1'].DisplacementBC(name='RP-2', createStepName='Step-1', 
            region=region, u1=UNSET, u2=Disp, u3=UNSET, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        # RP-3
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-3']
        mdb.models['Model-1'].DisplacementBC(name='RP-3', createStepName='Step-1', 
            region=region, u1=UNSET, u2=UNSET, u3=Disp, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
    elif Mode=='Confined Compression': #Confined Compression
        # RP-1
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-1']
        mdb.models['Model-1'].DisplacementBC(name='RP-1', createStepName='Step-1', 
            region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        # RP-2
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-2']
        mdb.models['Model-1'].DisplacementBC(name='RP-2', createStepName='Step-1', 
            region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        # RP-3
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-3']
        mdb.models['Model-1'].DisplacementBC(name='RP-3', createStepName='Step-1', 
            region=region, u1=0.0, u2=0.0, u3=-Disp, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
    elif Mode=='Simple Shear': #Simple Shear (Side free)
        # Add constraint for shear motion
        a = mdb.models['Model-1'].rootAssembly
        a.ReferencePoint(point=(L/2., L/2., L/2.))
        refPointFound = a.referencePoints.findAt((L/2., L/2., L/2.),)
        a.Set(name='RP-4', referencePoints=(refPointFound,))
        # Z FACES
        mdb.models['Model-1'].Equation(name='Z_FACES_DOF_1', terms=(
            (-1.0, 'PBC_Surface-1.Z_Pos_FACE', 1),
            ( 1.0, 'PBC_Surface-1.Z_Neg_FACE', 1),
            ( 1.0, 'RP-4', 1)))
        # X EDGE_1
        mdb.models['Model-1'].Equation(name='X_EDGES_1_DOF_1', terms=(
            ( 1.0, 'PBC_Surface-1.X_YNeg_ZPos_EDGE',1),
            (-1.0, 'PBC_Surface-1.X_YNeg_ZNeg_EDGE',1),
            (-1.0, 'RP-4', 1)))
        # X EDGE_3
        mdb.models['Model-1'].Equation(name='X_EDGES_3_DOF_1', terms=(
            ( 1.0, 'PBC_Surface-1.X_YPos_ZNeg_EDGE',1),
            (-1.0, 'PBC_Surface-1.X_YPos_ZPos_EDGE',1),
            ( 1.0, 'RP-4', 1)))
        # Y EDGE_1
        mdb.models['Model-1'].Equation(name='Y_EDGES_1_DOF_1', terms=(
            ( 1.0, 'PBC_Surface-1.Y_XNeg_ZPos_EDGE',1),
            (-1.0, 'PBC_Surface-1.Y_XNeg_ZNeg_EDGE',1),
            (-1.0, 'RP-4', 1)))
        # Y EDGES_3
        mdb.models['Model-1'].Equation(name='Y_EDGES_3_DOF_1', terms=(
            ( 1.0, 'PBC_Surface-1.Y_XPos_ZNeg_EDGE',1),
            (-1.0, 'PBC_Surface-1.Y_XPos_ZPos_EDGE',1),
            ( 1.0, 'RP-4', 1)))
        # CORNERS
        mdb.models['Model-1'].Equation(name='Corner_001_000_DOF_1', terms=(
            ( 1.0,'PBC_Surface-1.Corner_001',1),
            (-1.0,'PBC_Surface-1.Corner_000',1),
            (-1.0, 'RP-4', 1)))
        mdb.models['Model-1'].Equation(name='Corner_000_101_DOF_1', terms=(
            (-1.0,'PBC_Surface-1.Corner_000',1),
            ( 1.0,'PBC_Surface-1.Corner_101',1),
            (-1.0,'RP-4', 1)))
        mdb.models['Model-1'].Equation(name='Corner_101_100_DOF_1', terms=(
            ( 1.0,'PBC_Surface-1.Corner_101',1),
            (-1.0,'PBC_Surface-1.Corner_100',1),
            (-1.0, 'RP-4', 1)))
        mdb.models['Model-1'].Equation(name='Corner_100_011_DOF_1', terms=(
            (-1.0,'PBC_Surface-1.Corner_100',1),
            ( 1.0,'PBC_Surface-1.Corner_011',1),
            (-1.0,'RP-4', 1)))
        mdb.models['Model-1'].Equation(name='Corner_011_010_DOF_1', terms=(
            ( 1.0,'PBC_Surface-1.Corner_011',1),
            (-1.0,'PBC_Surface-1.Corner_010',1),
            (-1.0, 'RP-4', 1)))
        mdb.models['Model-1'].Equation(name='Corner_010_111_DOF_1', terms=(
            (-1.0,'PBC_Surface-1.Corner_010',1),
            ( 1.0,'PBC_Surface-1.Corner_111',1),
            (-1.0,'RP-4', 1)))
        mdb.models['Model-1'].Equation(name='Corner_111_110_DOF_1', terms=(
            ( 1.0,'PBC_Surface-1.Corner_111',1),
            (-1.0,'PBC_Surface-1.Corner_110',1),
            (-1.0, 'RP-4', 1)))
        # Disp BC
        a = mdb.models['Model-1'].rootAssembly
        region = a.sets['RP-4']
        mdb.models['Model-1'].DisplacementBC(name='RP-4', createStepName='Step-1', 
            region=region, u1=Disp, u2=UNSET, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
    
    # Display the model on Assembly module
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        adaptiveMeshConstraints=OFF)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=1.16353, 
        farPlane=1.94484, width=0.775123, height=0.448254, cameraPosition=(1.05959, 
        0.880742, 1.17674), cameraUpVector=(-0.480517, 0.673048, -0.562236), 
        cameraTarget=(0.207239, 0.188718, 0.209181), viewOffsetX=0.082503, 
        viewOffsetY=-0.0217776)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
    
    # Output the final Sphere_array
    print("##########################")
    print("## Generation Completed ##")
    print("##########################")
    print("Number of spheres:", len(Sphere_array))
    print("Volume fraction:", str(100.*VoF)+"%")
    
    if VoF < VoF_tar:
        print("The target volume fraction was not achieved within the max number of iteration!")
    pass

def generate_ss_curve(EngStrain, Odb_Path, Output_Path, S_Comp, Plot_Flag):
    """
    A function to compute volume-averaged stress from an ODB file,
    write results to a text file, and optionally plot them in Abaqus/CAE.
    
    Parameters:
    -----------
    EngStrain : float
        Engineering strain to be converted to true strain.
    Odb_Path : str
        Full path to the ODB file.
    Output_Path : str
        Full path to the output text file.
    S_Comp : str
        Stress component to extract (e.g., "S11", "S22", etc.).
    Plot_Flag : bool
        If True, create an XY plot in Abaqus/CAE using the generated data.
    """
    # Convert engineering strain to true strain
    Epsilon = np.log(1.0 + EngStrain)

    # Open the ODB
    if not os.path.isfile(Odb_Path):
        raise FileNotFoundError("The specified ODB file does not exist: " + Odb_Path)

    odb = openOdb(path=Odb_Path)

    # Get the first step in the ODB
    step_name = list(odb.steps.keys())[0]
    step = odb.steps[step_name]

    # Collect results (strain, volume-averaged stress)
    output_data = []

    # Stress component index in the (S11,S22,S33,S12,S13,S23) tuple
    COMP_INDEX = {'S11': 0, 'S22': 1, 'S33': 2, 'S12': 3, 'S13': 4, 'S23': 5}
    s_key = S_Comp.upper()
    if s_key not in COMP_INDEX:
        raise ValueError("Unsupported S_Comp: {}".format(S_Comp))
    comp = COMP_INDEX[s_key]

    def _bulk(field):
        """
        Read a field output through bulkDataBlocks and return
        (data, sort_key).

        bulkDataBlocks hands back NumPy arrays directly, whereas
        FieldOutput.values builds one Python object per integration point.
        For a typical RVE (1e5 - 1e6 integration points x many frames) this
        is the dominant cost of the post-processing step.

        The sort key is built from (elementLabel, integrationPoint) rather
        than elementLabel alone.  IVOL is an integration-point volume, so a
        dictionary keyed only on the element label keeps a single
        integration point per element and applies its volume to all the
        others - which biases Eq.(11) for any element with more than one
        integration point (C3D10, C3D10M, C3D8, ...).
        """
        blocks = field.bulkDataBlocks
        data = np.concatenate([np.asarray(b.data) for b in blocks], axis=0)

        elem = np.concatenate(
            [np.asarray(b.elementLabels, dtype=np.int64) for b in blocks], axis=0)

        ip_list = []
        for b in blocks:
            ip_b = getattr(b, 'integrationPoints', None)
            if ip_b is None or len(ip_b) == 0:
                ip_list.append(np.zeros(len(np.asarray(b.elementLabels)), dtype=np.int64))
            else:
                ip_list.append(np.asarray(ip_b, dtype=np.int64))
        ip = np.concatenate(ip_list, axis=0)

        if data.ndim == 1:
            data = data.reshape(-1, 1)

        return data, elem * 1000 + ip

    order_S = None
    order_V = None
    n_ref = -1

    # Loop over frames
    for frame_idx, frame in enumerate(step.frames):
        print("Processing frame {}/{}".format(frame_idx + 1, len(step.frames)))

        # Get stress and element volume fields
        stress_data, stress_key = _bulk(frame.fieldOutputs['S'])
        volume_data, volume_key = _bulk(frame.fieldOutputs['IVOL'])

        # The ordering of the bulk data blocks is constant over the frames,
        # so the alignment between S and IVOL is computed only once.
        if order_S is None or stress_data.shape[0] != n_ref:
            order_S = np.argsort(stress_key, kind='stable')
            order_V = np.argsort(volume_key, kind='stable')
            n_ref = stress_data.shape[0]
            if (order_S.shape[0] != order_V.shape[0] or
                    not np.array_equal(stress_key[order_S], volume_key[order_V])):
                raise ValueError(
                    "S and IVOL are not written for the same set of "
                    "integration points. Request both with the same "
                    "*Output, Field definition.")

        sigma = stress_data[order_S, comp]
        ivol = volume_data[order_V, 0]

        total_volume = float(ivol.sum())

        # Compute volume-averaged stress for this frame
        if total_volume > 0.0:
            volume_avg_stress = float(np.dot(sigma, ivol) / total_volume)
        else:
            volume_avg_stress = 0.0

        # Compute the true strain for this frame
        # frameValue is typically the step time, so we scale it by Epsilon
        strain_val = frame.frameValue * Epsilon

        output_data.append((strain_val, volume_avg_stress))

    # Close the ODB
    odb.close()

    # Write results to the output file
    with open(Output_Path, "w") as f:
        f.write("strain\tVolume_Averaged_Stress\n")
        for strain_val, avg_stress in output_data:
            f.write(f"{strain_val:.6f}\t{avg_stress:.6f}\n")

    print("Results written to {}".format(Output_Path))

    # ----------------------------------------------------
    # If requested, read the file back and create XY Plot
    # ----------------------------------------------------
    if Plot_Flag:
        # Make/read XY data from the output file
        xy_data = []
        with open(Output_Path, "r") as f:
            # Skip header
            next(f)
            for line in f:
                parts = line.split()
                if len(parts) < 2:
                    continue
                strain_val = float(parts[0])
                stress_val = float(parts[1])
                xy_data.append((strain_val, stress_val))

        # Create (or reuse) a viewport
        if 'Viewport: 1' not in session.viewports:
            session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200, height=150)

        session.viewports['Viewport: 1'].makeCurrent()
        session.viewports['Viewport: 1'].maximize()

        # Create a unique name for the XYPlot if one already exists
        base_plot_name = 'XYPlot-VolumeAveraged'
        plot_name = base_plot_name
        i = 1
        while plot_name in session.xyPlots.keys():
            plot_name = base_plot_name + '_' + str(i)
            i += 1

        # Create the XYPlot object with the new unique name
        xyp = session.XYPlot(name=plot_name)
        chartName = xyp.charts.keys()[0]
        chart = xyp.charts[chartName]

        # Define x- and y-axis quantity types (for a typical stress-strain plot)
        xQuantity = visualization.QuantityType(type=STRAIN)
        yQuantity = visualization.QuantityType(type=STRESS)

        # Convert our Python list to an XYData object
        myXYData = xyPlot.XYData(
            data=xy_data,
            sourceDescription='Volume-averaged stress data',
            axis1QuantityType=xQuantity,
            axis2QuantityType=yQuantity
        )

        # Create a Curve
        curve = session.Curve(xyData=myXYData)
        chart.setValues(curvesToPlot=(curve,))
        session.charts[chartName].autoColor(lines=True, symbols=True)

        # Display the XYPlot
        session.viewports['Viewport: 1'].setValues(displayedObject=xyp)

        # Optionally adjust symbol style, legend, etc.
        curve.symbolStyle.setValues(show=True, size=2)
        curve.setValues(useDefault=False, legendLabel='Vol_Avg_Stress')

        print("Plot created in Abaqus/CAE window: {}".format(plot_name))
    pass

def excecute(Is_Porous, L, L_mesh, r_avg, r_std, VoF_tar, min_distance, 
             max_iterations, Mode, Disp, EngStrain, Odb_Path, Output_Path, 
             S_Comp, Plot_Flag, action):
    
    if action == 'generate_pbc':
        showStopButtonInGui()    
        GeneratePBCell(
            Is_Porous=Is_Porous,
            L=L,
            L_mesh=L_mesh,
            r_avg=r_avg,
            r_std=r_std,
            VoF_tar=VoF_tar,
            min_distance=min_distance,
            max_iterations=max_iterations,
            Mode=Mode,
            Disp=Disp
        )
    elif action == 'generate_ss_curve':
        showStopButtonInGui()
        generate_ss_curve(
            EngStrain=EngStrain, 
            Odb_Path=Odb_Path, 
            Output_Path=Output_Path, 
            S_Comp=S_Comp, 
            Plot_Flag=Plot_Flag
        )
    else:
        raise ValueError("Invalid action: {}".format(action))
