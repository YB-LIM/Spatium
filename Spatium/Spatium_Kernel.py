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
"""

import numpy as np
from itertools import product
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

class VoxelGrid:
    def __init__(self, L, min_distance, resolutionFactor=2):
        """
        :param L: Size of the cell
        :param min_distance: Minimum center-to-center distance between spheres
        :param resolutionFactor: Increase to get finer voxel spacing.
                                (Voxel size = min_distance / resolutionFactor)
        """
        # Define the voxel size
        self.voxel_size = float(min_distance) / float(resolutionFactor)
        # Number of voxels along each dimension
        self.nx = int(np.ceil(L / self.voxel_size))
        self.ny = int(np.ceil(L / self.voxel_size))
        self.nz = int(np.ceil(L / self.voxel_size))
        
        # Occupancy array; False means "free," True means "occupied"
        self.occupied = np.zeros((self.nx, self.ny, self.nz), dtype=bool)
        self.L = L

    def _xyz_to_ijk(self, x, y, z):
        """
        Convert continuous coordinates (x,y,z) in [0,L] to integer voxel indices (i,j,k).
        We'll clamp to [0, nx-1], etc., just for safety in edge cases.
        """
        i = int(np.floor(x / self.voxel_size))
        j = int(np.floor(y / self.voxel_size))
        k = int(np.floor(z / self.voxel_size))
        i = max(0, min(self.nx - 1, i))
        j = max(0, min(self.ny - 1, j))
        k = max(0, min(self.nz - 1, k))
        return i, j, k

    def is_occupied(self, center, radius):
        """
        Check if the sphere with given (center, radius) intersects any already-occupied voxel.
        We do a bounding-box over the voxel grid and see if any voxel in that region is True.
        """
        # We have to handle periodic offsets in x,y,z
        cx, cy, cz = center
        # Because we do periodic boundary checks in the main code,
        # we only consider direct bounding box in [0, L].
        # If center is outside [0,L], wrap it in [0,L] for voxel check:
        cx_period = cx % self.L
        cy_period = cy % self.L
        cz_period = cz % self.L
        
        # bounding box in voxel coords
        r_margin = radius + 0.5 * self.voxel_size
        x_min = max(0.0, cx_period - r_margin)
        x_max = min(self.L, cx_period + r_margin)
        y_min = max(0.0, cy_period - r_margin)
        y_max = min(self.L, cy_period + r_margin)
        z_min = max(0.0, cz_period - r_margin)
        z_max = min(self.L, cz_period + r_margin)

        i_min, j_min, k_min = self._xyz_to_ijk(x_min, y_min, z_min)
        i_max, j_max, k_max = self._xyz_to_ijk(x_max, y_max, z_max)

        # Loop through that voxel region to see if any is True
        subarray = self.occupied[i_min:i_max+1, j_min:j_max+1, k_min:k_max+1]
        return np.any(subarray)  # True if any voxel is occupied

    def mark_voxels_for_sphere(voxel_grid, sphere_center, sphere_radius, L, voxel_size):
        """
        Marks voxels as occupied for the sphere and its periodic images.
    
        Parameters:
        - voxel_grid: 3D numpy array representing the voxel grid.
        - sphere_center: Center of the sphere (x, y, z).
        - sphere_radius: Radius of the sphere.
        - L: Size of the domain.
        - voxel_size: Size of each voxel.
        """
        shifts = [-L, 0, L]  # For periodic images
        x, y, z = sphere_center
        r = sphere_radius
    
        for dx in shifts:
            for dy in shifts:
                for dz in shifts:
                    # Adjust for periodic images
                    x_new = x + dx
                    y_new = y + dy
                    z_new = z + dz
    
                    # Compute voxel indices for the bounding box of the sphere
                    x_min = int(max(0, (x_new - r) // voxel_size))
                    x_max = int(min(voxel_grid.shape[0] - 1, (x_new + r) // voxel_size))
                    y_min = int(max(0, (y_new - r) // voxel_size))
                    y_max = int(min(voxel_grid.shape[1] - 1, (y_new + r) // voxel_size))
                    z_min = int(max(0, (z_new - r) // voxel_size))
                    z_max = int(min(voxel_grid.shape[2] - 1, (z_new + r) // voxel_size))
    
                    # Mark voxels as occupied
                    voxel_grid[x_min:x_max + 1, y_min:y_max + 1, z_min:z_max + 1] = True

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


class Octree:
    """A simple wrapper around the root node."""
    def __init__(self, L, capacity=4, max_depth=10):
        self.root = OctreeNode(
            x_min=-0.2*L, x_max=1.2*L,
            y_min=-0.2*L, y_max=1.2*L,
            z_min=-0.2*L, z_max=1.2*L,
            capacity=capacity,
            max_depth=max_depth
        )

    def insert_sphere(self, center, radius):
        self.root.insert(center, radius)

    def query_spheres(self, center, radius):
        return self.root.query(center, radius)

# Create a global (or external) voxelGrid reference. If you prefer, you can
# make it a class variable or pass it around. As long as we don't modify
# GeneratePBCell, we do it here.
voxelGrid = None  # Will be assigned before calling GeneratePBCell

def is_overlapping(Sphere_array_unused, center, radius, L, min_distance, octree):
    global voxelGrid
    if voxelGrid is not None:
        # First, check if the bounding volume is occupied in the voxel grid
        if voxelGrid.is_occupied(center, radius):
            return True

    # Query the octree for nearby spheres
    search_radius = radius + min_distance + L/8.0
    nearby_spheres = octree.query_spheres(center, search_radius)

    for (center2, radius2) in nearby_spheres:
        # Minimum image logic for periodic boundaries
        dx = center[0] - center2[0]
        dx -= L * np.round(dx / L)
        dy = center[1] - center2[1]
        dy -= L * np.round(dy / L)
        dz = center[2] - center2[2]
        dz -= L * np.round(dz / L)

        dist2 = dx*dx + dy*dy + dz*dz
        radii_sum = radius + radius2 + min_distance
        if dist2 < radii_sum * radii_sum:
            return True

    # If we got here => no overlap with any existing spheres
    # Mark the voxel space as occupied so no subsequent sphere tries
    # to place in the same region.
    if voxelGrid is not None:
        voxelGrid.mark_occupied(center, radius)

    return False

def GeneratePBCell(Is_Porous, L, L_mesh, r_avg, r_std, VoF_tar, min_distance, 
                   max_iterations, Mode, Disp):
    #################################################################
    # Generate numpy array that contains center position and radius #
    #################################################################
    # Show stop button
    showStopButtonInGui()
    VoF = 0.0
    Sphere_array = np.empty((0, 4))

    # -------------------------------------------------------
    # Initialize the Octree for overlap checks
    # -------------------------------------------------------
    # Increase capacity or max_depth if many spheres
    octree = Octree(L, capacity=8, max_depth=10)

    iterations = 0

    while VoF <= VoF_tar and iterations < max_iterations:
        iterations += 1

        # Generate sphere using random PDF
        center_X = L * np.random.rand()
        center_Y = L * np.random.rand()
        center_Z = L * np.random.rand()
        center = np.array([center_X, center_Y, center_Z])
        radius = np.random.normal(r_avg, r_std)

        # Skip negative or zero radii
        if radius <= 0:
            continue

        # Check if the sphere overlaps with existing spheres using octree
        if is_overlapping(Sphere_array, center, radius, L, min_distance, octree):
            continue  # Discard this sphere and try again

        # Adjust the radius to account for min_distance when checking boundaries
        effective_radius = radius + min_distance / 2

        # Determine if the sphere overlaps with any boundaries
        Is_inside = (
            effective_radius < center_X < L - effective_radius and
            effective_radius < center_Y < L - effective_radius and
            effective_radius < center_Z < L - effective_radius
        )

        # Boundary overlap checks
        Is_X_Left = center_X - effective_radius < 0
        Is_X_Right = center_X + effective_radius > L
        Is_Y_Front = center_Y - effective_radius < 0
        Is_Y_Back = center_Y + effective_radius > L
        Is_Z_Lower = center_Z - effective_radius < 0
        Is_Z_Upper = center_Z + effective_radius > L

        # Function to check overlapping portion
        def check_overlap_portion(center_coord, L_, r_, is_left, is_right):
            if is_left:
                d = center_coord
            elif is_right:
                d = L_ - center_coord
            else:
                return False  # Does not overlap in this direction

            if d <= 0 or r_ - d <= 0:
                return True  # Invalid, exclude
            Tol_ = 0.4
            ratio = d / (r_ - d)
            if ratio < Tol_ or (1 / ratio) < Tol_:
                return True  # Overlapping portion is too small
            return False  # Overlapping portion is acceptable

        overlaps_too_small = False
        if Is_X_Left or Is_X_Right:
            if check_overlap_portion(center_X, L, radius, Is_X_Left, Is_X_Right):
                overlaps_too_small = True

        if Is_Y_Front or Is_Y_Back:
            if check_overlap_portion(center_Y, L, radius, Is_Y_Front, Is_Y_Back):
                overlaps_too_small = True

        if Is_Z_Lower or Is_Z_Upper:
            if check_overlap_portion(center_Z, L, radius, Is_Z_Lower, Is_Z_Upper):
                overlaps_too_small = True

        if overlaps_too_small:
            continue  # Discard this sphere and try again

        # ------------------------------------------------------------------
        # Case A: Sphere is fully inside; just add it & insert into octree
        # ------------------------------------------------------------------
        if Is_inside:
            new_row = np.array([[center_X, center_Y, center_Z, radius]])
            Sphere_array = np.vstack((Sphere_array, new_row))
            VoF += (4.0 / 3.0) * np.pi * (radius ** 3) / (L ** 3)

            # Insert into octree
            octree.insert_sphere(center, radius)

        else:
            # ------------------------------------------------------------------
            # Case B: Sphere crosses boundary => create mirrored spheres
            # ------------------------------------------------------------------
            from itertools import product
            new_rows = []
            new_rows.append([center_X, center_Y, center_Z, radius])  # original

            # Generate shifts for each axis based on boundary overlaps
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

            # Create mirrored spheres based on shift combinations
            shift_combinations = list(product(x_shifts, y_shifts, z_shifts))
            shift_combinations.remove((0, 0, 0))  # remove original shift
            for dx, dy, dz in shift_combinations:
                x_new = center_X + dx
                y_new = center_Y + dy
                z_new = center_Z + dz
                new_rows.append([x_new, y_new, z_new, radius])

            # Remove duplicates (if the same center was produced more than once)
            unique_new_rows = []
            for row in new_rows:
                if not any(np.allclose(row[:3], r[:3]) for r in unique_new_rows):
                    unique_new_rows.append(row)

            # Check for overlaps with existing spheres in octree
            overlaps = False
            for row in unique_new_rows:
                c_ = np.array(row[:3])
                r_ = row[3]
                if is_overlapping(Sphere_array, c_, r_, L, min_distance, octree):
                    overlaps = True
                    break

            if overlaps:
                continue  # Discard entire set and try again

            # If no overlaps, add them all
            Sphere_array = np.vstack((Sphere_array, unique_new_rows))
            # Increase volume only once for the real sphere
            VoF += (4.0 / 3.0) * np.pi * (radius ** 3) / (L ** 3)

            # Insert all new centers into octree
            for row in unique_new_rows:
                c_ = np.array(row[:3])
                r_ = row[3]
                octree.insert_sphere(c_, r_)

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
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD, 
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
        createStepName='Step-1', variables=('S', 'EVOL'), numIntervals=10, 
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
