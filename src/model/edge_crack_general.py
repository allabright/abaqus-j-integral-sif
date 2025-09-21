from abaqusConstants import *
import regionToolset

from src.model.create_model import CreateValidationModel

class CreateValidationModelEdgeCrack(CreateValidationModel):

    def __init__(self):
        super(CreateValidationModelEdgeCrack, self).__init__()
        
    def define_part_geometry(self):
        """
        Define a dictionary to hold the part geometry, calculating the positions of each corner while keeping
        the crack tip at the origin.
        """

        self.part_geometry = {
            "top_left_corner_coords": (-self.crack_length_mm, self.part_height_mm / 2),
            "top_right_corner_coords": (self.part_width_mm - self.crack_length_mm, self.part_height_mm / 2),
            "bottom_left_corner_coords": (-self.crack_length_mm, -self.part_height_mm / 2),
            "bottom_right_corner_coords": (self.part_width_mm - self.crack_length_mm, -self.part_height_mm / 2),
            "crack_start_coords": (-self.crack_length_mm, self.crack_tip_coordinates[1]),
            "crack_end_coords": self.crack_tip_coordinates}

        return self
    
    def create_part(self):
        """
        Create the sketch for the overall part, and extrude it into a shell or a solid based on the dimension.
        """

        sketch = self.model.ConstrainedSketch(name='Part_Sketch', sheetSize=500.0)
        sketch.rectangle(point1=self.part_geometry["bottom_left_corner_coords"], point2=self.part_geometry["top_right_corner_coords"])

        if self.dimensions == "2d":
            self.part = self.model.Part(name="Part_" + self.analysis_name, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
            self.part.BaseShell(sketch=sketch)

        elif self.dimensions == "3d":
            self.part = self.model.Part(name="Part_" + self.analysis_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
            self.part.BaseSolidExtrude(sketch=sketch, depth=self.part_thickness_mm/2)

        del self.model.sketches['Part_Sketch']

        return self
    
    def calculate_horizontal_partition_heights(self):
        """
        The model is partitioned horizontally in a number of sections defined in the configuration. This method calculates
        the heights of each partition. The partitions may be biased so that they are closer together the closer they are to the
        crack tip.
        """

        base_height = self.crack_length_mm / 2.0
        self.partition_heights = []
        total_ratio = sum(self.horizontal_partition_bias**i for i in range(self.horizontal_partition_count))
        step = (self.part_height_mm / 2.0 - base_height) / total_ratio

        for i in range(self.horizontal_partition_count):
            self.partition_heights.append(base_height + step * sum(self.horizontal_partition_bias**j for j in range(i + 1)))

        return self

    def partition_part_top_to_bottom(self):
        """
        Create the main horizontal partitions for the part, using horizontal lines. This partitions the face rather than the solid.
        A separate method is used to partition the 3D solid once the face has been partitioned.
        """

        line_sketch = self.model.ConstrainedSketch(name='LineSketch', sheetSize=500.0, transform=self.transform)

        line_sketch.Line(point1=(-self.crack_length_mm, 0.0), point2=(self.crack_length_mm, 0.0))

        for partition in self.partition_heights:
            line_sketch.Line(point1=(-self.crack_length_mm, partition), point2=(self.part_width_mm-self.crack_length_mm, partition))
            line_sketch.Line(point1=(-self.crack_length_mm, -partition), point2=(self.part_width_mm-self.crack_length_mm, -partition))

        checkPoint = (0.0, 0.0, self.z)
        faceList = self.part.faces.findAt((checkPoint,))
        if not faceList:
            raise ValueError("No face found at point {}. Check geometry.".format(checkPoint))
        main_face = faceList[0]

        self.part.PartitionFaceBySketch(faces=main_face, sketch=line_sketch)
        del self.model.sketches['LineSketch']

        return self
    
    def partition_edges(self):
        """
        Create partitions for the crack. One partition runs the length of the crack, and the other runs the length of the crack again,
        starting at the crack tip. Therefore we have two line partitions, both starting at the crack tip, and moving in the direction of crack
        growth and opposite (x and -x) for a length equal to the crack length.
        """

        # Create two datum planes - one at the crack tip and one ahead of the crack tip by one crack length.
        datum_plane_crack_tip = self.part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
        datum_plane_forward_of_crack_tip = self.part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=self.crack_length_mm)

        datum_planes = (datum_plane_crack_tip, datum_plane_forward_of_crack_tip)
        x = self.crack_length_mm / 2.0
        y_coords = (-self.partition_heights[0], 0.0, self.partition_heights[0])
        z_coords = (self.tolerance_mm, -self.tolerance_mm)

        for datum_plane in datum_planes:

            entities_to_partition = []

            if self.dimensions == "2d":
                for y in y_coords:
                    entity_to_partition = self.part.edges.findAt(((x, y, 0.0),))[0]
                    entities_to_partition.append(entity_to_partition)
                try:
                    self.part.PartitionEdgeByDatumPlane(edges=entities_to_partition, datumPlane=self.part.datums[datum_plane.id])
                except:
                    continue

            elif self.dimensions == "3d":
                
                for y in y_coords:
                    for z in z_coords:
                        entity_to_partition = self.part.faces.findAt(((x, y, z),))[0]
                        entities_to_partition.append(entity_to_partition)
                try:
                    self.part.PartitionFaceByDatumPlane(faces=entities_to_partition, datumPlane=self.part.datums[datum_plane.id])
                except:
                    continue
        
        return self
    
    def create_sets_for_top_and_bottom_edges(self):
        """
        Create sets for the top and bottom edges / faces, and the bottom left corner, to be used to apply the loading and boundary conditions.
        """

        if self.dimensions == "2d":
            top_edge = self.part.edges.findAt(((0.0, self.part_height_mm / 2, 0.0),))
            bottom_left_corner = self.part.vertices.findAt(((-self.crack_length_mm, -self.part_height_mm / 2, 0.0),))
            bottom_edge = self.part.edges.findAt(((0.0, -self.part_height_mm / 2, 0.0),))
            top_surface = self.part.Surface(name='Top_Surface', side1Edges=top_edge)

            self.part.Set(name='Top_Edge', edges=top_edge)
            self.part.Set(name='Bottom_Left_Corner', vertices=bottom_left_corner)
            self.part.Set(name='Bottom_Edge', edges=bottom_edge)

        elif self.dimensions == "3d":
            top_face = self.part.faces.findAt(((0.0, self.part_height_mm / 2.0, 0.0),))
            bottom_left_edge = self.part.edges.findAt(((-self.crack_length_mm, -self.part_height_mm / 2, self.z / 2.0),))
            bottom_face = self.part.faces.findAt(((0.0, -self.part_height_mm / 2.0, self.z / 2.0),))
            top_surface = self.part.Surface(name='Top_Surface', side1Faces=top_face)

            self.part.Set(name='Top_Face', faces=top_face)
            self.part.Set(name='Bottom_Left_Edge', edges=bottom_left_edge)
            self.part.Set(name='Bottom_Face', faces=bottom_face)

        return self
    
    def create_crack_seam_set_geometry(self):
        """
        Define a crack seam within Abaqus, allowing the crack nodes to separate under load and allow the crack to open.
        """

        x_min, x_max = -self.crack_length_mm - self.tolerance_mm, 0.0 + self.tolerance_mm
        y_min, y_max = -self.tolerance_mm, self.tolerance_mm       

        if self.dimensions == "2d":
            z_min, z_max = -self.tolerance_mm, self.tolerance_mm
            crack_edges = self.part.edges.getByBoundingBox(
                xMin=x_min, xMax=x_max,
                yMin=y_min, yMax=y_max,
                zMin=z_min, zMax=z_max)
            
            self.part.Set(name='CrackSeam', edges=crack_edges)
            crack_region = regionToolset.Region(edges=crack_edges)
            self.part.engineeringFeatures.assignSeam(regions=crack_region)

        elif self.dimensions == "3d":
            z_min, z_max = -self.part_thickness_mm - self.tolerance_mm, self.part_thickness_mm + self.tolerance_mm
            crack_edges = self.instance.faces.getByBoundingBox(
                xMin=x_min, xMax=x_max,
                yMin=y_min, yMax=y_max,
                zMin=z_min, zMax=z_max)
                        
            self.assembly.Set(name='CrackSeam', faces=crack_edges)
            crack_region = regionToolset.Region(faces=crack_edges)
            self.assembly.engineeringFeatures.assignSeam(regions=crack_region)

        return self
        
    def create_crack_seam_set_nodes(self):
        """
        Create a set containing all of the nodes that make up the crack seam.
        """
        
        node_labels = []

        if self.dimensions == "2d":
            geom_set =  self.part.sets['CrackSeam']
            for edge in geom_set.edges:
                nodes = edge.getNodes()
                node_labels.extend([node.label for node in nodes])

        elif self.dimensions == "3d":
            geom_set =  self.assembly.sets['CrackSeam']
            nodes = geom_set.nodes
            self.assembly.Set(name="CrackSeam_Nodes",  nodes=nodes)

        # Create a node set on the part
        if node_labels:
            self.part.SetFromNodeLabels(name='CrackSeam_Nodes', nodeLabels=node_labels)
        else:
            print("No nodes found in 'CrackSeam'!")

        return self

    def apply_top_pressure_load(self):
        """
        Apply the negative pressure load (i.e. tensile stress) to the top face/edge of the part.
        """

        top_region = self.instance.surfaces["Top_Surface"]

        try:
            self.model.Pressure(name='Pressure_Load_Top',
                                createStepName=self.step_name,
                                region=top_region,
                                magnitude=self.applied_stress_mpa,
                                amplitude=UNSET)
        except Exception as e:
            self.print_abaqus("Failed to apply pressure load: " + str(e))

        return self
    
    def apply_bottom_fixed_boundary_condition(self):
        """
        Apply the boundary conditions to fix the lower left corner of the part, and prevent displacement of the lower
        surface in the y-direction (i.e. fixed and roller support respectively).
        """

        if self.dimensions == "2d":
            bottom_left = self.instance.sets['Bottom_Left_Corner']
            bottom_region = self.instance.sets['Bottom_Edge']
        elif self.dimensions == "3d":
            bottom_left = self.instance.sets['Bottom_Left_Edge']
            bottom_region = self.instance.sets['Bottom_Face']
            
        try:
            self.model.DisplacementBC(name='Pinned_Support_Bottom_Left',
                                      createStepName=self.step_name,
                                      region=bottom_left,
                                      u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)

            
            self.model.DisplacementBC(name='Roller_Support_Bottom',
                                      createStepName=self.step_name,
                                      region=bottom_region,
                                      u2=0.0)
        except Exception as e:
            self.print_abaqus("Failed to create boundary condition: " + str(e))

        return self
    
    def seed_crack_tip_region(self):
        """
        Seed the region around the crack tip - effectively four lines starting at the crack tip, in the x, -x, y, and -y directions.
        This allows the mesh to be refined as it gets closer to the crack tip and maximises performance.
        """

        self.left_edges = []
        self.right_edges = []
        self.top_edges = []
        self.bottom_edges = []

        for x, left_right in (
            (-self.crack_length_mm / 2.0, self.left_edges),
            (self.crack_length_mm / 2.0, self.right_edges)):
            for y in (-self.partition_heights[0], 0.0, self.partition_heights[0]):
                for z in (-self.z, 0.0, self.z):
                    edge = self.part.edges.findAt(((x, y, z),))
                    left_right.append(edge[0])

        for y, top_bottom in (
            (self.partition_heights[0] / 2.0, self.top_edges),
            (-self.partition_heights[0] / 2.0, self.bottom_edges)):
            for x in (-self.crack_length_mm, self.crack_length_mm):
                for z in (-self.z, 0.0, self.z):
                    edge = self.part.edges.findAt(((x, y, z),))
                    if edge:
                        top_bottom.append(edge[0])

        for edges in (self.right_edges, self.left_edges, self.top_edges, self.bottom_edges):
            self.part.seedEdgeByBias(SINGLE, end1Edges=edges, ratio=self.crack_element_bias, number=self.crack_element_count)

        return self
    
    def seed_rest_of_part_regular_partitioning(self):
        """
        Seed the rest of the part with a coarse mesh, after excluding the edges near the crack tip that have already had the finer
        mesh seeds applied.
        """

        allEdges = self.part.edges[:]

        other_edges = []

        if self.left_edges and self.right_edges and self.top_edges and self.bottom_edges:
            for edge in allEdges:
                if edge not in self.left_edges and \
                    edge not in self.right_edges and \
                        edge not in self.top_edges and \
                            edge not in self.bottom_edges:
                    
                    if self.dimensions == "2d":
                        other_edges.append(edge)

                    elif self.dimensions == "3d":
                        if edge not in self.through_thickness_edges:
                            other_edges.append(edge)       
        else:
            other_edges = [ed for ed in allEdges]

        if other_edges:
            self.part.seedEdgeBySize(edges=other_edges, size=self.coarse_seed_size_mm, deviationFactor=0.1)

        return self

    def move_crack_tip_nodes(self):
        """
        Move the nodes of the crack tip elements to the quarter-point position. This was only implemented for the
        2D model due to time. Implementation of this for 3D is available in:
        "On the use of quarter-point tetrahedral finite elements in linear elastic fracture mechanics"
        (M. Nejati, A. Paluszny, R. W. Zimmerman, 2015) and can be implemented here in the future.

        """

        crack_tip_element_side_lengths = []
        crack_tip_node = self.part.nodes.getClosest((self.crack_tip_coordinates,))[0]
        crack_tip_elements = self.part.nodes.getFromLabel(crack_tip_node.label).getElements()

        crack_tip_edge_nodes = []
        for elem in crack_tip_elements:
            for edge in elem.getElemEdges():
                edge_nodes = edge.getNodes()
                if crack_tip_node in edge_nodes:
                    nodes = [n.label for n in edge_nodes]
                    crack_tip_edge_nodes.append(sorted(nodes))
        
        unique_label_sets = list(set(tuple(lst) for lst in crack_tip_edge_nodes))
        unique_label_sets = [list(t) for t in unique_label_sets]

        for node_set in unique_label_sets:
            nodes = [self.part.nodes.getFromLabel(label) for label in node_set]
            
            # Identify the crack tip node in this set.
            A = None
            for n in nodes:
                if n.label == crack_tip_node.label:
                    A = n
                    break
            if A is None:
                # If for some reason the crack tip node isn't in the set, skip it.
                continue
            
            # The other two nodes in the set.
            candidates = [n for n in nodes if n.label != crack_tip_node.label]
            
            # Determine which candidate is farther from the crack tip (A).
            distances = [sum((n.coordinates[i] - A.coordinates[i])**2 for i in range(3)) for n in candidates]
            crack_tip_element_side_lengths.append(max(distances))

            far_node = candidates[distances.index(max(distances))]
            mid_node = candidates[distances.index(min(distances))]

            # Compute the new mid node coordinates as the quarter point between A and the far node.
            A_coords = A.coordinates
            B_coords = far_node.coordinates
            M_new = tuple(A_coords[i] + 0.25 * (B_coords[i] - A_coords[i]) for i in range(3))
                        
            # Update only the mid node's coordinates.
            self.part.editNode(nodes=(mid_node,),
                                coordinate1=M_new[0],
                                coordinate2=M_new[1],
                                coordinate3=M_new[2])
            
        return self
