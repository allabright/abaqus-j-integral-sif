from abaqusConstants import *
from src.model.edge_crack_general import CreateValidationModelEdgeCrack


class CreateValidationModelEdgeCrack3D(CreateValidationModelEdgeCrack):

    def __init__(self,
                 analysis_name = "Edge_Crack_3D",
                 dimensions = "3d"):
        
        super(CreateValidationModelEdgeCrack3D, self).__init__()
        self.analysis_name = analysis_name
        self.dimensions = dimensions
   
    def mirror_geometry(self):
        """
        Mirror the geometry to create two cells. This is necessary to keep the crack tip at the origin when using a 3D model,
        as Abaqus does not allow for a mid-point extrusion.
        """

        checkPoint = (0.0, 0.0, 0.0)
        faceList = self.part.faces.findAt((checkPoint,))
        if not faceList:
            raise ValueError("No face found at point {}. Check geometry.".format(checkPoint))
        main_face = faceList[0]
        self.part.Mirror(mirrorPlane=main_face, keepOriginal=True)       

        return self
        
    def partition_model(self, mirror=False):
        """
        Finds radial and arc edges, filters out boundary edges, and sweeps them downward through the full thickness
        to create a pie-slice-shaped partitioned structure. This is performed on the full model, after mirroring.
        """

        x_min, x_max = -self.crack_length_mm - self.tolerance_mm, self.part_width_mm - self.crack_length_mm + self.tolerance_mm
        y_min, y_max = -self.part_height_mm / 2.0 - self.tolerance_mm, self.part_height_mm / 2.0 + self.tolerance_mm
        
        if mirror:
            z_min = -self.z - self.tolerance_mm
            z_max = 0.0
            z_gen = -self.z
        else:
            z_min = 0.0
            z_max = self.z + self.tolerance_mm
            z_gen = self.z

        # Select candidate edges
        candidate_edges = self.part.edges.getByBoundingBox(
            xMin=x_min, xMax=x_max,
            yMin=y_min, yMax=y_max,
            zMin=z_min, zMax=z_max
        )

        cellMidPoint = (0.0, 0.0, z_gen)
        cell = self.part.cells.findAt((cellMidPoint,))
        
        internal_edges = []
        
        for edge in candidate_edges:
            x, y, z = edge.pointOn[0]
            if x != -self.crack_length_mm and x != self.part_width_mm - self.crack_length_mm and \
               y != self.part_height_mm / 2.0 and y != -self.part_height_mm / 2.0 and \
               z == z_gen:
                internal_edges.append(edge)
        
        direction_edge = self.part.edges.findAt(
            ((-self.crack_length_mm, self.part_height_mm / 2.0, z_gen / 2.0),)
        )

        self.part.PartitionCellBySweepEdge(
            cells=cell,
            edges=tuple(internal_edges),
            sweepPath=direction_edge[0]
        )

        return self

    def seed_through_thickness_edges(self):
        """
        Seed the through-thickness edges near the crack tip, in the z-direction.
        """

        self.through_thickness_edges = []

        for edge in self.part.edges[:]:

            # Get the vertices of each edge of the part, along with their starts and ends.
            vertex_1 = edge.getVertices()[0]
            vertex_2 = edge.getVertices()[1]

            vertex_1_coords_x = self.part.vertices[vertex_1].pointOn[0][0]
            vertex_1_coords_y = self.part.vertices[vertex_1].pointOn[0][1]

            vertex_2_coords_x = self.part.vertices[vertex_2].pointOn[0][0]
            vertex_2_coords_y = self.part.vertices[vertex_2].pointOn[0][1]

            # If the edge is within the central partition of the model, and is directly in the z-direction,
            # register it as a through-thickness edge to seed.
            if abs(vertex_1_coords_y) <= abs(self.partition_heights[1]):
                if abs(vertex_1_coords_x) <= self.crack_length_mm:
                    if vertex_1_coords_x == vertex_2_coords_x and vertex_1_coords_y == vertex_2_coords_y:
                        self.through_thickness_edges.append(edge)

        self.part.seedEdgeByNumber(self.through_thickness_edges, self.through_thickness_element_count)

        return self
        
if __name__ == '__main__':

    run = CreateValidationModelEdgeCrack3D()

    (run        
        .parse_arguments()
        .read_configuration_data()
        .calculate_horizontal_partition_heights()
        .set_z_and_transform()
        .create_model()
        .create_load_step()
        .set_field_outputs()
        .define_part_geometry()
        .create_part()
        .create_material()
        .partition_part_top_to_bottom()
        .mirror_geometry()
        .partition_model(mirror=False)
        .partition_model(mirror=True)
        .partition_edges()
        .create_section()
        .assign_section() 
        .define_mesh_options()     
        .seed_crack_tip_region()
        .seed_through_thickness_edges()
        .seed_rest_of_part_regular_partitioning()
        .create_sets_for_top_and_bottom_edges()
        .create_assembly()
        .create_instance()
        .create_crack_seam_set_geometry()
        .mesh_part()
        .create_crack_seam_set_nodes()
        .apply_top_pressure_load()
        .apply_bottom_fixed_boundary_condition()
        .create_job()
        .submit_job()
        .save_model()
    )