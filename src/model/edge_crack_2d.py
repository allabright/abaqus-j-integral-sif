import math
from abaqusConstants import *
from src.model.edge_crack_general import CreateValidationModelEdgeCrack


class CreateValidationModelEdgeCrack2D(CreateValidationModelEdgeCrack):

    def __init__(self,
                 analysis_name = "Edge_Crack_2D",
                 dimensions = "2d"
                 ):
        
        super(CreateValidationModelEdgeCrack2D, self).__init__()
        self.analysis_name = analysis_name
        self.dimensions = dimensions        
              
    def calculate_circle_radii(self):
        """
        Calculate the radii to use when using the circular partitioning scheme.
        """
        self.circle_radii = []
        self.circle_inner_radius = self.crack_length_mm * self.circle_inner_radius_factor
        self.circle_outer_radius = self.crack_length_mm * self.circle_outer_radius_factor

        if self.number_of_circles <= 1:
            self.circle_radii.append(self.circle_inner_radius)

        else:    
            for i in range(self.number_of_circles):
                s = float(i) / (self.number_of_circles - 1)
                sbias = s ** self.circle_bias
                # Interpolate between inner and outer
                r = self.circle_inner_radius + (self.circle_outer_radius - self.circle_inner_radius) * sbias

                self.circle_radii.append(r)
        
        return self

    def partition_circle_around_crack(self):
        """
        Partition the set of circles around the crack tip, when using the circular partitioning scheme.
        """

        circleSketch = self.model.ConstrainedSketch(name='CrackTipCircleSketch', sheetSize=500.0, transform=self.transform)

        for radius in self.circle_radii:
            
            circleSketch.CircleByCenterPerimeter(
                center=(0.0, 0.0),
                point1=(radius, 0)
            )

        top_face = self.part.faces.getClosest(coordinates=((0.0, self.tolerance_mm, self.part_thickness_mm / 2.0), ))
        bottom_face = self.part.faces.getClosest(coordinates=((0.0, -self.tolerance_mm, self.part_thickness_mm / 2.0), ))

        partition_faces = (top_face[0][0], bottom_face[0][0])

        self.part.PartitionFaceBySketch(faces=partition_faces, sketch=circleSketch)

        del self.model.sketches['CrackTipCircleSketch']

        return self
       
    def partition_radially_around_crack(self):
        """
        Partition the radial spokes within the circles around the crack tip, when using the circular partitioning scheme.
        """

        spokeSketch = self.model.ConstrainedSketch(name='SpokeSketch', sheetSize=500.0, transform=self.transform)       

        for i in range(self.number_of_spokes):
            angle = 2.0 * math.pi * (float(i) / self.number_of_spokes)
            x_end = self.circle_outer_radius * math.cos(angle)
            y_end = self.circle_outer_radius * math.sin(angle)

            if y_end != 0:
                spokeSketch.Line(point1=(0.0, 0.0), point2=(x_end, y_end))
                spokeSketch.Line(point1=(0.0, 0.0), point2=(x_end, -y_end))

        top_face = self.part.faces.getClosest(coordinates=((0.0, self.tolerance_mm, self.part_thickness_mm / 2.0), ))
        bottom_face = self.part.faces.getClosest(coordinates=((0.0, -self.tolerance_mm, self.part_thickness_mm / 2.0), ))

        partition_faces = (top_face[0][0], bottom_face[0][0])

        self.part.PartitionFaceBySketch(faces=partition_faces, sketch=spokeSketch)
        del self.model.sketches['SpokeSketch']

        return self
    
    def seed_arc_segments(self):
        """
        Seed the arc segments (the edges along each circle between two adjacent spokes).
        """

        self.arc_edges = []
       
        for r in self.circle_radii:
            # Consider each pair of consecutive spokes
            for j in range(self.number_of_spokes):
                # Define an arc section swept along the circle
                angle_start = 2.0 * math.pi * j / self.number_of_spokes
                angle_end   = 2.0 * math.pi * (j + 1) / self.number_of_spokes
                mid_angle = 0.5 * (angle_start + angle_end)
                
                # Coordinates for midpoint on this arc
                x_mid = r * math.cos(mid_angle)
                y_mid = r * math.sin(mid_angle)
                
                try:
                    found_edges = self.part.edges.findAt(((x_mid, y_mid, self.z),))       
                    for edge in found_edges:
                        self.arc_edges.append(edge)
                    if self.dimensions == "3d":
                        found_edges_2 = self.part.edges.findAt(((x_mid, y_mid, -self.z),))       
                        for edge in found_edges_2:
                            self.arc_edges.append(edge)
                except Exception as e:
                    print("Could not find arc edge at r={r:.3f}, angles=({angle_start:.3f},{angle_end:.3f}): {e}")

        # Seed each arc edge with a fixed number of elements
        if self.arc_edges:
            self.part.seedEdgeByNumber(edges=self.arc_edges, number=self.circle_arc_elements)

        return self
    
    def seed_radial_segments(self):
        """
        Seed the radial segments (spokes) created when using the circular partitioning scheme.
        """

        self.radial_edges = []

        for spoke_index in range(self.number_of_spokes):
            angle = 2.0 * math.pi * spoke_index / self.number_of_spokes
            cosA = math.cos(angle)
            sinA = math.sin(angle)
            
            for i in range(len(self.circle_radii)):

                if i == 0:
                    r_i = 0.0
                else:
                    r_i = self.circle_radii[i-1]
                r_j = self.circle_radii[i]
                
                r_mid = 0.5 * (r_i + r_j)  
                # Cartesian coordinates of that midpoint
                x_mid = r_mid * cosA
                y_mid = r_mid * sinA
                
                try:
                    found_edges = self.part.edges.findAt(((x_mid, y_mid, self.z),))
                    for edge in found_edges:
                        self.radial_edges.append(edge)

                    if self.dimensions == "3d":
                        found_edges_2 = self.part.edges.findAt(((x_mid, y_mid, -self.z),))
                        for edge in found_edges_2:
                            self.radial_edges.append(edge)
                except Exception as e:
                    print("Could not find edge for spoke {spoke_index}, segment {i}: {e}")

        # Seed each segment edge with a fixed number of elements
        if self.radial_edges:
            self.part.seedEdgeByNumber(edges=self.radial_edges, number=self.circle_radius_elements)
        
        return self
    
    def seed_rest_of_part_circular_partitioning(self):
        """
        When using the circular partitioning scheme, seed the rest of the part, excluding the arc segments and
        radial segments which have already been seeded.
        """

        allEdges = self.part.edges[:]

        if self.radial_edges and self.arc_edges:
            otherEdges = [ed for ed in allEdges if ed not in self.radial_edges and ed not in self.arc_edges]

        if otherEdges:
            self.part.seedEdgeBySize(edges=otherEdges, size=self.coarse_seed_size_mm, deviationFactor=0.1)

        return self
    
if __name__ == '__main__':

    run = (
        CreateValidationModelEdgeCrack2D()
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
            
            .create_section()
            .assign_section()
            .define_mesh_options()    
    )
            

    if run.partitioning_method == "regular":
        (run
            .partition_edges()
            .seed_crack_tip_region()
            .seed_rest_of_part_regular_partitioning()
        )

    elif run.partitioning_method == "circular":
        (run
            .calculate_circle_radii()
            .partition_radially_around_crack()
            .partition_circle_around_crack()
            .seed_radial_segments()
            .seed_arc_segments()
            .seed_rest_of_part_circular_partitioning()
        )

    run.mesh_part()

    if run.use_quarter_point_elements:
        run.move_crack_tip_nodes()

    (run
        .create_crack_seam_set_geometry()
        .create_crack_seam_set_nodes()
        .create_sets_for_top_and_bottom_edges()
        .create_assembly()
        .create_instance()
        .apply_top_pressure_load()
        .apply_bottom_fixed_boundary_condition()
        .create_job()
        .submit_job()
        .save_model()
     )
