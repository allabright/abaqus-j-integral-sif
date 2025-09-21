from abaqus import mdb
from abaqusConstants import *
import ast
import ConfigParser
import mesh
import os
import sys

class CreateValidationModel(object):

    def __init__(self):
        pass

    def parse_arguments(self):

        self.dimensions = sys.argv[-4]
        self.model_name = sys.argv[-3]
        self.run_dir = sys.argv[-2]
        self.model_path = sys.argv[-1]

        os.chdir(self.run_dir)
            
        return self
    
    def read_configuration_data(self):

        config_path = os.path.join(self.run_dir, "config.ini")
        config = ConfigParser.ConfigParser()
        config.read(config_path)

        # Read geometry data.
        self.part_height_mm = float(config.get("Geometry", "part_height_mm"))
        self.part_width_mm = float(config.get("Geometry", "part_width_mm"))
        self.part_thickness_mm = float(config.get("Geometry", "part_thickness_mm"))
        self.crack_length_mm = float(config.get("Geometry", "crack_length_mm"))
        self.crack_tip_coordinates = ast.literal_eval(config.get("Geometry", "crack_tip_coordinates"))

        if self.dimensions == "2d":
            self.crack_tip_coordinates = self.crack_tip_coordinates[:2]
            
        self.tolerance_mm = float(config.get("Geometry", "tolerance_mm"))

        # Read loading data.
        self.applied_stress_mpa = -abs(float(config.get("Loading", "applied_stress_mpa")))

        # Read material data.
        self.material_name = config.get("Material", "material_name")
        self.material_E = float(config.get("Material", "material_E"))
        self.material_nu = float(config.get("Material", "material_nu"))

        # Read modelling data.
        self.partitioning_method = config.get("Modelling", "partitioning_method")
        self.step_name = config.get("Modelling", "step_name")
        self.run_job = config.getboolean("Modelling", "run_job")
        self.use_quarter_point_elements = config.getboolean("Modelling", "use_quarter_point_elements")

        # Read partitioning data.
        self.horizontal_partition_count = int(config.get("Partitioning", "horizontal_partition_count"))
        self.horizontal_partition_bias = float(config.get("Partitioning", "horizontal_partition_bias"))

        if self.dimensions == "2d" and self.partitioning_method == "circular":
            self.circular_partition_count = int(config.get("Partitioning", "circular_partition_count"))
            self.spoke_partition_count = int(config.get("Partitioning", "spoke_partition_count"))
            self.circular_partition_bias = float(config.get("Partitioning", "circular_partition_bias"))
            self.circle_inner_radius_factor = float(config.get("Partitioning", "circle_inner_radius_factor"))
            self.circle_outer_radius_factor = float(config.get("Partitioning", "circle_outer_radius_factor"))
            
        # Read meshing data.
        self.crack_element_count = int(config.get("Meshing", "crack_element_count"))
        self.crack_element_bias = float(config.get("Meshing", "crack_element_bias"))
        self.coarse_seed_size_mm = float(config.get("Meshing", "coarse_seed_size_mm"))

        if self.dimensions == "2d" and self.partitioning_method == "circular":
            self.circle_arc_element_count = int(config.get("Meshing", "circle_arc_element_count"))
            self.circle_spoke_element_count = int(config.get("Meshing", "circle_spoke_element_count"))

        elif self.dimensions == "3d":
            self.through_thickness_element_count = int(config.get("Meshing", "through_thickness_element_count"))

        return self

    def print_abaqus(self, title=None, spacers=False, spacer_line = "-" * 120):
        if spacers == True:
            print >> sys.__stdout__, spacer_line
        if title is not None:
            print >> sys.__stdout__, title
            if spacers == True:
                print >> sys.__stdout__, spacer_line
        return self
    
    def set_z_and_transform(self):
        """
        Set the general z-dimension to be used for the model. For the 2D model this is 0, for the 3D model it is half of the part
        thickness, as the model is created in two halves. The transform variable is used to ensure that the crack tip stays on the
        origin, even in the z-direction.
        """

        if self.dimensions == "2d":
            self.z = 0.0
            self.transform = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
        elif self.dimensions == "3d":
            self.z = self.part_thickness_mm / 2.0
            self.transform = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, self.part_thickness_mm/2)

        return self

    def create_model(self):
        """
        Create the model in Abaqus using the name defined in the configuration parameters.
        """
        self.model = mdb.Model(name="Model_" + self.analysis_name)
        if 'Model-1' in mdb.models.keys():
            del mdb.models['Model-1']

        return self
    
    def create_section(self):
        """
        Create the section to be used for the model. This is a solid section for both the 2D and 3D models.
        """
        if self.dimensions == "2d":
            self.model.HomogeneousSolidSection(
                name='PlateSection',
                material=self.material_name,
                thickness=self.part_thickness_mm
            )
        elif self.dimensions == "3d":
            self.model.HomogeneousSolidSection(name='SolidSection', material=self.material_name)

        return self
    
    def assign_section(self):
        """
        Assign the defined section to the model. The 2D model assigns it to faces, and the 3D model assigns it to cells.
        """

        if self.dimensions == "2d":

            self.part.Set(name='AllFaces', faces=self.part.faces[:])
            self.part.SectionAssignment(
                region=self.part.sets['AllFaces'],
                sectionName='PlateSection',
                offset=0.0,
                offsetType=MIDDLE_SURFACE,
                offsetField='',
                thicknessAssignment=FROM_SECTION
            )

        elif self.dimensions == "3d":

            self.part.Set(cells=self.part.cells[:], name='WholeBlock')
            self.part.SectionAssignment(
                region=self.part.sets['WholeBlock'],
                sectionName='SolidSection')

        return self
    
    def create_material(self):
        """
        Define the material using the parameters defined in the configuration file.
        """

        self.model.Material(name=self.material_name)
        self.model.materials[self.material_name].Elastic(table=((self.material_E, self.material_nu),))

        return self
    
    def define_mesh_options(self):
        """
        Define the mesh options for the 2D and 3D elements.
        """

        if self.dimensions == "2d":
            elem_type_tri = mesh.ElemType(elemCode=CPS6, elemLibrary=STANDARD)
            self.part.setMeshControls(regions=self.part.faces[:], elemShape=TRI, technique=FREE)
            self.part.setElementType(regions=(self.part.faces[:],), elemTypes=(elem_type_tri,))

        elif self.dimensions == "3d":
            elem_type_tet = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
            self.part.setMeshControls(regions=self.part.cells[:], elemShape=TET, technique=FREE)
            self.part.setElementType(regions=(self.part.cells[:],), elemTypes=(elem_type_tet,))

        return self
    
    def mesh_part(self):

        if self.dimensions == "2d":
            self.part.generateMesh()
        elif self.dimensions == "3d":
            self.assembly.generateMesh(regions=((self.instance),))

        return self
    
    def create_load_step(self):

        self.model.StaticStep(name=self.step_name, previous="Initial")

        return self
    
    def create_assembly(self):

        self.assembly = self.model.rootAssembly
        self.assembly.DatumCsysByDefault(CARTESIAN)

        return self
    
    def create_instance(self):

        if self.dimensions == "2d":
            self.instance = self.assembly.Instance(name="Instance_" + self.analysis_name, part=self.part, dependent=ON)
        elif self.dimensions == "3d":
            self.instance = self.assembly.Instance(name="Instance_" + self.analysis_name, part=self.part, dependent=OFF)

        return self
    
    def set_field_outputs(self):

        if 'F-Output-1' in self.model.fieldOutputRequests.keys():
            foRequest = self.model.fieldOutputRequests['F-Output-1']
            foRequest.setValues(variables=('S', 'EE', 'U', 'RF', 'SENER'))
        else:
            foRequest = self.model.FieldOutputRequest(name='F-Output-1',
                                                createStepName=self.step_name,
                                                variables=('S', 'EE', 'U', 'RF'))
            
        return self
    
    def create_job(self):

        self.job = mdb.Job(
            name = self.analysis_name,
            model = "Model_" + self.analysis_name,
            description = 'Static Strength Analysis Job.'
            )
        
        return self
    
    def submit_job(self):

        if self.run_job == True:
            self.job.submit()
            self.job.waitForCompletion()

        return self
    
    def save_model(self):

        mdb.saveAs(pathName=self.model_path)

        return self
        