from collections import defaultdict
import ast
import ConfigParser
import json
import os
import sys
from odbAccess import *


class ExportResults:
    """
    The purpose of this class is to export results from an Abaqus run. The end result is a JSON file with a key for each element.
    Each element then contains data on the nodes which connect to it, including the coordinates and displacements of those nodes.
    The stress data and strain energy data for each element is also stored.
    """
    def __init__(self,
                 elements = defaultdict(dict),
                 nodes = defaultdict(dict),
                 crack_tip_nodes = defaultdict(dict),
                 shear_stress_factor = 1.0, # Abaqus provides shear stresses already in tensorial form.
                 shear_strain_factor = 2.0 # Abaqus provides shear strains in engineering form. We need to divide by 2 to convert to tensorial form.
                 ):
        self.print_abaqus("Exporting Data from Abaqus", spacers=True)
        self.elements = elements
        self.nodes = nodes
        self.crack_tip_nodes = crack_tip_nodes
        self.shear_stress_factor = shear_stress_factor
        self.shear_strain_factor = shear_strain_factor

    def parse_arguments(self):

        self.print_abaqus("Parsing command line arguments")

        self.dimensions = sys.argv[-4]
        self.model_name = sys.argv[-3]
        self.run_dir = sys.argv[-2]
        self.odb_path = sys.argv[-1]

        os.chdir(self.run_dir)

        return self
    
    def read_configuration_data(self):

        config_path = os.path.join(self.run_dir, "config.ini")

        config = ConfigParser.ConfigParser()
        config.read(config_path)

        self.crack_tip_coordinates = ast.literal_eval(config.get("Geometry", "crack_tip_coordinates"))

        if self.dimensions == "2d":
            self.crack_tip_coordinates = self.crack_tip_coordinates[:2]
            
        self.tolerance_mm = float(config.get("Geometry", "tolerance_mm"))

        return self

    def load_model_parameters(self):

        self.print_abaqus("Loading model parameters")

        self.instance_name = "INSTANCE_" + self.model_name.upper()

        self.odb = openOdb(self.odb_path)
        self.assembly = self.odb.rootAssembly
        self.instance = self.odb.rootAssembly.instances[self.instance_name]
        self.frame = self.odb.steps["Load_Step"].frames[-1]

        return self

    def get_crack_tip_nodes(self):
        """
        The modelling script creates a set for the crack tip nodes. This is retrieved here.
        """

        self.print_abaqus("Retrieving crack tip nodes")

        self.crack_tip_nodes = []

        if self.dimensions == "2d":
            crack_seam = self.instance.nodeSets["CRACKSEAM_NODES"]
            for node in crack_seam.nodes:
                self.crack_tip_nodes.append(node.label)

        elif self.dimensions == "3d":
            crack_seam_nodes = self.assembly.nodeSets["CRACKSEAM_NODES"]
            for node_set in crack_seam_nodes.nodes:
                for node in node_set:
                    self.crack_tip_nodes.append(node.label)

        return self

    def get_element_stresses_and_strains(self):
        """
        This function gets the stress and strain energy field outputs from the Abaqus ODB file,
        and saves them to the element data dictionary.
        """

        self.print_abaqus("Getting element stress and strain data")

        fields = [
            ("S", "stress_tensor", self.shear_stress_factor),
            ("EE", "strain_tensor", self.shear_strain_factor)]
        
        for field_name, storage_key, factor in fields:

            field = self.frame.fieldOutputs[field_name]
            field_values = list(field.getSubset(position=INTEGRATION_POINT).values)

            for value in field_values:
                elem_id = int(value.elementLabel)
                if storage_key not in self.elements[elem_id]:
                    self.elements[elem_id][storage_key] = {}
                    
                if factor:
                    self.elements[elem_id][storage_key][int(value.integrationPoint)] = self._convert_stress_strain_data(value.data, factor)
                else:
                    self.elements[elem_id][storage_key][int(value.integrationPoint)] = value.data

        return self

    def get_node_coordinates_and_displacements(self):
        """
        This function gets the coordinates of each node, along with the displacements, from the Abaqus ODB file. They are saved to the node data dictionary.
        """

        self.print_abaqus("Getting node coordinate and displacement data")

        displacement_field = list(self.frame.fieldOutputs["U"].values)

        for node in self.instance.nodes:
            self.nodes[int(node.label)]["coordinates_original"] = self._convert_coordinate_displacement_data(node.coordinates)

        for displacement in displacement_field:
            self.nodes[int(displacement.nodeLabel)]["displacement"] = self._convert_coordinate_displacement_data(displacement.data)

        return self

    def get_element_and_node_connectivity(self):
        """
        This function creates a list for each element, containing the nodes which are connected to that element.
        """

        self.print_abaqus("Getting element and node connectivity")

        # N1, N2, and N3 are the corner nodes.
        # N4 is the node between N1 and N2, N5 is the node between N2 and N3, and N6 is the node between N1 and N3.
        for element in self.instance.elements:
            self.elements[int(element.label)]["connected_nodes"] = list(element.connectivity)

        for element in self.elements:
            for node in self.elements[element]["connected_nodes"]:
                if "connected_elements" not in self.nodes[node]:
                    self.nodes[node]["connected_elements"] = []
                self.nodes[node]["connected_elements"].append(element)

        for node in self.nodes:
            self.nodes[node]["connected_elements"] = sorted(self.nodes[node]["connected_elements"])

        return self
    
    def save_output_data_to_json(self):
        """
        Save the processed results data to a JSON file inside the run directory.
        """

        self.print_abaqus("Saving exported data to JSON file", spacers=True)

        self.output_data = {
            "crack_tip_nodes": self.crack_tip_nodes,
            "element_data": self.elements,
            "node_data": self.nodes
        }
        filename = "Export_" + self.model_name + ".json"

        output_data_path = os.path.join(self.run_dir, filename)

        with open(output_data_path, "w") as f:
            json.dump(self.output_data, f, indent=4)

        return self

    def print_abaqus(self, title=None, spacers=False, spacer_line = "-" * 120):

        if spacers == True:
            print >> sys.__stdout__, spacer_line

        if title is not None:
            print >> sys.__stdout__, title

            if spacers == True:
                print >> sys.__stdout__, spacer_line

        return self

    def _convert_stress_strain_data(self, data, factor):
        """
        Convert stresses and strains from the 1D array provided by Abaqus into a 2D array.
        Convert the strains from engineering to tensorial form.
        """

        if self.dimensions == "2d":
            return [
                [float(data[0]),          float(data[3] / factor)],
                [float(data[3] / factor), float(data[1])         ]]

        elif self.dimensions == "3d":
            return [
                [float(data[0]),          float(data[3] / factor), float(data[4] / factor)],
                [float(data[3] / factor), float(data[1]),          float(data[5] / factor)],
                [float(data[4] / factor), float(data[5] / factor), float(data[2])         ]]
        
        return 0

    def _convert_coordinate_displacement_data(self, data):
        """
        Get the coordinates and displacements relative to the crack tip. Usually the crack tip is at the origin in any case.
        """

        if self.dimensions == "2d":
            return [float(data[0] - self.crack_tip_coordinates[0]),
                    float(data[1] - self.crack_tip_coordinates[1])]

        elif self.dimensions == "3d":
            return [float(data[0] - self.crack_tip_coordinates[0]),
                    float(data[1] - self.crack_tip_coordinates[1]),
                    float(data[2] - self.crack_tip_coordinates[2])]
        
        return 0

if __name__ == "__main__":

    export = (
        ExportResults()
            .parse_arguments()
            .read_configuration_data()
            .load_model_parameters()
            .get_crack_tip_nodes()
            .get_element_stresses_and_strains()
            .get_node_coordinates_and_displacements()
            .get_element_and_node_connectivity()
            .save_output_data_to_json()
    )
