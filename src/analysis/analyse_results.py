import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

import numpy as np

from .mixins.domain_methods_mixin import DomainMethodsMixin
from .mixins.element_methods_mixin import ElementMethodsMixin
from .mixins.node_methods_mixin import NodeMethodsMixin
from .mixins.shape_function_methods_mixin import ShapeFunctionMethodsMixin

class AnalyseResults(DomainMethodsMixin,
                     ElementMethodsMixin,
                     NodeMethodsMixin,
                     ShapeFunctionMethodsMixin):

    def __init__(self,
                 dimensions,
                 input_data,
                 run_dir
): 
        self.dimensions = dimensions
        self.input_data = input_data
        self.run_dir = run_dir

        if dimensions == "2d":
            self.integration_point_count = 3
            self.integration_point_weights = [1/6, 1/6, 1/6]
        
        elif dimensions == "3d":
            self.integration_point_count = 4
            self.integration_point_weights = [1/24, 1/24, 1/24, 1/24]
        

    def analyse_results(self):

        self.elements = self.input_data["element_data"].copy()
        self.nodes = self.input_data["node_data"].copy()

        logging.info("Collecting Parameters")
        (self
             .read_configuration_data()
             .filter_elements_and_nodes_by_integral_domain()
             .calculate_minimum_element_size()
             .calculate_displaced_nodal_coordinates()
             .get_element_nodal_coordinates()
             .convert_element_stress_and_strain_tensors()
             .define_integration_domains())
        
        logging.info("Computing Shape Functions")
        self.compute_shape_functions_and_derivatives()

        logging.info("Calculating Element Values")
        (self.calculate_element_integration_point_coordinates()
             .calculate_integration_point_jacobians()
             .calculate_element_strain_energy_densities()
             .calculate_element_displacement_gradients())
        
        logging.info("Performing Calculations")
        (self
             .calculate_element_weights()
             .calculate_domain_j_integrals()
             .calculate_domain_stress_intensity_factors()
             .calculate_analytical_values())
        
        self.output_data_dict = {
            "configuration": self.configuration,
            "stress_state": self.stress_state,
            "approx_min_element_size_mm": self.min_element_size_mm,
            "weight_function_type": self.weight_function_type,
            "applied_stress_mpa": np.abs(self.applied_stress_mpa),
            "crack_length_mm": self.crack_length_mm,
            "number_of_domains": self.number_of_domains,
            "minimum_domain_mm": np.round(self.domain_r_min, 4),
            "maximum_domain_mm": np.round(self.domain_r_max, 4),
            "results": {}
        }

        for domain_id, domain_data in self.domains.items():
            self.output_data_dict["results"][domain_id] = {}
            self.output_data_dict["results"][domain_id]["domain_r_min"] = domain_data["r_min"]
            self.output_data_dict["results"][domain_id]["domain_r_max"] = domain_data["r_max"]

            for param, unit, solution in ["j_integral", "", self.j_analytical], ["stress_intensity_factor", "MPa mm^1/2", self.sif_analytical]:
                val = domain_data[f'{param}']
                error = abs(((val - solution) / solution) * 100)
                self.output_data_dict["results"][domain_id][f"{param}_value_analytical"] = np.round(solution, 4)
                self.output_data_dict["results"][domain_id][f"{param}_value_calculated"] = np.round(val, 4)
                self.output_data_dict["results"][domain_id][f"{param}_error_percentage"] = np.round(error, 4)

        print("-" * 100)

        self.convert_results_to_dataframe()

        print("Analysis complete.")
        return self