import ast
import configparser
import logging
import os
import time
import numpy as np
import pandas as pd
from functools import wraps
logging.basicConfig(level=logging.INFO, format='%(message)s')

class SharedMethodsMixin:

    @staticmethod
    def timeit(func):
        """
        Function to give the time taken to perform each step of the calculation.
        """
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.perf_counter()
            result = func(*args, **kwargs)
            elapsed_time = time.perf_counter() - start_time
            print(f"COMPLETED: {elapsed_time:.3f} seconds")
            print()
            return result
        return wrapper

    def print_spacers(self, title=None, spacers=False, spacer_line = "-" * 120):
            """
            Pretty print using spacers.
            """
            if spacers:
                print(spacer_line)
            if title:
                print(title)
                if spacers:
                    print(spacer_line)
            return 0
    
    def read_configuration_data(self):
        """
        Read in the configuration data from the config.ini file.
        """

        config_path = os.path.join(self.run_dir, "config.ini")

        config = configparser.ConfigParser()
        config.read(config_path)

        # Read analysis data.
        self.configuration = config.get("Analysis", "configuration")
        self.stress_state = config.get("Analysis", "stress_state")
        self.weight_function_type = config.get("Analysis", "weight_function_type")
        self.number_of_domains = int(config.get("Analysis", "number_of_domains"))
        self.domain_r_min_factor = float(config.get("Analysis", "domain_r_min_factor"))
        self.domain_r_max_factor = float(config.get("Analysis", "domain_r_max_factor"))

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
    
    @timeit
    def calculate_analytical_values(self):
        """
        Calculate an analytical value for the stress intensity factor using the formula obtained from
        'Research on the stress intensity factor of crack on 'the finite width plate with edge damage based on the finite element method'
        (Hong Yuan & Hao Zhang, 2023).
        Use the stress intensity factor to calculate a J-integral for the configuration.
        """
        logging.info("Calculating analytical values for the SIF and the J-integral.")

        a = self.crack_length_mm
        b = self.part_width_mm

        # Shape factor formula Beta for an edge-crack in a finite plate.
        F = 1.12 - (0.23 * a / b) + (10.6 * (a / b) ** 2) - (21.7 * (a / b) ** 3) + (30.4 * (a / b) ** 4)

        self.sif_analytical = F * np.abs(self.applied_stress_mpa) * np.sqrt(np.pi * a)
        self.j_analytical = (self.sif_analytical ** 2) / self.E_prime

        return self
    
    @timeit
    def calculate_minimum_element_size(self):
        """
        Calculate the minimum element size near the crack tip. Required for the mesh sensitivity study.
        """
        # Uniform mesh case: bias of 1 means no growth.
        if self.crack_element_bias == 1:
            return self.crack_length_mm / float(self.crack_element_count)
        
        # Calculate growth factor r from the bias ratio and number of elements.
        r = self.crack_element_bias ** (1.0 / (self.crack_element_count - 1))
        
        # Calculate the minimum element size using the geometric series sum formula.
        self.min_element_size_mm = np.round(self.crack_length_mm * (r - 1) / (r**self.crack_element_count - 1), 3)

        return self
    
    def convert_results_to_dataframe(self):
        """
        Convert the results of the tool to a pandas DataFrame, so that charts can be more easily created.
        """        
        rows = []
        for domain_id_str, domain_data in self.output_data_dict["results"].items():
            
            row = {
                "configuration": self.output_data_dict["configuration"],
                "stress_state": self.output_data_dict["stress_state"],
                "weight_function_type": self.output_data_dict["weight_function_type"],
                "crack_length_mm": self.output_data_dict["crack_length_mm"],
                "minimum_element_size_mm": self.min_element_size_mm,
                "domain_id": int(domain_id_str),
                "domain_r_min": domain_data["domain_r_min"],
                "domain_r_max": domain_data["domain_r_max"],
                "j_integral_value_analytical": domain_data["j_integral_value_analytical"],
                "j_integral_value_calculated": domain_data["j_integral_value_calculated"],
                "j_integral_error_percentage": domain_data["j_integral_error_percentage"],
                "stress_intensity_factor_value_analytical": domain_data["stress_intensity_factor_value_analytical"],
                "stress_intensity_factor_value_calculated": domain_data["stress_intensity_factor_value_calculated"],
                "stress_intensity_factor_error_percentage": domain_data["stress_intensity_factor_error_percentage"],
            }
            rows.append(row)
        
        self.output_data_df = pd.DataFrame(rows)
        self.output_data_df.sort_values(by="domain_id", inplace=True)
        self.output_data_df.reset_index(drop=True, inplace=True)
        
        return self
