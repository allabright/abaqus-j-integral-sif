import logging
import numpy as np
from .shared_methods_mixin import *
logging.basicConfig(level=logging.INFO, format='%(message)s')

class ElementMethodsMixin(SharedMethodsMixin):

    @SharedMethodsMixin.timeit
    def get_element_nodal_coordinates(self):
        """
        Get a set of all the nodal coordinates for each element.
        """
        logging.info("Retrieving the nodal coordinates for each element.")

        failed_elements = []

        for element_id, element_data in self.elements.items():

            # Store both the original and displaced nodal coordinates.
            nodal_coords = []
            nodal_coords_deformed = []

            for node in element_data["connected_nodes"]:
                try:
                    nodal_coords.append(self.nodes[node]["coordinates_original"])
                    nodal_coords_deformed.append(self.nodes[node]["coordinates_deformed"])
                except:
                    failed_elements.append(element_id)

            element_data["nodal_coordinates_original"] = nodal_coords
            element_data["nodal_coordinates_deformed"] = nodal_coords_deformed

        for element in failed_elements:
            self.elements.pop(element, None)

        return self
    
    @SharedMethodsMixin.timeit
    def convert_element_stress_and_strain_tensors(self):
        """
        The stress and strain tensors in the export JSON file are stored as dictionaries due to the format.
        However, for calculations, it is more convenient to have them as arrays so they can be used with NumPy.
        """

        logging.info("Converting the element stress and strain tensors from dictionaries to arrays.")
        
        # Number of integration points per element.
        n_ip = len(self.integration_point_weights)

        for elem_id, elem_data in self.elements.items():

            for param in ["stress_tensor", "strain_tensor"]:
                tensor = []

                for ip in range(n_ip):
                    tensor.append(elem_data[param][str(ip+1)])

                elem_data[param] = tensor

        return self

    @SharedMethodsMixin.timeit
    def calculate_element_integration_point_coordinates(self):
        """
        Get the coordinates of the integration points for each element, as an array.
        """
        logging.info("Calculating the integration point coordinates for each element.")

        for element in self.elements:
            nodal_coords = np.array(self.elements[element][f"nodal_coordinates_{self.configuration}"])
            self.elements[element][f"integration_point_coordinates"] = np.dot(self.shape_functions_num, nodal_coords).tolist()

        return self

    @SharedMethodsMixin.timeit
    def calculate_element_strain_energy_densities(self):
        """
        Computes the strain energy density W for each element at each integration point.
        """
        logging.info("Calculating the strain energy density for each element.")

        for elem_id, elem_data in self.elements.items():

            # Get the lists of stress and strain tensors, one per integration point.
            stress_list = elem_data["stress_tensor"]
            strain_list = elem_data["strain_tensor"]
            n_ip = len(stress_list)
            W_list = []

            for ip in range(n_ip):
                stress_arr = np.array(stress_list[ip])
                strain_arr = np.array(strain_list[ip])

                # There is the option available here to convert from engineering to tensorial strain.
                # However, this is done during the data export.
                strain_factor = 1.0

                if self.dimensions == "2d":
                    strain_arr_factored = np.array([
                        [strain_arr[0][0]                , strain_arr[0][1]/2 / strain_factor],
                        [strain_arr[1][0] / strain_factor, strain_arr[1][1]                  ]
                    ])
                elif self.dimensions == "3d":
                    strain_arr_factored = np.array([
                        [strain_arr[0][0]                , strain_arr[0][1] / strain_factor, strain_arr[0][2] / strain_factor],
                        [strain_arr[1][0] / strain_factor, strain_arr[1][1]                , strain_arr[1][2] / strain_factor],
                        [strain_arr[2][0] / strain_factor, strain_arr[2][1] / strain_factor, strain_arr[2][2]                ]
                    ])
                # Compute the strain energy density at this integration point.
                W = 0.5 * (np.tensordot(stress_arr, strain_arr_factored))
                W_list.append(W)
                
            self.elements[elem_id]["strain_energy_density"] = W_list

        return self

    @SharedMethodsMixin.timeit
    def calculate_integration_point_jacobians(self):
        """
        Compute the Jacobian matrix and determinant for each integration point of each element
        using the shape function derivatives.
        """
        logging.info("Calculating the integration point Jacobians for each element.")

        for element_id, element_data in self.elements.items():
            nodal_coords = np.array(element_data[f"nodal_coordinates_{self.configuration}"])
            jacobians = np.tensordot(self.shape_function_derivatives_num, nodal_coords, axes=([1], [0]))
            det_jacobians = np.linalg.det(jacobians)

            element_data[f"jacobian"] = jacobians.tolist()
            element_data[f"jacobian_determinant"] = det_jacobians.tolist()

        return self
    
    @SharedMethodsMixin.timeit
    def calculate_element_displacement_gradients(self):
        """
        Calculate the displacement gradients for each integration point of each element, using the
        shape function derivatives, Jacobians, and nodal displacements.
        """
        logging.info("Calculating the displacement gradient for each element.")

        for elem_id, elem_data in self.elements.items():
            n_ip = len(self.integration_point_weights)
            gradients_list = []
            displacements = np.array([self.nodes[n][f"displacement"] for n in elem_data["connected_nodes"]])
            jacobians = np.array(elem_data[f"jacobian"])
            
            for ip in range(n_ip):
                invJ = np.linalg.inv(jacobians[ip])
                global_derivs = (invJ @ np.array(self.shape_function_derivatives_num)[ip].T).T
                gradients_list.append((displacements.T @ global_derivs).tolist())
            
            elem_data[f"displacement_gradients"] = gradients_list
        
        return self

    @SharedMethodsMixin.timeit
    def calculate_element_weights(self):
        """
        Calculate the element weight function value and derivative for each integration point of each element,
        for each annular domain.
        """
        logging.info("Calculating weights for each element based on distance to the crack tip.")
        crack_tip_coordinates = np.array(self.crack_tip_coordinates)
        dim = int(self.dimensions[0])
        
        for elem_id, elem_data in self.elements.items():
            ip_coords = elem_data.get(f"integration_point_coordinates", [])
            weight_vals = {domain: [] for domain in self.domains}
            weight_grads = {domain: [] for domain in self.domains}
            
            for ip_coord in ip_coords:
                ip_coord = np.array(ip_coord)
                # Calculate the distance from the crack tip to the integration point, in the x-y plane only.
                r_vec = ip_coord[:2] - crack_tip_coordinates[:2]
                r = np.linalg.norm(r_vec)
                
                for domain in self.domains:
                    r_min = self.domains[domain]["r_min"]
                    r_max = self.domains[domain]["r_max"]

                    # If the integration point is inside or outside the domain, gradient is zero, so no contribution.
                    if r <= r_min:
                        q = 1.0
                        grad_q = np.zeros(dim)
                    elif r >= r_max:
                        q = 0.0
                        grad_q = np.zeros(dim)
                    else:
                        # Normalize the distance through the domain. s=0 is at the inside edge of the domain, s=1 is at
                        # the outside edge of the domain.
                        s = (r - r_min) / (r_max - r_min)
                        if self.weight_function_type == "linear":
                            q = 1 - s
                            dq_dr = -1.0 / (r_max - r_min)
                        elif self.weight_function_type == "polynomial":
                            q = 1 - 3 * s**2 + 2 * s**3
                            dq_dr = (-6 * s + 6 * s**2) / (r_max - r_min)
                        elif self.weight_function_type == "gaussian":
                            sigma = 0.4
                            q = np.exp(-(s / sigma) ** 2)
                            dq_dr = (-2 * s / (sigma ** 2) * np.exp(-(s / sigma) ** 2)) / (r_max - r_min)
                        else:
                            raise ValueError(f"Unknown weight_function_type: {self.weight_function_type}")

                        # Get the gradient of the weight function at the integration point.
                        grad_q = dq_dr * (r_vec / r)

                    weight_vals[domain].append(q)
                    weight_grads[domain].append(grad_q.tolist())

            elem_data[f"weight_function_values"] = weight_vals
            elem_data[f"weight_function_gradients"] = weight_grads

        return self

    @SharedMethodsMixin.timeit
    def calculate_domain_j_integrals(self):
        """
        Calculate the J-integral for each domain.
        """
        logging.info("Calculating J-integral for each integral domain.")

        for domain_id, domain_data in self.domains.items():

            domain_elements = domain_data["elements"]

            # Total J-integral for the domain.
            j_value = 0

            for elem_id in domain_elements:
                elem_data = self.elements[elem_id]
                stress_tensor = np.array(elem_data["stress_tensor"])
                displacement_grads = np.array(elem_data[f"displacement_gradients"])
                strain_energy_density = np.array(elem_data["strain_energy_density"])
                weight_grads = elem_data[f"weight_function_gradients"][domain_id]
                n_ip = len(elem_data[f"integration_point_coordinates"])

                for ip in range(n_ip):

                    # Given as 2x2 arrays. xx = [0, 0], xy = [0, 1], yx = [1, 0], yy = [1, 1]
                    sigma = stress_tensor[ip]               
                    grad_u = displacement_grads[ip]    
                    W = strain_energy_density[ip]
                    dq_dx, dq_dy = weight_grads[ip][0], weight_grads[ip][1]
                    detJ_ip = elem_data[f"jacobian_determinant"][ip]
                    w_ip = self.integration_point_weights[ip]

                    term_x = sigma[0, 0] * grad_u[0, 0] + sigma[0, 1] * grad_u[1, 0]
                    term_y = sigma[1, 0] * grad_u[0, 0] + sigma[1, 1] * grad_u[1, 0]
                    
                    # Integrand for the individual integration point.
                    integrand = ((term_x - W) * dq_dx) + (term_y * dq_dy)
                    
                    # Multiply by detJ and integration point weight to account for the relative area/volume.
                    j_value += integrand * detJ_ip * w_ip

            domain_data[f"j_integral"] = j_value

        return self

    @SharedMethodsMixin.timeit
    def calculate_domain_stress_intensity_factors(self):
        """
        Converts the J-integral value into the Stress Intensity Factor (SIF).
        """
        logging.info("Calculating the stress intensity factor for each integral domain.")
        
        # Compute effective modulus E'
        if self.stress_state == "plane_stress":
            self.E_prime = self.material_E
        elif self.stress_state == "plane_strain":
            self.E_prime = self.material_E / (1 - self.material_nu**2)

        for domain in self.domains:
            self.domains[domain][f"stress_intensity_factor"] = np.sqrt((self.domains[domain][f"j_integral"]) * self.E_prime)

        return self
