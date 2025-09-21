import logging
import numpy as np
from .shared_methods_mixin import *
logging.basicConfig(level=logging.INFO, format='%(message)s')

class DomainMethodsMixin(SharedMethodsMixin):

    @SharedMethodsMixin.timeit
    def define_integration_domains(self):
        """
        Defines multiple annular integration domains based on the distance from the crack tip.
        Each domain is stored in a separate dictionary, with elements allowed in multiple domains.
        """
        logging.info("Defining integration domains around the crack tip.")
        
        crack_tip = np.array(self.crack_tip_coordinates)
    
        # Create an array containing the outer radii of each annular domain.
        radii = np.linspace(self.domain_r_min, self.domain_r_max, self.number_of_domains + 1)
        
        # Build the domains, storing the inner and outer radii of each domain.
        self.domains = {}
        for i in range(self.number_of_domains):
            r_min = radii[i]
            r_max = radii[i + 1]
            self.domains[i + 1] = {
                "r_min": r_min,
                "r_max": r_max,
                "elements": []
            }

        # For each element, if any nodal coordinate lies in the domain, include that element in the domain.
        for elem_id, elem_data in self.elements.items():
            for node_coords in elem_data["nodal_coordinates_original"]:
                distance = np.linalg.norm(node_coords[:2] - crack_tip[:2])
                for i in range(self.number_of_domains):
                    r_min = radii[i]
                    r_max = radii[i + 1]
                    if r_min <= distance <= r_max:
                        self.domains[i + 1]["elements"].append(elem_id)

        # Remove duplicate elements to ensure each element only appears once in each domain.
        for domain in self.domains:
            self.domains[domain]["elements"] = list(set(self.domains[domain]["elements"]))

        return self
    
    @SharedMethodsMixin.timeit
    def filter_elements_and_nodes_by_integral_domain(self, extension_factor=2.0):
        """
        Filters elements that are within a specified distance from the crack tip.
        
        args:
            extension_factor (float): The extension of the maximum domain size, for which elements are captured.
                                      i.e. an extension factor of 2.0 means all nodes/elements within double the radius
                                      of the domain are captured.
        """
        logging.info("Filtering elements and nodes to keep only those within the maximum integral domain radius.")

        # Get the actual minimum and maximum domain radii, rather than relative to the crack tip.
        self.domain_r_min = self.domain_r_min_factor * self.crack_length_mm
        self.domain_r_max = self.domain_r_max_factor * self.crack_length_mm

        domain_elements, domain_nodes = [], []
        filtered_elements, filtered_nodes = {}, {}

        # Extend beyond the actual maximum integration domain by some factor to make sure that all nodes are captured.
        r_max = self.domain_r_max * extension_factor

        # If a node falls within the maximum domain, capture it, along with any elements it is connected to.
        for node_id, node_data in self.nodes.items():
            node_r = np.linalg.norm(np.array(node_data["coordinates_original"]) - self.crack_tip_coordinates)
            if  node_r <= r_max:
                domain_nodes.append(int(node_id))
                for element in node_data["connected_elements"]:
                    domain_elements.append(int(element))

        domain_elements = sorted(list(set(domain_elements)))
        domain_nodes = sorted(list(set(domain_nodes)))

        for element in domain_elements:
            filtered_elements[element] = self.elements[str(element)].copy()

        for node in domain_nodes:
            filtered_nodes[node] = self.nodes[str(node)].copy()

        self.elements = filtered_elements
        self.nodes = filtered_nodes

        return self
