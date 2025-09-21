import logging
import numpy as np
from .shared_methods_mixin import *
logging.basicConfig(level=logging.INFO, format='%(message)s')

class NodeMethodsMixin(SharedMethodsMixin):

    @SharedMethodsMixin.timeit
    def calculate_displaced_nodal_coordinates(self):
        """
        Calculate the displaced coordinates for each node. Abaqus only gives the original coordinates along with the displacement vectors.
        """
        logging.info("Calculating the displaced coordinates for each node.")

        for node in self.nodes:
            self.nodes[node]["coordinates_deformed"] = (np.array(self.nodes[node]["coordinates_original"]) + np.array(self.nodes[node]["displacement"])).tolist()

        return self
