import logging
import numpy as np
import sympy as sp
from .shared_methods_mixin import *
logging.basicConfig(level=logging.INFO, format='%(message)s')

class ShapeFunctionMethodsMixin(SharedMethodsMixin):

    @SharedMethodsMixin.timeit
    def compute_shape_functions_symbolically_cps6(self):
        """
        Generate the symbolic shape functions for a CPS6 (6-node triangular) element,
        as defined in the Abaqus theory manual, section 3.2.6.
        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/stm/default.htm?startat=ch03s02ath64.html
        """
        logging.info("Computing symbolic shape functions for triangular CPS6 elements.")

        x, y = sp.symbols('x y')
        L1 = x
        L2 = y
        L3 = (1 - x - y)

        N1 = -(L3 * (1 - (2 * L3)))
        N2 = -(L1 * (1 - (2 * L1)))
        N3 = -(L2 * (1 - (2 * L2)))
        N4 = 4 * L1 * L3
        N5 = 4 * L1 * L2
        N6 = 4 * L2 * L3

        self.shape_functions_sym = [N1, N2, N3, N4, N5, N6]

        return self

    @SharedMethodsMixin.timeit
    def compute_shape_functions_symbolically_c3d10(self):

        """
        Generate the symbolic shape functions for a C3D10 (10-node tetrahedral) element,
        as defined in the Abaqus theory manual, section 3.2.6.
        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/stm/default.htm?startat=ch03s02ath64.html
        """
        logging.info("Computing symbolic shape functions for tetrahedral C3D10 elements.")

        x, y, z = sp.symbols('x y z')

        w1 = y * ((2 * y) - 1)
        w2 = z * ((2 * z) - 1)
        w3 = (1 - x - y - z) * (1 - (2 * x) - (2 * y) - (2 * z))
        w4 = x * ((2 * x) - 1)
        w5 = 4 * y * z
        w6 = (4 * z) * (1 - x - y - z)
        w7 = (4 * y) * (1 - x - y - z)
        w8 = 4 * x * y
        w9 = 4 * x * z
        w10 = (4 * x) * (1 - x - y - z)

        self.shape_functions_sym = [w1, w2, w3, w4, w5, w6, w7, w8, w9, w10]

        return self
    
    @SharedMethodsMixin.timeit
    def compute_shape_function_derivatives_symbolically(self):
        """
        Compute the symbolic derivatives of the shape functions with respect to the natural coordinates.
        For each shape function w, this function computes the partial derivatives dw/dx and dw/dy (and dw/dz for 3D elements).
        """
        logging.info("Computing symbolic shape function derivatives.")

        x, y, z = sp.symbols('x y z')
        derivative_vars = (x, y) if self.dimensions == "2d" else (x, y, z)

        self.shape_function_derivatives_sym = []

        for w in self.shape_functions_sym:
            w_derivs = []
            for var in derivative_vars:
                w_derivs.append(sp.diff(w, var))
            self.shape_function_derivatives_sym.append(w_derivs)

        return self

    @SharedMethodsMixin.timeit
    def compute_shape_functions_and_derivatives_numerically(self):
        """
        Numerically evaluate the symbolic shape functions and their derivatives for an element
        at given natural coordinate points.
        """
        logging.info("Performing numerical evaluation of shape functions and shape function derivatives.")

        if self.dimensions == "2d":
            get_sub_dict = lambda point: {'x': point[0], 'y': point[1]}
            dim = 2
        elif self.dimensions == "3d":
            get_sub_dict = lambda point: {'x': point[0], 'y': point[1], 'z': point[2]}
            dim = 3
        else:
            raise ValueError("Unsupported dimensions: {}".format(self.dimensions))
        
        self.shape_functions_num = []
        self.shape_function_derivatives_num = []

        for point in self.natural_points:
            sub_dict = get_sub_dict(point)

            shape_function_single = [float(function.subs(sub_dict)) for function in self.shape_functions_sym]
            shape_function_derivative_single = [
                [float(derivative[i].subs(sub_dict)) for i in range(dim)]
                for derivative in self.shape_function_derivatives_sym]
            
            self.shape_functions_num.append(shape_function_single)
            self.shape_function_derivatives_num.append(shape_function_derivative_single)
                
        return self

    def compute_shape_functions_and_derivatives(self):

        """
        Compute the shape function values and derivatives at predefined natural integration points,
        as defined in section 4.1 of the Code_Aster documentation (https://code-aster.org/doc/v12/en/man_r/r3/r3.01.01.pdf).
        """
    
        if self.dimensions == "2d":

            a = 1.0 / 6.0
            b = 2.0 / 3.0

            self.natural_points = np.array([
                (a, a),
                (b, a),
                (a, b)
            ])
            self.compute_shape_functions_symbolically_cps6()

        elif self.dimensions == "3d":

            a = (5 - np.sqrt(5)) / 20.0
            b = (5 + 3 * np.sqrt(5)) / 20.0

            self.natural_points = np.array([
                (a, a, a),
                (a, a, b),
                (a, b, a),
                (b, a, a)
            ])
            self.compute_shape_functions_symbolically_c3d10()

        self.compute_shape_function_derivatives_symbolically()
        self.compute_shape_functions_and_derivatives_numerically()

        return self
