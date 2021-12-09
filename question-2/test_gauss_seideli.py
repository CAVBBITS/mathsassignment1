from unittest import TestCase
from  gauss_seidel import *
import numpy as np

class Test(TestCase):

    def test_guass_seidel_iterate_function_2(self):
        gaussJacobi = Gausssiedel(3, 0, 0.001)
        matrix = np.array([86, 20, 33.0, 3.0, 91, 496, 53, 97, 19, 88, 280, 95.0, 83., 33., 8., 124.]).reshape(4, 4)
        rhs_vector = np.array([63.0, 24.0, 90.0, 29.0]).reshape(4, 1)
        x_vector   = gaussJacobi.guass_siedel_iterate_function(matrix,rhs_vector)
        print(f"x_vector={x_vector}")

    def test_guass_seidel_iterate_function_3(self):
        gaussJacobi = Gausssiedel(3, 0, 0.001)
        matrix = np.array([2, 0, 0.0, 0.0, 0, 2, 0, 0, 0, 0, 3, 0.0, 0., 0., 0., 2.]).reshape(4, 4)
        rhs_vector = np.array([1,1,1,1.0]).reshape(4, 1)
        x_vector   = gaussJacobi.guass_siedel_iterate_function(matrix,rhs_vector)
        print(f"x_vector={x_vector}")

