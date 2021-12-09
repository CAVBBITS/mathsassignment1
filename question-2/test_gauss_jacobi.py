from unittest import TestCase
from  gauss_jacobi import *
import numpy as np

class Test(TestCase):

    def test_create_random_matrix(self):
        gaussJacobi = GaussJacobi (3,0,0.001)
        matrix = gaussJacobi.create_random_matrix(4)
        print(f"matrix={matrix}")

    def test_create_rhs_value_vector(self):
        gaussJacobi = GaussJacobi (3,0,0.001)
        rhs_value = gaussJacobi.create_rhs_value_vector(4)
        print(f"matrix={rhs_value}")

    def test_matrix_norm1(self):
        gaussJacobi = GaussJacobi (3,0,0.001)
        matrix      = np.array([86,20,33.0,3.0,91,96,53,97,19,88,80,95.0,83., 33.,  8.,  4.]).reshape(4,4)
        rhs_vector  = np.array([63.0, 24.0, 90.0, 29.0]).reshape(4,1)
        norm        = gaussJacobi.get_matrix_fobo_norm(matrix)
        rhs_norm    = gaussJacobi.get_matrix_fobo_norm(rhs_vector)
        print(f"norm={norm} and rhs_norm={rhs_norm}")
        assert(norm==264.947)
        assert(rhs_norm==116.129)

    def test_matrix_norm2(self):
        gaussJacobi = GaussJacobi (3,0,0.001)
        matrix      = np.array([86,20,33.0,3.0,91,96,53,97,19,88,80,95.0,83., 33.,  8.,  4.]).reshape(4,4)
        isDominent  = gaussJacobi.is_nn_matrix_to_diagnolly_dominant(matrix)
        print(f"isDominent={isDominent}")
        assert(isDominent==False)

    def test_matrix_norm3(self):
        gaussJacobi = GaussJacobi(3, 0, 0.001)
        matrix = np.array([86, 20, 33.0, 3.0, 91, 496, 53, 97, 19, 88, 280, 95.0, 83., 33., 8., 124.]).reshape(4, 4)
        isDominent = gaussJacobi.is_nn_matrix_to_diagnolly_dominant(matrix)
        print(f"isDominent={isDominent}")
        assert (isDominent == True)

    def test_guass_jacobi_iterate_function(self):
        gaussJacobi = GaussJacobi(3, 0, 0.001)
        matrix = np.array([86, 20, 33.0, 3.0, 91, 496, 53, 97, 19, 88, 280, 95.0, 83., 33., 8., 124.]).reshape(4, 4)
        rhs_vector = np.array([63.0, 24.0, 90.0, 29.0]).reshape(4, 1)
        x_vector   = gaussJacobi.guass_jacobi_iterate_function(matrix,rhs_vector)
        print(f"x_vector={x_vector}")

    def test_guass_jacobi_iterate_function_2(self):
        gaussJacobi = GaussJacobi(3, 0, 0.001)
        matrix = np.array([2, 0, 0.0, 0.0, 0, 2, 0, 0, 0, 0, 3, 0.0, 0., 0., 0., 2.]).reshape(4, 4)
        rhs_vector = np.array([1,1,1,1.0]).reshape(4, 1)
        x_vector   = gaussJacobi.guass_jacobi_iterate_function(matrix,rhs_vector)
        print(f"x_vector={x_vector}")

    def test_guass_seidel_iterate_function_2(self):
        gaussJacobi = Gausssiedel(3, 0, 0.001)
        matrix = np.array([86, 20, 33.0, 3.0, 91, 496, 53, 97, 19, 88, 280, 95.0, 83., 33., 8., 124.]).reshape(4, 4)
        rhs_vector = np.array([63.0, 24.0, 90.0, 29.0]).reshape(4, 1)
        x_vector   = gaussJacobi.guass_jacobi_iterate_function(matrix,rhs_vector)
        print(f"x_vector={x_vector}")

    def test_guass_seidel_iterate_function_2(self):
        gaussJacobi = Gausssiedel(3, 0, 0.001)
        matrix = np.array([2, 0, 0.0, 0.0, 0, 2, 0, 0, 0, 0, 3, 0.0, 0., 0., 0., 2.]).reshape(4, 4)
        rhs_vector = np.array([1,1,1,1.0]).reshape(4, 1)
        x_vector   = gaussJacobi.guass_jacobi_iterate_function(matrix,rhs_vector)
        print(f"x_vector={x_vector}")

