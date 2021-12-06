from unittest import TestCase
import numpy as np
from gauss_elimination2 import *

class Test(TestCase):



    def test_partitial_pivot(self):
        ran = range(0,10,-1)
        for i in ran:
            print(i)

        a_matrix = np.array([[1.0, 2.0, 7.0], [3.0, 4.0, 8.0]],dtype=float)
        gaussObject = GaussElimination(3)
        a_matrix = gaussObject.partial_pivot(a_matrix, 0)
        print(a_matrix)

    def test_elimination(self):
        aug_matrix  =   np.array([[3,4,8],[1,2,7]],dtype=float)
        gaussObject = GaussElimination(2)
        aug_matrix  = gaussObject.elimination(aug_matrix)
        print(aug_matrix)

    def test_elimination_many_solution(self):
        aug_matrix  =   np.array([[1,2,7],[1,2,8]],dtype=float)
        gaussObject = GaussElimination()
        aug_matrix  = gaussObject.elimination(aug_matrix)
        print(aug_matrix)

    def test_back_substituion(self):
        aug_matrix  =   np.array([[1,2,7],[0,0.667,4.333]],dtype=float)
        gaussObject = GaussElimination(2)
        x           =   gaussObject.back_substitution(aug_matrix)
        print(x)

    def test_back_substituion_2(self):
        aug_matrix  =   np.array([[3,4,8,9],[0,1,3,7],[0,0,1,2]],dtype=float)
        gaussObject = GaussElimination(4)
        x           =   gaussObject.back_substitution(aug_matrix)
        print(x)

    def test_guass_elimination(self):
        aug_matrix  =   np.array([[1,2,7],[3,4,8]],dtype=float)
        gaussObject = GaussElimination(4)
        x           =  gaussObject.guass_elimination(aug_matrix)
        print(x)

    def test_guass_elimination_error_test(self):
        aug_matrix  =   np.array([[1,2,7],[1,2,8]],dtype=float)
        gaussObject = GaussElimination(4)
        x           = gaussObject.guass_elimination(aug_matrix)
        print(x)

    def test_guass_elimination_3(self):
        aug_matrix  =   np.array([[1,2,-1,6],[3,8,9,10],[2,-1,2,-2]],dtype=float)
        gaussObject =   GaussElimination(2)
        x           =   gaussObject.guass_elimination(aug_matrix)
        print(x)

    def test_guass_elimination_4(self):
        aug_matrix  =   np.array([[1,2,3,18],[2,1,-4,-30],[-5,8,17,96]],dtype=float)
        gaussObject =   GaussElimination()
        x           =   gaussObject.guass_elimination(aug_matrix)
        print(x)
