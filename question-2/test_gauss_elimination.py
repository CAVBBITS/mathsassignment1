from unittest import TestCase
import numpy as np
from gauss_elimination import *

class Test(TestCase):
    def test_elimination(self):
        aug_matrix  =   np.array([[1,2,7],[3,4,8]])
        n           = len(aug_matrix) - 1
        print(f"length ={n}")
        #print(f"matrix={aug_matrix}")
        aug_matrix  = elimination(aug_matrix,n)
        print(aug_matrix)

    def test_elimination(self):
        aug_matrix  =   np.array([[1,2,7],[1,2,8]])
        n           = len(aug_matrix) - 1
        print(f"length ={n}")
        #print(f"matrix={aug_matrix}")
        aug_matrix  = elimination(aug_matrix,n)
        print(aug_matrix)

    def test_back_substituion(self):
        aug_matrix  =   np.array([[1,2,7],[0,-2,-13]])
        n           =   len(aug_matrix) - 1
        x           =   back_substituion(aug_matrix,n)
        print(x)

    def test_guass_elimination(self):
        aug_matrix  =   np.array([[1,2,7],[3,4,8]])
        x           =   guass_elimination(aug_matrix)
        print(x)

    def test_guass_elimination2(self):
        aug_matrix  =   np.array([[1,2,7],[1,2,8]])
        x           =   guass_elimination(aug_matrix)
        print(x)


