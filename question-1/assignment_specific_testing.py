from unittest import TestCase
import numpy as np
from gauss_elimination2 import *

class Test(TestCase):

    def test_guass_elimination_with_random_matrix(self):
        for n in range(100,1100,100):
            ## compute precision is up to 4 decimals and input matrix is also 4 decimal randam numbers
            ## without partial pivot by default
            gaussObject = GaussElimination(4, 4, False)
            ## random matrix with 100*101
            a_matrix = gaussObject.create_random_aug_matrix(n)
            xpp      = gaussObject.guass_elimination(a_matrix)
            pp_add   = gaussObject.additions_cnt
            pp_mul   = gaussObject.multiplications_cnt
            pp_div   = gaussObject.divisions_cnt
            #print(f"{n},{pp_add},{pp_mul},{pp_div}")
            # print(f"x_values={xpp}")
            # print(f"xx={gaussObject}")
            gaussObject.apply_partial_pivot = True
            gaussObject.reset_counts()
            xnpp = gaussObject.guass_elimination(a_matrix)
            print(f"n,pp.add,pp.mul,pp.div,nop.add,nod.mul,nod.dev")
            print(f"{n},{pp_add},{pp_mul},{pp_div},{gaussObject.additions_cnt},{gaussObject.multiplications_cnt},{gaussObject.divisions_cnt}")
