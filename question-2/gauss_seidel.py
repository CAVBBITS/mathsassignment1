# This modules consists functions to implement guass elimination with and without pivoting.
# It also counts the number of additions, multiplications
# and divisions performed during guassian elimination.
import random

import numpy as np



class ERROR_IN_DIAGNALLY_DOMINANT_CONVERSION(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class MATRIX_IS_NOT_DIAGONALLY_DOMINANT(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class GaussSeidel:
    """
        GaussElimination is the class which initializes parameters.
        The constructor expects element_precision and  compute_precision.  Default values for both 5
        element_precision: used while creating the ramdom matrix
        compute_precision : used while calculating the matrix

    """
    def __init__(self, compute_precision=5, element_precision=5, compare_precision=3,apply_rows_interchange=True):
        self.precision                  = compute_precision
        self.element_precision          = element_precision
        self.compare_precision          = compare_precision
        self.apply_rows_interchange     = apply_rows_interchange
        self.additions_cnt         = 0
        self.multiplications_cnt   = 0
        self.divisions_cnt         = 0
        self.number_of_iterations  = 0

    def __str__(self) -> str:
        output = f""" additions_cnt={self.additions_cnt}
                      multiplications_cnt={self.multiplications_cnt}
                      divisions_cnt={self.divisions_cnt}
                      apply_partial_pivot={self.apply_rows_interchange}
                      precision={self.precision}
                      element_precision={self.element_precision}
                      compare_precision={self.compare_precision}
                      number_of_iterations={self.number_of_iterations}
                    """
        return (output)

    def compare_vectors(self, v1,v2)->bool:

        if len(v1)!=len(v2) :
            return False

        for i in range (0,len(v1)):
            if round(v1[i,0],self.compare_precision) != round(v2[i,0],self.compare_precision) :
                return False
        return True

    def get_matrix_fobo_norm(self,matrix)->float:
        """
        get_matrix_fobo_norm calculates the fobo norm
        :param matrix:
        :return:
        """
        rows = len(matrix)
        cols = len(matrix[0])
        norm =0.0
        for i in range(0, rows):
            for j in range(0,cols):
                norm += pow(matrix[i,j],2)
        norm = pow(norm,0.5)
        return norm


    def is_nn_matrix_to_diagnolly_dominant(matrix: np.array) -> bool:
        """
         is_nn_matrix_to_diagnolly_dominant checks whether given matrix is diagonally dominant
         Throws an exception if the given matrix is not nXn matrix
        :param matrix:
        :return:
        """
        rows_size = len(matrix)
        col_szie = len(matrix[0])
        if (rows_size != col_szie):
            raise ERROR_IN_DIAGNALLY_DOMINANT_CONVERSION(f"""Matrix must be a square matrix. 
             Given matrix {matrix} is not square matrix rows_size = {rows_size} and col_size={col_szie}""")

        for i in range(0, rows_size):
            row_magnitude_sum = 0.0
            for j in range(0, col_szie):
                if i != j:
                    row_magnitude_sum += abs(matrix[i, j])
            if (row_magnitude_sum > matrix[i, i]):
                print(f""" Given matrix {matrix} is not diagnolly dominant at row {i} and rows values {matrix[i]}""")
                return False
        return True


    def guass_seidel_iterate_function(self, matrix, rhs_vector):
        """
         guass_seidel_iterate_function, applies the iterate model calculates the final x values
         The iterations continue till the current and previous x values are same with defined precision
        :param matrix:
        :param rhs_vector:
        :return:
        """
        n           = len(matrix)
        prev_vals   = np.zeros(n, dtype='float32')

        while True:
            curr_vals = np.zeros(n, dtype='float32')
            for i in range(0,n):
                coproduct = 0.0
                for j in range(0,n):
                    if(i!=j):
                        ## mostly the change between seidel and gauss is substitution of current iteration values itself
                        coproduct += (matrix[i,j])*curr_vals[i]
                        self.additions_cnt+=1
                        self.multiplications_cnt+=1
                curr_vals[i] = round((rhs_vector[i,0] - coproduct)/matrix[i,i],self.precision)
                self.additions_cnt+=1
                self.divisions_cnt+=1
            norm = self.get_matrix_fobo_norm(curr_vals)
            print(f""" seidel for iteration {self.number_of_iterations} norm is = {norm} and matrix = {curr_vals} """)
            self.number_of_iterations += 1

            if self.compare_vectors(prev_vals,curr_vals):
                return curr_vals
            prev_vals = curr_vals[i,i]

        ## mostly the program doesn reach this statements
        return curr_vals;





    def guass_seidel(self, matrix:np.array,rhs_vector:np.array)->np.array:
        """
        guass_seidel, this first checks whether the matrix is diagonally dominant or not
        and if the matrix is diagonally dominant then it calculates values
        :param matrix:
        :param rhs_vector:
        :return:
        """
        rows            = len(matrix)
        cols            = len(matrix[0])
        rhs_vector_len  = len(rhs_vector)
        if (cols  != rows):
            raise " Not a valid matrix.  Matrix must be in nxn matrix"
        if(rows != rhs_vector_len):
            raise f"""Matrix and rhs vector (RHS scalar values) should be of same number of rows 
                        matrix rows={rows}  rhs_ventor_len={rhs_vector_len}"""

        # declare few variables to process the algorthm

        try:
            if self.is_nn_matrix_to_diagnolly_dominant(matrix) == False:
                raise MATRIX_IS_NOT_DIAGONALLY_DOMINANT(f"Matrix {matrix} is not diagonally dominant")

            matrix_fobo_norm = self.get_matrix_fobo_norm(matrix)
            print(f""" Seidel input matrix={matrix} and matrix_fobo_norm={matrix_fobo_norm}""")
            final_x_vector  = self.guass_seidel_iterate_function(matrix,rhs_vector)
            return final_x_vector
        except Exception as e:
            print(f"Exception {e}")

        return None


    def create_random_matrix(self,n:int)->np.array:
        """
        This creates a random matrix with n - rows   and n+1 columns with given precision
        :param n:
        :return:
        """
        data        = []
        for x in range(0,n*n):
            data.append(round(random.uniform(0, 100),self.element_precision))
        matrix    = np.array(data).reshape(n,n)
        return matrix

    def create_rhs_value_vector(self, n:int)->np.array:
        """

        :param n:
        :return:
        """
        data = []
        for x in range(0, n):
            data.append(round(random.uniform(0, 100), self.element_precision))
        matrix = np.array(data)
        return matrix
