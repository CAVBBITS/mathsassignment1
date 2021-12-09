# This modules consists functions to implement guass elimination with and without pivoting.
# It also counts the number of additions, multiplications
# and divisions performed during guassian elimination.
import numpy as np
import random

class ERROR_IN_PARTITIAL_PIVOT_STEPS(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class ERROR_IN_ELEMENTATION_STEPS(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ERROR_IN_BACK_SUBSTITUTION_STEPS(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class GaussElimination:
    """
        GaussElimination is the class which initializes parameters.
        The constructor expects element_precision and  compute_precision.  Default values for both 5
        element_precision: used while creating the ramdom matrix
        compute_precision : used while calculating the matrix

    """
    def __init__(self, compute_precision=5,element_precision=5,apply_partial_pivot=True):
        self.precision = compute_precision
        self.element_precision = element_precision
        self.apply_partial_pivot = apply_partial_pivot
        self.additions_cnt         = 0
        self.multiplications_cnt   = 0
        self.divisions_cnt         = 0
    def reset_counts(self):
        self.additions_cnt         = 0
        self.multiplications_cnt   = 0
        self.divisions_cnt         = 0

    def __str__(self) -> str:
        output = f""" additions_cnt={self.additions_cnt}
                      multiplications_cnt={self.multiplications_cnt}
                      divisions_cnt={self.divisions_cnt}
                      apply_partial_pivot={self.apply_partial_pivot}
                      precision={self.precision}
                      element_precision={self.element_precision}
                    """
        return (output)

    def partial_pivot(self,a_matrix:np.array, head_row_index:int)->np.array:
        """
        partial_pivot methods applies row pivot on a_matrix by comparing curr_row_index with subsequent rows
        and swaps the rows if the a_matrix[curr_row_index,column_index] with a_matrix[max_row_index,column_index]
        :param a_matrix:
        :param head_row_index:
        :return:
        """
        max_value_row_index = head_row_index
        ## This is same as row_index because, the guass elimination works with nxn matrix
        head_col_index     = head_row_index
       # print(f"head_row_index={head_row_index} and length = {len(a_matrix)} range = {range(head_row_index + 1, 1, len(a_matrix)+1)}")
        for next_row_index in range(head_row_index + 1, len(a_matrix),1):
        #    print(f" printing next_row_index ={next_row_index} ")
            if a_matrix[max_value_row_index, head_col_index] < a_matrix[next_row_index, head_col_index]:
                max_value_row_index = next_row_index

        if max_value_row_index>head_row_index :
            ## print(f"head_row_index={head_row_index} and max_value_row_index={max_value_row_index}")
            ## need to create a new array instead of refering the index.  Refering the index is like pass by reference and
            ## that would not swap
            temp                            =   np.array(a_matrix[max_value_row_index])
            a_matrix[max_value_row_index]   =   a_matrix[head_row_index]
            a_matrix[head_row_index]        =   temp
            #
            print("partial pivot is applied")

        if a_matrix[max_value_row_index, head_col_index] == 0:
            raise  ERROR_IN_PARTITIAL_PIVOT_STEPS(f" Multiple solutions exist as augmented matrix = {a_matrix} "
                                                  "row {max_value_rc_index}  and {curr_row_index} is zero")
        return a_matrix

    def elimination(self,a_matrix):
        """
        :param a_matrix:
        :param n:
        :return:
        """
        ## for each row until n-1 rows.  This loop is n-1 rows because this is the row compared with the next row
        ## hence next row counter is till n.  However, in python range is not inclusive and hence need to give till n


        for i in range(0, len(a_matrix) - 1, 1):
            cc = i

            if(self.apply_partial_pivot):
                a_matrix = self.partial_pivot(a_matrix,i)
            for nr in range (i+1, len(a_matrix), 1):
                coefficient = round(a_matrix[nr][cc] / a_matrix[i][cc],self.precision)
                self.divisions_cnt+=1
                for cc1 in range(cc, len(a_matrix[0]), 1):
                    a_matrix[nr, cc1] = round(a_matrix[nr, cc1] - (coefficient * a_matrix[i, cc1]),self.precision)
                    self.multiplications_cnt+=1
                    self.additions_cnt+=1


        if a_matrix[len(a_matrix) - 1, len(a_matrix) - 1] == 0:
            raise ERROR_IN_ELEMENTATION_STEPS(f" Matrix nxn cell is 0 and hence solution is not possible.  "
                                              f" Matrix = {a_matrix}"
                                              f" aug_matrix[{len(a_matrix) - 1},{len(a_matrix) - 1}] = "
                                              f" {a_matrix[len(a_matrix) - 1, len(a_matrix) - 1]}")


        return a_matrix


    def back_substitution(self,a_matrix):
        """
        back substitution finds the last row's x(n) value and uses the substitues that value in the previous equation to find
        previous rows x value i,e x(n-1)
        :param a_matrix:  matrix after bring the matrix into upper triangle matrix, this is augmented matrix
        :param n: this is aug_matrix rows -1 because array index starts from 0
        :return:
        """
        n       = len(a_matrix)
        x       = np.zeros(n, dtype='float32')
        x[n-1]  = round(a_matrix[n-1, n] / a_matrix[n-1, n-1],self.precision)
        self.divisions_cnt+=1
        ## remember x[n-1] already filled so start from n-2


        for i in range(n-2, -1,-1):
            coproduct = 0
            for j in range (i+1, n):
                coproduct = round(coproduct + a_matrix[i, j] * x[j],self.precision)
                self.multiplications_cnt += 1
                self.additions_cnt+1
                # print(f"coproduct={coproduct} i={i} j={j} x[j={j}] = {x[j]} and a_matrix[i={i}, j={j}] = {a_matrix[i, j]}")
            x[i] = round((1 / a_matrix[i, i]) * (a_matrix[i, n] - coproduct),self.precision)
            self.divisions_cnt += 1
            self.multiplications_cnt += 1
            self.additions_cnt + 1
        return x



    def guass_elimination(self,aug_matrix:np.array)->np.array:
        """
        guass_elimination first implements eliminates the row elements to get the upper triangle format
        it then applies back substitution and solves the equation.  It is expected that the input matrix
        is nxn matrix and augumented matrix is nxn+1 matrix
        :param aug_matrix:  Augmented matrix is Matrix + R.H.S column vector.
        :param aug_mat_rows: This is same as input matrix number of rows
        :param aug_mat_cols: This is augmented i.,e  input may have n columns then add RHS vector and it becomes n+1 columns
        """
        aug_mat_rows = len(aug_matrix)
        aug_mat_cols = len(aug_matrix[0])
        if (aug_mat_cols - 1 != aug_mat_rows):
            raise " Not a valid matrix.  Matrix must be in nxn and augumented matrix must be nxn+1"
        # declare few variables to process the algorthm

        try:
            aug_matrix      = self.elimination(aug_matrix)
            output_x_vector = self.back_substitution(aug_matrix)
            return output_x_vector
        except Exception as e:
            print(f"Exception {e}")

        return None


    def create_random_aug_matrix(self,n:int)->np.array:
        """
        This creates a random matrix with n - rows   and n+1 columns with given precision
        :param n:
        :return:
        """
        data        = []
        for x in range(0,n*(n+1)):
            data.append(round(random.uniform(0, 100),self.element_precision))
        a_matrix    = np.array(data).reshape(n,n+1)
        return a_matrix

