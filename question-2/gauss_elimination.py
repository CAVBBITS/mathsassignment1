# This modules consists functions to implement Gauss elimination with and without pivoting.
# It also counts the number of additions, multiplications
# and divisions performed during Gaussian elimination.
import numpy as np

class ERROR_IN_ELEMENTATION_STEPS(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ERROR_IN_BACK_SUBSTITUTION_STEPS(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def elimination(aug_matrix,n):
    """

    :param aug_matrix:
    :param n:
    :return:
    """
    ## for each row until n-1 rows.  This loop is n-1 rows because this is the row compared with the next row
    ## hence next row counter is till n.  However, in python range is not inclusive and hence need to give till n


    for k in range(0, n):
        m = k
        ## Here k starts from i+1 and ends till n because it gets compared with previous row
        for j in range(k + 1, n+1):
            ## This is for pivoting.  Generally previous row's (1,1) or (2,2) cell absolute value should be less than
            ## Next rows (2,1)  or (3,2) cells absolute value
            if (abs(aug_matrix[m, k]) < abs(aug_matrix[j, k])):
                m = j
        ## If one of the row's diagnal element is zero then there wont be unique solutions so quit
        print(f" at k={k} m={m} and j ={j} time aug_matrix={aug_matrix}")
        if aug_matrix[m, k] == 0:
            raise ERROR_IN_ELEMENTATION_STEPS(f" No unique solution exist as aug_matrix({m},{k}) is zero")
        ## Exchange the rows k and m
        temp_row = aug_matrix[m]
        aug_matrix[m] = aug_matrix[k]
        aug_matrix[k] = temp_row

        if aug_matrix[n, n] == 0:
            raise ERROR_IN_ELEMENTATION_STEPS(f" No unique solution exist as aug_matrix({n},{n}) is zero")
        print(f"aug_matrix={aug_matrix}")
        for l in range(k +1, n+1):
            coefficent_lk = aug_matrix[l, k] / aug_matrix[k, k]
            print(f"coefficent_lk={coefficent_lk} for row={k} ")
            for p in range(k , n + 2):
                print(f"for row {l} col {p} aug_matrix[l, p]={aug_matrix[l, p]}")
                aug_matrix[l, p] = aug_matrix[l, p] - (coefficent_lk * aug_matrix[k, p])
        if aug_matrix[n, n] == 0:
            raise ERROR_IN_ELEMENTATION_STEPS(f" No unique solution exist as aug_matrix({n},{n}) is zero")

    return aug_matrix


def back_substituion(aug_matrix,n):
    """
    back substitution finds the last row's x(n) value and uses the substitues that value in the previous equation to find
    previous rows x value i,e x(n-1)
    :param aug_matrix:  matrix after bring the matrix into upper triangle matrix, this is augmented matrix
    :param n: this is aug_matrix rows -1 because array index starts from 0
    :return:
    """
    x       =   np.zeros(n+1, dtype='float32')
    print(f"x={x}")
    x[n]    =   aug_matrix[n,n+1]/aug_matrix[n,n]
    for i in range(n-1, -1,-1):
        for j in range (i+1, n+1):
            coproduct = aug_matrix[i,j] * x[j]

        x[i] = (1/aug_matrix[i,i]) * (aug_matrix[i,n+1] -  coproduct)
    return x


def guass_elimination(aug_matrix:np.array)->np.array:
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

    n = aug_mat_rows-1 ## as matrix index starts from 0
    try:

        aug_matrix      = elimination(aug_matrix, n)
        output_x_vector = back_substituion(aug_matrix, n)
        return output_x_vector
    except Exception as e:
        print(f"Exception {e}")

    return None