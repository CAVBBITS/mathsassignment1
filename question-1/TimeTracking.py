# This programs finds the approximate time this computer takes for a single addition
# by adding first power(10,6) positive integers using a for loop and dividing the
# time taken by power(10,6). Similarly find the approximate time taken for a
# single multiplication and division. Report the result obtained in the
# form of a table.

# This program uses time module of python and calculates the difference between start time and end time
# before and after the operation

import time
from enum import Enum

class MathOperations(Enum):
    ADD         = 1
    SUBTRACT    = 2
    MULTIPLY    = 3
    DEVIDE      = 4



def trackSingleMathOperation(inta, intb, action):
    if action == MathOperations.ADD:
        startTime = time.time_ns()
        intc      = inta + intb
        endTime    = time.time_ns()
        print(f""" Time taken for action {inta} + {intb} = {intc} is {endTime-startTime} nano seconds """ )
    elif action == MathOperations.SUBTRACT:
        startTime = time.time_ns()
        intc = inta - intb
        endTime = time.time_ns()
        print(f""" Time taken for action {inta} - {intb} = {intc} is {endTime - startTime} nano seconds """)
    elif action == MathOperations.MULTIPLY:
        startTime = time.time_ns()
        intc = inta * intb
        endTime = time.time_ns()
        print(f""" Time taken for action {inta} * {intb} = {intc} is {endTime - startTime} nano seconds """)
    elif action == MathOperations.DEVIDE:
        startTime = time.time_ns()
        intc = inta / intb
        endTime = time.time_ns()
        print(f""" Time taken for action {inta} / {intb} = {intc} is {endTime - startTime} nano seconds """)

def trackMathOperationsForInput(inta, intb, action, base, pw):
    if action == MathOperations.ADD:
        startTime = time.time()
        for i in range(1, pow(base, pw)):
            intc = inta + intb
        endTime = time.time()
        print(f""" endTime={endTime} and startTime ={startTime} """)
        print(f""" Time taken for running {inta} + {intb} = {intc} for  {base}^{pw} times is {endTime - startTime}  seconds """)
    elif action == MathOperations.SUBTRACT:
        startTime = time.time()
        for i in range(1, pow(base, pw)):
            inta - intb
        endTime = time.time()
        print(
            f""" Time taken for running {inta} - {intb} for 10^6 times is {endTime - startTime}  seconds """)
    elif action == MathOperations.MULTIPLY:
        startTime = time.time()
        for i in range(1, pow(base, pw)):
            inta * intb
        endTime = time.time()
        print(
            f""" Time taken for running {inta} * {intb} for 10^6 times is {endTime - startTime}  seconds """)
    elif action == MathOperations.DEVIDE:
        startTime = time.time()
        for i in range(1, pow(base, pw)):
            inta / intb
        endTime = time.time()
        print(
            f""" Time taken for running {inta} / {intb} for 10^6 times is {endTime - startTime}  seconds """)

if __name__ == '__main__':
    trackSingleMathOperation(10,20,MathOperations.ADD)
    trackSingleMathOperation(10, 20, MathOperations.SUBTRACT)
    trackSingleMathOperation(10, 20, MathOperations.MULTIPLY)
    trackSingleMathOperation(10, 20, MathOperations.DEVIDE)

    trackMathOperationsForInput(10,20,MathOperations.ADD,10,6)
    trackMathOperationsForInput(10, 20, MathOperations.SUBTRACT,10,6)
    trackMathOperationsForInput(10, 20, MathOperations.MULTIPLY,10,6)
    trackMathOperationsForInput(10, 20, MathOperations.DEVIDE,10,6)
    #trackMathOperationsForInput(10, 6)
