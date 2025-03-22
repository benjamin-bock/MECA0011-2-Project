import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from laplacian import create_system

def main():

    dom_1 = np.loadtxt('CL/1-dom.txt', dtype = int)
    cl_1 = np.loadtxt('CL/1-cl.txt', dtype = float)
    num_1 = np.loadtxt('CL/1-num.txt', dtype = int)

    A, B = create_system(dom_1, num_1, cl_1)

    print(A)
    print(B)
    

    return A, B
