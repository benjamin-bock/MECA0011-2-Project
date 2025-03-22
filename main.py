import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

import laplacian as lp

def main():

    dom_1 = np.loadtxt('CL/1-dom.txt', dtype = int)
    cl_1 = np.loadtxt('CL/1-cl.txt', dtype = float)
    num_1 = np.loadtxt('CL/1-num.txt', dtype = int)

    A, B = lp.create_system(dom_1, num_1, cl_1)
    print(f"La matrice A vaut \n {A} \n")
    print(f"La matrice B vaut \n {B} \n")
    
    X = lp.solve_system(A, B)
    print(f"La matrice X vaut \n {X} \n")
    
    return X
