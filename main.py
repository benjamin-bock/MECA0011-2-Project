import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from laplacian import create_system

def main():

    dom_1 = np.loadtxt('CL\1-dom.txt', dtype = int)
    cl_1 = np.loadtxt('CL\1-cl.txt', dtype = int)
    num_1 = np.loadtxt('CL\1-num.txt', dtype = int)

    create_system(dom_1, num_1, cl_1)