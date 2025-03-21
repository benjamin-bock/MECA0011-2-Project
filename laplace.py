import numpy as np
import scipy as sc
from getCoeff import getCoeff as gc
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def laplace(n, cl, f):

    my_array = np.loadtxt('CL\1-num.txt', dtype = int)