import numpy as np
import scipy as sc
from getCoeff import getCoeff as gc
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

def Mat(h, cl, dom, num):
    data = []
    row = []
    cols = []
    nb_noeud = (dom[-1] - dom[0])/h
    B = np.zeros((nb_noeud))
    for i in range(nb_noeud):
        j, a, b = gc(num_left, num_right, num_down, num_up, num_cent, type_cent, cl_cent)
        for j_col and a_val in zip(j,a):
            rows.append(i)
            cols.append(j_col -1)
            data.append(a_val)
            B[i] = b
    matrix = csr_matrix((data, (row, row)), shape=(nb_noeud, nb_noeud))
    return matrix, B