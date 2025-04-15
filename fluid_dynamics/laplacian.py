import numpy as np
import scipy as sp


from .getCoeff import getCoeff as gc

def create_system(dom_1, num_1, cl_1):
    """
    Creates the system of equations for the Laplacian operator using sparse matrix techniques.

    Parameters
    ----------
    dom_1 : matrix
        Matrix representing the domain, where each element indicates whether the node is
        an interior point (1), a Dirichlet boundary (2), or outside the domain (0).
    num_1 : matrix
        Matrix containing the node numbers for each point in the domain.
    cl_1 : matrix
        Matrix of boundary conditions, where each element corresponds to the Dirichlet
        value at that node if it is a boundary point.

    Returns
    -------
    A_sparse : csr_matrix
        Sparse matrix representing the system of equations derived from the finite difference
        discretization of the Laplacian operator.
    B : array
        Right-hand side vector of the system of equations, incorporating boundary conditions.
    """
    num_max = int(np.max(num_1))
    num_min = int(np.min(num_1[num_1 != 0.0]))
    num = int(num_max - num_min + 1)

    
    rows = []
    cols = []
    data = []

    
    B = np.zeros(num)

    for i in range(num_min, num_max + 1):
        num_cent = i
        u, v = np.where(num_1 == num_cent)

        num_left = num_1[u, v - 1]
        num_right = num_1[u, v + 1]
        num_down = num_1[u + 1, v]
        num_up = num_1[u - 1, v]
        type_cent = dom_1[u, v]
        cl_cent = cl_1[u, v]

        j, a, b_val = gc(num_left, num_right, num_down, num_up, num_cent, type_cent, cl_cent)

        
        j_flat = j.ravel()
        a_flat = a.ravel()

        
        for jj, aa in zip(j_flat, a_flat):
            rows.append(i - num_min)
            cols.append(jj - num_min)
            data.append(aa)

        B[i - num_min] = b_val

    
    A_sparse = sp.sparse.csc_matrix((data, (rows, cols)), shape=(num, num)).tocsr()
    return A_sparse, B



def min_exclude_zero(arr):
    """
    Finds the minimum value in an array excluding zeros.

    Parameters
    ----------
    arr : array-like
        Input array containing node numbers, including zeros for non-domain points.

    Returns
    -------
    min_val : int
        The minimum non-zero value in the array. Returns None if all elements are zero.
    """
    non_zero = arr[arr != 0.0]
    if len(non_zero) == 0:
        return None  
    return np.min(non_zero)

def solve_system(A, B):
    """
    Solves the sparse linear system AX = B using efficient sparse matrix methods.

    Parameters
    ----------
    A : array
        Sparse matrix of the system of equations derived from the finite difference
        discretization of the Laplacian operator.
    B : array
        Right-hand side vector incorporating boundary conditions.

    Returns
    -------
    X : array
        Solution vector representing the field variable (e.g., stream function or potential)
        at each node in the computational domain.
    """
    
    X = sp.sparse.linalg.spsolve(A, B)
    return X

