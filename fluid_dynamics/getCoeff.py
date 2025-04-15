import numpy as np 


def getCoeff(num_left, num_right, num_down, num_up, num_cent, type_cent, cl_cent):
    """
    Computes the coefficients for a finite difference scheme in a computational grid.

    This function calculates the coefficients for the matrix equation arising from
    the discretization of a partial differential equation using the finite difference
    method. It takes into account the type of the central node (whether it's a
    Dirichlet boundary condition or an interior point) and constructs the
    corresponding row of the coefficient matrix and the right-hand side vector.

    Parameters
    ----------
    num_left, num_right, num_down, num_up, num_cent : int
        Grid node indices in the specified directions relative to the central node.
    type_cent : int
        Type of the central node (1 for interior point, 2 for Dirichlet boundary).
    cl_cent : float
        Dirichlet boundary condition value for the central node if type_cent is 2.

    Returns
    -------
    j : numpy.ndarray
        Column vector containing the column indices of the non-zero coefficients.
    a : numpy.ndarray
        Column vector containing the non-zero coefficients of the matrix row.
    b : float
        Value of the right-hand side term in the equation.

    Raises
    ------
    ValueError
        If type_cent is 0, which is not a valid node type in the computational domain.
    """
    
    if type_cent == 0:
        raise ValueError("Invalid value for type_cent: 0 is not in the domain of calcul.")
    
    
    num_left = int(num_left)
    num_right = int(num_right)
    num_down = int(num_down)
    num_up = int(num_up)
    num_cent = int(num_cent)
    
    if type_cent == 2:
        j = np.array([num_cent]).reshape(-1,1)
        a = np.array([1]).reshape(-1,1)
        b = cl_cent
    
    if type_cent == 1:
        j = np.array([num_left, num_right, num_down, num_up, num_cent]).reshape(-1,1)  
        a = np.array([1, 1, 1, 1, -4]).reshape(-1,1)
        b = 0
    
    return j, a, b