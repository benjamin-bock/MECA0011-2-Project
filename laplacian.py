#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 18:53:03 2025

@author: benjamin-bock
"""

def create_system(dom_1, num_1, cl_1):
    """
    Create the system of equations for the Laplacian operator.
    
    Parameters
    ----------
    dom_1 : matrix
        Matrix of the domain.
    num_1 : matrix
        Matrix of the number of the nodes in the domain.
    cl_1 : list
        List of the boundary conditions.
        
    Returns
    -------
    A : array
        The matrix of the system of equations.
    b : array
        The right-hand side of the system of equations.
    """
    import numpy as np
    
    # Create the matrix A
    A = np.zeros((num_1, num_1))
    
    # Create the vector b
    b = np.zeros(num_1)
    
    # Loop over the nodes
    for i in range(num_1):
        # Check if the node is on the boundary
        if i in cl_1:
            # Set the diagonal element to 1
            A[i, i] = 1
            # Set the right-hand side to the boundary condition
            b[i] = 0
        else:
            # Set the diagonal element to -4
            A[i, i] = -4
            # Set the off-diagonal elements to 1
            if i - 1 >= 0:
                A[i, i - 1] = 1
            if i + 1 < num_1:
                A[i, i + 1] = 1
            if i - dom_1[1] >= 0:
                A[i, i - dom_1[1]] = 1
            if i + dom_1[1] < num_1:
                A[i, i + dom_1[1]] = 1
    
    return A, b