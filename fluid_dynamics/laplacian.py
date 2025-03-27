#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 18:53:03 2025

@author: benjamin-bock
"""

import numpy as np
import scipy as sp

from fluid_dynamics.getCoeff import getCoeff as gc
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
    B : array
        The right-hand side of the system of equations.
    """
    
    
    # Create the matrix A
    num_max = int(np.max(num_1)) #Numéro de noeud max
    num_min = int(min_exclude_zero(num_1)) #Numéro de noeud min 
    num = int(num_max - num_min + 1) #Nombre de noeuds
    A = np.zeros((num, num)) #Matrice A de taille num x num
    
    # Create the vector b
    B = np.zeros((num)) #Vecteur b de taille num
    
    # Loop over the nodes
    for i in range(num_min, num_max + 1):

        num_cent = i #Numéro du noeud central
        u, v = np.where(num_1 == num_cent) #Coordonnées du noeud central

        num_left = num_1[u, v - 1] #Numéro du noeud à gauche
        num_right = num_1[u, v + 1] #Numéro du noeud à droite
        num_down = num_1[u + 1, v] #Numéro du noeud en bas
        num_up = num_1[u - 1, v] #Numéro du noeud en haut
        type_cent = dom_1[u, v] #Type du noeud central
        cl_cent = cl_1[u, v] #Condition limite du noeud central

        j, a, b = gc(num_left, num_right, num_down, num_up, num_cent, type_cent, cl_cent)
        A[i - num_min, j - num_min] = a
        B[i - num_min] = b
    return A, B

def min_exclude_zero(arr):
    non_zero = arr[arr != 0.0]
    if len(non_zero) == 0.0:
        return None  # ou une autre valeur par défaut
    return np.min(non_zero)

def solve_system(A, B):
    """
    Solve the system of equations AX = B.
    
    Parameters
    ----------
    A : array
        The matrix of the system of equations.
    B : array
        The right-hand side of the system of equations.
        
    Returns
    ------- 
    X : array 
        The solution of the system of equations. Either \phi or \psi wether it is a potential or a stream function.
    """
    # Convertir A en format CSR
    A_csr = sp.sparse.csr_matrix(A)
    X = sp.sparse.linalg.spsolve(A_csr, B)
    return X
