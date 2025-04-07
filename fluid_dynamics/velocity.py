from tools.deriv import deriv as d
import numpy as np
import scipy as sp


def min_exclude_zero(arr):
    non_zero = arr[arr != 0.0]
    if len(non_zero) == 0.0:
        return None  # ou une autre valeur par défaut
    return np.min(non_zero)

def velocity(X, dom, h):
    # Initialize velocity arrays with same shape as domain
    u = np.zeros_like(X, dtype=float)
    v = np.zeros_like(X, dtype=float)
    
    # Calculate velocities at each point in the domain
    for i in range(1, X.shape[0]-1):
        for j in range(1, X.shape[1]-1):
            if dom[i,j] != 0:  # Only calculate for points in the domain
                # Get node types for derivatives
                type_left = dom[i,j-1]
                type_right = dom[i,j+1] 
                type_up = dom[i-1,j]
                type_down = dom[i+1,j]
                type_cent = dom[i,j]
                
                # Get stream function values
                psi_left = X[i,j-1]
                psi_right = X[i,j+1]
                psi_up = X[i-1,j] 
                psi_down = X[i+1,j]
                psi_cent = X[i,j]
                
                # Calculate velocity components using derivatives of stream function
                # u = ∂ψ/∂y
                u[i,j] = d(psi_up, psi_cent, psi_down,
                          type_up, type_cent, type_down, h)
                
                # v = -∂ψ/∂x
                v[i,j] = -d(psi_left, psi_cent, psi_right,
                           type_left, type_cent, type_right, h)
    
    return u, v