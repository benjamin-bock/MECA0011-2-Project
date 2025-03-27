from tools.deriv import deriv as d
import numpy as np
import scipy as sp


def min_exclude_zero(arr):
    non_zero = arr[arr != 0.0]
    if len(non_zero) == 0.0:
        return None  # ou une autre valeur par défaut
    return np.min(non_zero)

def velocity(X, dom, cl,num,h ):
    col = X.shape[1]
    row = X.shape[0]
    num_max = int(np.max(num)) #Numéro de noeud max
    num_min = int(min_exclude_zero(num)) #Numéro de noeud min 
    num = int(num_max - num_min + 1)
    V = np.zeros_like(X, dtype=float)
    for i in range(1, X.shape[0]-1):
        num_cent = i
        u, v = np.where(num == num_cent)
        num_left = num[u, v - 1]
        num_right = num[u, v + 1]
        num_down = num[u + 1, v]
        num_up = num[u - 1, v]
        type_cent = dom[u, v]
        type_left = dom[u, v - 1]
        type_right = dom[u, v + 1]
        cl_cent = cl[u, v]
        V[u, v] = d(X[u, v - 1], X[u, v], X[u, v + 1], type_left,type_cent, type_right, h)
    return V