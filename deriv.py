import numpy as np
import scipy as sp

def deriv(f_left, f_c, f_right, type_left, type_c, type_right, h) :
    if type_c == 0:
        raise ValueError("Invalid value for type_cent: 0 is not in the domain of calcul.")
    if type_c == 1:
        v = (f_right - f_left) / (2*h)
    if type_c == 2:
        if type_left == 0:
            v = (f_right - f_c) / h
        if type_right == 0:
            v = (f_c - f_left) / h
    return v