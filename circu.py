import numpy as np
import scipy as sp

def circu(u, v, x, y):
    c = 0.0
    n = len(x)
    
    for i in range(n):
        j = (i + 1) % n  
        
        dx = x[j] - x[i]
        dy = y[j] - y[i]
    
        u_avg = 0.5 * (u[i] + u[j])
        v_avg = 0.5 * (v[i] + v[j])
        
        c += u_avg * dx + v_avg * dy
    return c