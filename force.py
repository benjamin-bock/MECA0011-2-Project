import numpy as np
import scipy as sp

def force(p,x,y):
    fx = 0.0
    fy = 0.0
    n = len(x)
    
    for i in range(n):
        j = (i + 1) % n 
        
        dx = x[j] - x[i]
        dy = y[j] - y[i]
        
        p_avg = 0.5 * (p[i] + p[j])
        
        #Aucune idéé de pourquoi il faut inverser les signes et les dx/dy, ehhh merceee le test du prof
        fx += p_avg * dy
        fy += p_avg * -dx
    return fx,fy