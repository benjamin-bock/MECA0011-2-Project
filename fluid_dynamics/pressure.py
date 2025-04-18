import numpy as np
from tools.constante import Q_out

def calculate_pressure(X, Y, U, V,dom) :
    
    C = Q_out
    rho = 1000
    g = 9.81
    
    p = np.zeros((X.shape[0], X.shape[1]))
    ampl = 0
    
    for i in range(X.shape[0]) :
        for j in range(X.shape[1]) :
                if(dom[i][j] == 0):
                     p[i][j] = np.nan
                else : 
                    ampl = (U[i][j])**2 + (V[i][j])**2 
                    p[i][j] = rho*g*(C- ampl/(2*g))
    
    return p