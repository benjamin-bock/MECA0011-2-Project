import numpy as np

def calculate_pressure(X, Y, U, V) :
    
    C = 0.001
    rho = 1000
    g = 9.81
    
    p = np.zeros((X.shape[0], X.shape[1]))
    ampl = 0
    
    for i in range(X.shape[0]) :
        for j in range(X.shape[1]) :
            if(U[i][j]==0 and V[i][j] == 0) :
                p[i][j] = np.nan
            else :
                ampl = (U[i][j])**2 + (V[i][j])**2 
                p[i][j] = rho*g*(C- ampl/(2*g))
    
    return p