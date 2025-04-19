import numpy as np
from tools.constante import Q_out

def calculate_pressure(X, Y, U, V,dom) :
    """
    Calculate the pressure of a domain.

    Parameters
    ----------
    X: array 2-D 
        The value of x-axis where the pressure is wanted
    Y: array 2-D 
        The value of y-axis where the pressure is wanted
    U: array 2-D 
        The value of the horizontal velocity where the pressure is wanted
    V: array 2-D 
        The value of the vertical velocity where the pressure is wanted
    dom : array 2-D
        The domain where the pressure is wanted

    Returns
    -------
    p : array 2-D
        The pressure in the domain.
    """
    
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