import numpy as np

def calculate_pressure(u, v, z, rho, g, C, p_ref, ref_point):

    velocity_sq = u**2 + v**2
    
    # Détermination de la constante C
    if C is None:
        if p_ref is None or ref_point is None:
            raise ValueError("Either provide C or (p_ref and ref_point)")
        i, j = ref_point
        C = (p_ref/(rho*g)) + z[i,j] + velocity_sq[i,j]/(2*g)
    
    # Calcul de la pression selon l'équation de Bernoulli
    pressure = rho * g * (C + z - velocity_sq/(2*g))
    
    return pressure