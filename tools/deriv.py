import numpy as np
import scipy as sp

def deriv(f_left, f_c, f_right, type_left, type_c, type_right, h) :
    if type_c == 0:
        raise ValueError("Invalid value for type_cent: 0 is not in the domain of calcul.")
    
    # Initialiser v à None pour détecter si aucun cas n'est traité
    v = None
    
    if type_c == 1:
        # Différence centrée pour les points intérieurs
        v = (f_right - f_left) / (2*h)
    elif type_c == 2:
        if type_left == 0:
            # Différence avant si point gauche hors domaine
            v = (f_right - f_c) / h
        elif type_right == 0:
            # Différence arrière si point droit hors domaine
            v = (f_c - f_left) / h
        else:
            # Différence centrée par défaut
            v = (f_right - f_left) / (2*h)
    
    # Vérifier que v a été définie
    if v is None:
        raise ValueError(f"No valid case found for types: left={type_left}, center={type_c}, right={type_right}")
    
    return v