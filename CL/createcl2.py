from tools.constante import Q_out
import numpy as np


def createCL2(dom, num, contour_obj):
    """
    Create the CL2 matrix based on the provided domain, number, and contour object matrices.
    
    Parameters:
    dom (numpy.ndarray): The domain matrix.
    num (numpy.ndarray): The number matrix.
    contour_obj (numpy.ndarray): The contour object matrix.
    
    Returns:
    numpy.ndarray: The created CL2 matrix.
    """
    # Initialize the CL2 matrix with zeros
    cl2 = np.zeros_like(dom, dtype=float)
    
    # Set values in the CL2 matrix based on the domain and contour object matrices
    
    h = 0.01
    ouverture = 21
    vitesse = Q_out  / (h * ouverture)
    
    # Q_out
    for i in range(1, 22):
        cl2[1, i] = vitesse * (i * h) + 1
        
    # Borne du domaine
    for i in range(2, 41, 1): # Horizontale supérieure
        cl2[i, 21] = cl2[1, 21]
    for i in range(21, 80, 1): # Montée gauche
        cl2[40, i] = cl2[1, 21]
    for i in range(1, 80): # Montée droite
        cl2[60, i] = cl2[1, 1]
    for i in range(2, 61): # Horizontale inférieure
        cl2[i, 1] = cl2[1, 1] 
        
        
        
    # Q_in
    for i in range(21):
        cl2[1, 79 + i] = cl2[1, 21] - vitesse * (i * h) / 2 # Ajoute + 1 pour voir le contour dans Spyder
        cl2[99, 79 + i] = cl2[1, 1] + vitesse * (i * h) / 2 
        
    # Borne supérieure du domaine
    for i in range(99):
        cl2[1+i, 99] = cl2[1, 99]
    
    # Borne inférieure de la barre du "J"
    for i in range(1, 41, 1):
        cl2[i, 79] = cl2[1, 21]
        
    for i in range(60, 100, 1):
        cl2[i, 79] = cl2[1, 1]
        
    # Contour obstacle
    for i in range(np.shape(contour_obj)[0]):
        cl2[contour_obj[i][0]][contour_obj[i][1]] = cl2[1, 99] # Moitié de Q_in
    
    
    
    
        
    np.savetxt("conditions_limites.txt", cl2, fmt="%.6f")
    return cl2


