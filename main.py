import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

import laplacian as lp

def main():
    
    try:
        dom_1 = np.loadtxt('CL/1-dom.txt', dtype = int)
        cl_1 = np.loadtxt('CL/1-cl.txt', dtype = float)
        num_1 = np.loadtxt('CL/1-num.txt', dtype = int)
    except:                
        print("Erreur lors de l'ouverture des fichiers")
        return None

    A, B = lp.create_system(dom_1, num_1, cl_1)
    print(f"A shape: {A.shape}, B shape: {B.shape}")
    if A.shape[0] != A.shape[1]:
        raise ValueError("Matrix A is not square!")
    if A.shape[0] != B.shape[0]:
        raise ValueError("A and B dimension mismatch!")
    if np.linalg.cond(A) > 1e12:
        print("Warning: Matrix A is ill-conditioned. Results may be inaccurate.")
    
    print(f"La matrice A vaut \n {A} \n")
    print(f"La matrice B vaut \n {B} \n")
    
    X = lp.solve_system(A, B)
    print(f"La matrice X vaut \n {X} \n")

    if X.ndim == 1 and dom_1.shape:
        sol_grid = np.zeros_like(dom_1, dtype=float)
        sol_grid[dom_1 == 1] = X  # Map solution back to domain
        plt.imshow(sol_grid, cmap='viridis')
        plt.colorbar()
        plt.title("Solution Visualization")
        plt.show()
    
    return X

if __name__ == "__main__":
    main()