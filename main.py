import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from fluid_dynamics import create_system, solve_system, calculate_pressure, velocity
from tools import circu, deriv, force
from tools.constante import rho, g, p_ref, ref_point, Q_out

def main():
    def canal_rectiligne():
        print("\nRésolution du canal rectiligne :\n")
        try:
            dom_1 = np.loadtxt('CL/1-dom.txt', dtype = int)
            cl_1 = np.loadtxt('CL/1-cl.txt', dtype = float)
            num_1 = np.loadtxt('CL/1-num.txt', dtype = int)
        except:                
            print("Erreur lors de l'ouverture des fichiers")
            return None

        A, B = create_system(dom_1, num_1, cl_1)
        print(f"A shape: {A.shape}, B shape: {B.shape}\n")
        if A.shape[0] != A.shape[1]:
            raise ValueError("Matrix A is not square!")
        if A.shape[0] != B.shape[0]:
            raise ValueError("A and B dimension mismatch!")
        if np.linalg.cond(A) > 1e12:
            print("Warning: Matrix A is ill-conditioned. Results may be inaccurate.")
        
        print(f"La matrice A vaut \n {A} \n")
        print(f"La matrice B vaut \n {B} \n")
        
        X = solve_system(A, B)
        print(f"La matrice X vaut \n {X} \n")

        if X.ndim == 1 and dom_1.shape:
            sol_grid = np.zeros_like(dom_1, dtype=float)
            for i in range(len(X)):
                sol_grid[np.where(num_1 == i + 1)] = X[i]  # Map solution back to domain
            
            plt.figure(figsize=(10, 8))
            plt.imshow(sol_grid, cmap='viridis')
            plt.colorbar()
            
            # Ajouter les numéros des noeuds
            for i in range(dom_1.shape[0]):
                for j in range(dom_1.shape[1]):
                    if num_1[i, j] != 0:  # Ne pas afficher les zéros
                        plt.text(j, i, str(num_1[i, j]), 
                                ha='center', va='center', 
                                color='white', fontweight='bold', fontsize=22)
            
            # Ajouter la grille
            plt.grid(True, which='major', color='black', linewidth=2)
            plt.xticks(np.arange(-.5, dom_1.shape[1], 1), [])
            plt.yticks(np.arange(-.5, dom_1.shape[0], 1), [])
            
            plt.title("Visualisation de la solution", fontsize=22)
            plt.show()

            # Calcul du champ de vitesse
            h = 0.5  # pas spatial pour le canal rectiligne (d'après le README)
            u, v = velocity(sol_grid, dom_1, h)
            
            # Visualisation du champ de vitesse
            plt.figure(figsize=(12, 8))
            plt.quiver(u, v)
            plt.title("Champ de vitesse")
            plt.show()

            # Création d'une grille de z (hauteur constante)
            z = np.zeros_like(sol_grid)
            for i in range(z.shape[0]):
                z[i, :] = 0
            
            pressure = calculate_pressure(u, v, z, rho, g, None, p_ref, ref_point)
            
            # Visualisation de la pression
            plt.figure(figsize=(12, 8))
            plt.imshow(pressure, cmap='RdBu')
            plt.colorbar(label='Pression (Pa)')
            plt.title("Distribution de pression")
            plt.show()

            # Calcul des forces sur les parois
            # Exemple pour une section verticale au milieu
            mid = sol_grid.shape[1] // 2
            x = np.array([i * h for i in range(sol_grid.shape[0])])
            y = np.zeros_like(x)
            p_section = pressure[:, mid]
            
            fx, fy = force(p_section, x, y)
            print("\nForces sur la section centrale:")
            print(f"Fx = {fx:.2f} N")
            print(f"Fy = {fy:.2f} N")
            
            # Calcul de la circulation
            # Exemple pour un contour rectangulaire au centre
            x_contour = [2*h, 2*h, 3*h, 3*h, 2*h]
            y_contour = [2*h, 3*h, 3*h, 2*h, 2*h]
            u_contour = [u[2,2], u[2,3], u[3,3], u[3,2], u[2,2]]
            v_contour = [v[2,2], v[2,3], v[3,3], v[3,2], v[2,2]]
            
            circulation = circu(u_contour, v_contour, x_contour, y_contour)
            print(f"\nCirculation autour du contour central: {circulation:.2f} m²/s")

        return X, u, v, pressure
    
    def canal_en_j():
        
        print("\n############################################################################\n")
        print("Résolution du canal en \"J\" \n")
        
        try:
            dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
            num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
            contour_obj_2 = np.loadtxt('CL/2-contourObj.txt', dtype=int)
        except:
            print("Erreur lors de l'ouverture des fichiers")
            return None
        
        h = 0.01  # pas spatial de 0.01 m

        # Définition des conditions aux limites
        cl_2 = np.zeros_like(dom_2, dtype=float)
        
        # Q_out
        for i in range(2, 21):
            cl_2[1, i] = Q_out/19
            
        # Q_in1 et Q_in2
        for i in range(80,99):
            cl_2[1, i] = (Q_out/2)/19
            cl_2[99, i] = (Q_out/2)/19
        
        A, B = create_system(dom_2, num_2, cl_2)
        print(f"A shape: {A.shape}, B shape: {B.shape}\n")
        if A.shape[0] != A.shape[1]:
            raise ValueError("Matrix A is not square!")
        if A.shape[0] != B.shape[0]:
            raise ValueError("A and B dimension mismatch!")
        if np.linalg.cond(A) > 1e12:
            print("Warning: Matrix A is ill-conditioned. Results may be inaccurate.")
        
        print(f"La matrice A vaut \n {A} \n")
        print(f"La matrice B vaut \n {B} \n")
        
        X = solve_system(A, B)
        print(f"La matrice X vaut \n {X} \n")

        if X.ndim == 1 and dom_2.shape:
            sol_grid = np.zeros_like(dom_2, dtype=float)
            for i in range(len(X)):
                sol_grid[np.where(num_2 == i + 1)] = X[i]  # Map solution back to domain
            
            plt.figure(figsize=(10, 8))
            plt.imshow(sol_grid, cmap='viridis')
            plt.colorbar()
            
            # Ajouter les numéros des noeuds
            for i in range(dom_2.shape[0]):
                for j in range(dom_2.shape[1]):
                    if num_2[i, j] != 0:  # Ne pas afficher les zéros
                        plt.text(j, i, str(num_2[i, j]), 
                                ha='center', va='center', 
                                color='white', fontweight='bold', fontsize=22)
            
            # Ajouter la grille
            plt.grid(True, which='major', color='black', linewidth=2)
            plt.xticks(np.arange(-.5, dom_2.shape[1], 1), [])
            plt.yticks(np.arange(-.5, dom_2.shape[0], 1), [])
            
            plt.title("Visualisation de la solution", fontsize=22)
            plt.show()

            # Calcul du champ de vitesse
            h = 0.5  # pas spatial pour le canal rectiligne (d'après le README)
            u, v = velocity(sol_grid, dom_2, h)
            
            # Visualisation du champ de vitesse
            plt.figure(figsize=(12, 8))
            plt.quiver(u, v)
            plt.title("Champ de vitesse")
            plt.show()

            # Création d'une grille de z (hauteur constante)
            z = np.zeros_like(sol_grid)
            for i in range(z.shape[0]):
                z[i, :] = i*h
            
            pressure = calculate_pressure(u, v, z, rho, g, None, p_ref, ref_point)
            
            # Visualisation de la pression
            plt.figure(figsize=(12, 8))
            plt.imshow(pressure, cmap='RdBu')
            plt.colorbar(label='Pression (Pa)')
            plt.title("Distribution de pression")
            plt.show()

            # Calcul des forces sur les parois
            # Exemple pour une section verticale au milieu
            mid = sol_grid.shape[1] // 2
            x = np.array([i * h for i in range(sol_grid.shape[0])])
            y = np.zeros_like(x)
            p_section = pressure[:, mid]
            
            fx, fy = force(p_section, x, y)
            print("\nForces sur la section centrale:")
            print(f"Fx = {fx:.2f} N")
            print(f"Fy = {fy:.2f} N")
            
            # Calcul de la circulation
            # Exemple pour un contour rectangulaire au centre
            x_contour = [2*h, 2*h, 3*h, 3*h, 2*h]
            y_contour = [2*h, 3*h, 3*h, 2*h, 2*h]
            u_contour = [u[2,2], u[2,3], u[3,3], u[3,2], u[2,2]]
            v_contour = [v[2,2], v[2,3], v[3,3], v[3,2], v[2,2]]
            
            circulation = circu(u_contour, v_contour, x_contour, y_contour)
            print(f"\nCirculation autour du contour central: {circulation:.2f} m²/s")

        return X, u, v, pressure

    X_rect, u_rect, v_rect, pressure_rect = canal_rectiligne()
<<<<<<< HEAD
    X_j, u_j, v_j, pressure_j = canal_en_j()
    dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
    num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
    contour_obj_2 = np.loadtxt('CL/2-contourObj.txt', dtype=int)
    cl_2 = np.loadtxt('CL/2-cl.txt', dtype = int )
    
    
    return X_rect, u_rect, v_rect, pressure_rect, X_j, u_j, v_j, pressure_j, dom_2, num_2, contour_obj_2, cl_2

if __name__ == "__main__":
    X_rect, u_rect, v_rect, pressure_rect, X_j, u_j, v_j, pressure_j, dom_2, num_2, contour_obj_2, cl_2 = main()
=======
    X_j, u_j, v_j, pressure_j, cl_2 = canal_en_j()
    
    return X_rect, u_rect, v_rect, pressure_rect, X_j, u_j, v_j, pressure_j, cl_2

if __name__ == "__main__":
    X_rect, u_rect, v_rect, pressure_rect, X_j, u_j, v_j, pressure_j, cl_2 = main()
    
    
    
    dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
    num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
    contour_obj_2 = np.loadtxt('CL/2-contourObj.txt', dtype=int)
>>>>>>> 92024a60e71d10a2e27a4836478c63feeae35812
