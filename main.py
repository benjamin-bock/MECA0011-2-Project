import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib.colors import LinearSegmentedColormap

from fluid_dynamics import create_system, solve_system, calculate_pressure, velocity
from tools import circu, deriv, force
from tools.constante import rho, g, p_ref, ref_point, Q_out
from CL.createcl2 import createCL2

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
        
        print(f"La matrice A vaut \n {A} \n")
        print(f"La matrice B vaut \n {B} \n")
        
        X = solve_system(A, B)
        print(f"La matrice X vaut \n {X} \n")

        if X.ndim == 1 and dom_1.shape:
            sol_grid = np.zeros_like(dom_1, dtype=float)
            for i in range(len(X)):
                sol_grid[np.where(num_1 == i + 1)] = X[i]  # Map solution back to domain
            
            plt.figure(figsize=(14, 12))
            img = plt.imshow(sol_grid, cmap='coolwarm', origin='upper', interpolation='nearest')
            cbar = plt.colorbar(img, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=12) 
            cbar.set_label('Valeur de la solution', fontsize=14)
            
            # Ajouter les numéros des noeuds
            for i in range(dom_1.shape[0]):
                for j in range(dom_1.shape[1]):
                    if num_1[i, j] != 0:  # Ne pas afficher les zéros
                        plt.text(j, i, str(num_1[i, j]), 
                                ha='center', va='center', 
                                color='white', fontweight='bold', fontsize=20)
            
            # Ajouter la grille
            plt.grid(True, which='both', linestyle='-', linewidth=1, color='black', alpha=0.3)
            plt.xticks(np.arange(-0.5, dom_1.shape[1], 1), [])
            plt.yticks(np.arange(-0.5, dom_1.shape[0], 1), [])
            plt.gca().set_facecolor('#f0f0f0')
            
            plt.title("Visualisation de la solution", fontsize=20, pad = 20)
            
            for spine in plt.gca().spines.values():
               spine.set_color('dimgray')
               spine.set_linewidth(1.5)
               
            plt.tight_layout(pad = 3)
            plt.show()

            h = 0.5  # pas spatial pour le canal rectiligne (d'après le README)
            u, v = velocity(sol_grid, dom_1, h)
            
            speed = np.sqrt(u**2 + v**2)
            # Visualisation du champ de vitesse
            plt.style.use('seaborn-v0_8-whitegrid') 
            fig, ax = plt.subplots(figsize=(14, 10), dpi=100)
            q = ax.quiver(u, v, speed, 
              cmap='inferno',  # Palette de couleurs
              scale=20,         # Échelle des flèches
              width=0.002,      # Largeur des flèches
              headwidth=5,      # Largeur des têtes de flèches
              alpha=0.8)        # Transparence

            # Ajout d'une barre de couleur
            cbar = fig.colorbar(q, ax=ax, label='Vitesse (m/s)', shrink=0.8)
            cbar.ax.tick_params(labelsize=12)

             # Personnalisation de l'affichage 
            ax.set_title(r'Champ de vitesse dans le canal ($h = 0.01$)', 
             fontsize=16, pad=20, fontfamily='serif')
            ax.set_xlabel('Position x (m)', fontsize=14, fontfamily='serif')
            ax.set_ylabel('Position y (m)', fontsize=14, fontfamily='serif')
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_aspect('equal')  # Conserve les proportions
            plt.tight_layout()
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
            cl_2 = np.loadtxt("CL/2-cl.txt",dtype = float)
        except:
            print("Erreur lors de l'ouverture des fichiers")
            return None
        
        h = 0.01  #m

        # Définition des conditions aux limites
        cl_2 = createCL2(dom_2, num_2, contour_obj_2)
        
        
        A, B = create_system(dom_2, num_2, cl_2)
        print(f"A shape: {A.shape}, B shape: {B.shape}\n")

        print(f"La matrice A vaut \n {A} \n")
        print(f"La matrice B vaut \n {B} \n")
        
        X = solve_system(A, B)
        print(f"La matrice X vaut \n {X} \n")
        
        colors = [(1, 1, 1)]  # Blanc pour 0
        colors.extend(plt.cm.viridis(np.linspace(0, 1, 255))) 
        custom_cmap = LinearSegmentedColormap.from_list('custom', colors)
        white_inferno = LinearSegmentedColormap.from_list(
        'white_inferno',
        [(0.0, 'white')] + [(i / 255.0, plt.cm.Reds(i)) for i in range(1, 256)])
        

        if X.ndim == 1 and dom_2.shape:
            sol_grid = np.zeros_like(dom_2, dtype=float)
            for i in range(len(X)):
                sol_grid[np.where(num_2 == i + 1)] = X[i]  # Map solution back to domain
           
            
           
            plt.figure(figsize=(14, 12))
            img = plt.imshow(sol_grid, cmap="PuOr", interpolation='nearest')
            cbar = plt.colorbar(img, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=12) 
            cbar.set_label('Valeur de la solution', fontsize=14)
           
            # Ajouter la grille
            plt.grid(True, which='both', linestyle='-', linewidth=0.5, color='black', alpha=0.3)
            plt.xticks(np.arange(-0.5, dom_2.shape[1], 1), [])
            plt.yticks(np.arange(-0.5, dom_2.shape[0], 1), [])
            plt.gca().set_facecolor('#f0f0f0')
           
            plt.title("Visualisation de la solution", fontsize=20, pad = 20)
           
            for spine in plt.gca().spines.values():
               spine.set_color('dimgray')
               spine.set_linewidth(1.5)
              
            plt.tight_layout(pad = 3)
            plt.show()
            

            # Calcul du champ de vitesse
            h = 0.01  # pas spatial pour le canal rectiligne (d'après le README)
            u, v = velocity(sol_grid, dom_2, h)
            
            speed = np.sqrt(u**2 + v**2)
            
            # Visualisation du champ de vitesse
            plt.style.use('seaborn-v0_8-whitegrid')
            fig, ax = plt.subplots(figsize=(14, 10), dpi=100)
            x, y = np.meshgrid(np.arange(u.shape[1]), np.arange(u.shape[0]))
            x, y = np.meshgrid(np.arange(u.shape[1]), np.arange(u.shape[0]))
            y = np.flipud(y)  # Inverse l'ordre des indices de y

            
            # Utilisation correcte de quiver
            q = ax.quiver(x, y, u, v, speed,
                          cmap='inferno',  # Palette de couleurs
                          scale=5,         # Échelle des flèches
                          width=0.001,      # Largeur des flèches
                          headwidth=5,      # Largeur des têtes de flèches
                          alpha=0.8)        # Transparence
            
            # Ajout d'une barre de couleur
            cbar = fig.colorbar(q, ax=ax, label='Vitesse (m/s)', shrink=0.8)
            cbar.ax.tick_params(labelsize=12)
            
            # Personnalisation de l'affichage
            ax.set_title(r'Champ de vitesse dans le canal ($h = 0.01$)',
                         fontsize=16, pad=20, fontfamily='serif')
            ax.set_xlabel('Position x (m)', fontsize=14, fontfamily='serif')
            ax.set_ylabel('Position y (m)', fontsize=14, fontfamily='serif')
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_aspect('equal')  # Conserve les proportions
            
            # Ajout de lignes de courant (optionnel)
            ax.streamplot(x, y, u, v, color='w', linewidth=0.5)
            
            plt.tight_layout()
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

        return X, u, v, pressure, cl_2

    #X_rect = canal_rectiligne()
    X_j  = canal_en_j()
    dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
    num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
    contour_obj_2 = np.loadtxt('CL/2-contourObj.txt', dtype=int)
    cl_2 = np.loadtxt('CL/2-cl.txt', dtype = float)
    
    
    return  X_j, dom_2, num_2, contour_obj_2, cl_2

if __name__ == "__main__":
    X_j, dom_2, num_2, contour_obj_2, cl_2 = main()
