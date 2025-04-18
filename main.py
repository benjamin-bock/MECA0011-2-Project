import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches

from fluid_dynamics import create_system, solve_system, calculate_pressure, velocity
from tools import circu 
from tools.constante import rho, g, p_ref, ref_point, Q_out
from CL.createcl2 import createCL2, createCL2q3

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
            sol_grid = np.full_like(num_1, np.nan, dtype=float)
            sol_grid[num_1 > 0] = X[num_1[num_1 > 0] - 1]
            
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
            
            # Préparation des grilles pour quiver
            ny, nx = sol_grid.shape
            x = np.arange(nx) * h
            y = np.arange(ny) * h
            Xg, Yg = np.meshgrid(x, y)
            Xr = Yg.T
            Yr = Xg.T
            u_rot = v.T
            v_rot = u.T

            plt.figure(figsize=(10, 7))
            plt.quiver(Xr, Yr, u_rot, v_rot, color='cornflowerblue')
            plt.xlabel("y [m]")
            plt.ylabel("x [m]")
            plt.title("Champ de vitesse")
            plt.axis("equal")
            plt.grid(True)
            plt.show()

            
            pressure = calculate_pressure(Xr,Yr,u_rot,v_rot)
            
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
            
            # Calcul de la circulation
            # Exemple pour un contour rectangulaire au centre
            x_contour = [2*h, 2*h, 3*h, 3*h, 2*h]
            y_contour = [2*h, 3*h, 3*h, 2*h, 2*h]
            u_contour = [u[2,2], u[2,3], u[3,3], u[3,2], u[2,2]]
            v_contour = [v[2,2], v[2,3], v[3,3], v[3,2], v[2,2]]
            
            circulation = circu(u_contour, v_contour, x_contour, y_contour)
            print(f"\nCirculation autour du contour central: {circulation:.2f} m²/s")

        return X, u, v, pressure
    
    def canal_en_j(plot_solution=False, plot_velocity=False, plot_pressure=False):
        print("\n############################################################################\n")
        print("Résolution du canal en \"J\" \n")
        
        try:
            dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
            num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
            contour_obj_2 = np.loadtxt('CL/2-contourObj.txt', dtype=int)
        except:
            print("Erreur lors de l'ouverture des fichiers")
            return None
        
        h = 0.01  # Pas spatial

        # Définition des conditions aux limites
        cl_2 = createCL2(dom_2, num_2, contour_obj_2)
        cl_2q3 = createCL2q3(dom_2, num_2, contour_obj_2)
        A, B = create_system(dom_2, num_2, cl_2)
        A3, B3 = create_system(dom_2, num_2, cl_2q3)

        print(f"A shape: {A.shape}, B shape: {B.shape}\n")

        psi = solve_system(A, B)
        psi3 = solve_system(A3, B3)
        print(f"La matrice $\psi$ vaut \n {psi} \n")
        blue = (0x93/255, 0x0F/255, 0xFF/255)   # Bleu foncé
        gray = (0xAB/255, 0x9F/255, 0x9F/255) # Gris plus foncé
        yellow = (0xFF/255, 0x73/255, 0x00/255) 
        cdict = {
        'red': [(0.0, blue[0], blue[0]),
            (0.5, gray[0], gray[0]),
            (1.0, yellow[0], yellow[0])],
        'green': [(0.0, blue[1], blue[1]),
              (0.5, gray[1], gray[1]),
              (1.0, yellow[1], yellow[1])],
        'blue': [(0.0, blue[2],blue[2]),
             (0.5, gray[2], gray[2]),
             (1.0, yellow[2], yellow[2])]
        }
        purple_gray_orange = LinearSegmentedColormap('purple_gray_orange', cdict)
        if psi.ndim == 1 and dom_2.shape:
            sol_grid = np.full_like(num_2, np.nan, dtype=float)
            sol_grid[num_2 > 0] = psi[num_2[num_2 > 0] - 1]

            # Affichage de la solution
            if plot_solution:
                plt.figure(figsize=(14, 12))
                img = plt.imshow(sol_grid, cmap="PuOr", interpolation='nearest')
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=12) 
                cbar.set_label('Valeur de la solution', fontsize=14)
                plt.title("Visualisation de la solution", fontsize=20, pad=20)
                plt.grid(True, which='both', linestyle='-', linewidth=0.5, color='black', alpha=0.3)
                plt.tight_layout(pad=3)
                plt.show()
            
            u, v = velocity(sol_grid, dom_2, h)
            u_filtered = np.full_like(u, np.nan, dtype=float)
            v_filtered = np.full_like(v, np.nan, dtype=float)
            u_filtered[dom_2 == 1] = u[dom_2 == 1]
            v_filtered[dom_2 == 1] = v[dom_2 == 1]
            ny, nx = sol_grid.shape
            x = np.arange(nx) * h
            y = np.arange(ny) * h
            Xg, Yg = np.meshgrid(x, y)
            Xr = Yg.T
            Yr = Xg.T
            u_rot = v_filtered.T[::-1]
            v_rot = u_filtered.T[::-1]
            # Calcul et affichage du champ de vitesse
            if plot_velocity:
                

                plt.figure(figsize=(10, 7))
                img = plt.imshow(u_rot, cmap=purple_gray_orange, interpolation='bilinear', extent=[0,1,0,1])
                plt.colorbar(img, fraction=0.046, pad=0.04)
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseHorizontale.pdf', dpi = 300, format = 'pdf')
                plt.show()


                plt.figure(figsize=(10, 7))
                img = plt.imshow(v_rot, cmap=purple_gray_orange, interpolation='bilinear', extent=[0,1,0,1])
                plt.colorbar(img, fraction=0.046, pad=0.04)
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseVerticale.pdf', dpi = 300, format = 'pdf')
                plt.show()

            # Calcul et affichage de la pression
            if plot_pressure:
                pressure = calculate_pressure(Xr, Yr, v.T, u.T, dom_2.T)
                pressure = pressure[::-1]

                fig, ax = plt.subplots(figsize=(12, 8))
                im = ax.imshow(pressure, cmap='plasma', interpolation='bilinear', extent= [0,1,0,1])
                plt.colorbar(im, label='Pression (Pa)')
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.savefig('Figures/Pressure.pdf', dpi = 300, format = 'pdf')
                plt.show()

        if psi3.ndim == 1 and dom_2.shape:
            sol_grid3 = np.full_like(num_2, np.nan, dtype=float)
            sol_grid3[num_2 > 0] = psi3[num_2[num_2 > 0] - 1]

            if plot_solution:
                plt.figure(figsize=(14, 12))
                img = plt.imshow(sol_grid3, cmap="PuOr", interpolation='nearest', extent= [0,1,0,1])
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=12) 
                cbar.set_label('Valeur de la solution', fontsize=14)
                plt.title("Visualisation de la solution", fontsize=20, pad=20)
                plt.grid(True, which='both', linestyle='-', linewidth=0.5, color='black', alpha=0.3)
                plt.tight_layout(pad=3)
                plt.xlabel("x")
                plt.ylabel("y")
                plt.show()
            
            u3, v3 = velocity(sol_grid3, dom_2, h)
            u3_filtered = np.full_like(u3, np.nan, dtype=float)
            v3_filtered = np.full_like(v3, np.nan, dtype=float)
            u3_filtered[dom_2 == 1] = u3[dom_2 == 1]
            v3_filtered[dom_2 == 1] = v3[dom_2 == 1]
            ny, nx = sol_grid3.shape
            x = np.arange(nx) * h
            y = np.arange(ny) * h
            Xg, Yg = np.meshgrid(x, y)
            Xr3 = Yg.T
            Yr3 = Xg.T
            u_rot3 = v3_filtered.T[::-1]
            v_rot3 = u3_filtered.T[::-1]


            # Calcul et affichage du champ de vitesse
            if plot_velocity:

                plt.figure(figsize=(10, 7))
                img = plt.imshow(u_rot3, cmap=purple_gray_orange, interpolation='bilinear', extent=[0,1,0,1])
                plt.colorbar(img, fraction=0.046, pad=0.04)
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseHorizontaleQ3.pdf', dpi = 300, format = 'pdf')
                plt.show()


                plt.figure(figsize=(10, 7))
                img = plt.imshow(v_rot3, cmap=purple_gray_orange, interpolation='bilinear', extent=[0,1,0,1])
                plt.colorbar(img, fraction=0.046, pad=0.04)
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseVerticaleQ3.pdf', dpi = 300, format = 'pdf')
                plt.show()

            # Calcul et affichage de la pression
            if plot_pressure:
                pressure = calculate_pressure(Xr3, Yr3, v3.T,u3.T, dom_2.T)
                pressure = pressure[::-1]

                fig, ax = plt.subplots(figsize=(12, 8))
                im = ax.imshow(pressure, cmap='plasma', interpolation='bilinear',extent= [0,1,0,1])
                plt.colorbar(im, label='Pression (Pa)')
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.savefig('Figures/PressureQ3.pdf', dpi = 300, format = 'pdf')
                plt.show()
            return psi

    #X_rect = canal_rectiligne()
    X_j  = canal_en_j(plot_pressure=True)
    dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
    num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
    contour_obj_2 = np.loadtxt('CL/2-contourObj.txt', dtype=int)
    cl2 = createCL2(dom_2,num_2,contour_obj_2)

    return  X_j, dom_2, num_2, contour_obj_2

if __name__ == "__main__":
    X_j, dom_2, num_2, contour_obj_2 = main()
