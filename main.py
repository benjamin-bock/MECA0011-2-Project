import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches

from fluid_dynamics import create_system, solve_system, calculate_pressure, velocity
from tools import circu 
from tools.constante import rho, g, p_ref, ref_point, Q_out
from CL.createcl2 import createCL2, createCL2q3, createCL2q5

# Définir les tailles de police par défaut pour tous les graphiques
plt.rcParams.update({
    'font.size': 16,            # Augmenté de 14 à 16
    'axes.labelsize': 20,       # Augmenté de 16 à 20
    'axes.titlesize': 24,       # Augmenté de 18 à 24
    'xtick.labelsize': 16,      # Augmenté de 14 à 16
    'ytick.labelsize': 16,      # Augmenté de 14 à 16
    'legend.fontsize': 16,      # Augmenté de 14 à 16
    'axes.titlepad': 25         # Augmenté de 20 à 25
})



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
            
            plt.figure(figsize=(16, 14))  # Augmenté de (14, 12)
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

            plt.figure(figsize=(12, 9))  # Augmenté de (10, 7)
            plt.quiver(Xr, Yr, u_rot, v_rot, color='cornflowerblue')
            plt.xlabel("y [m]")
            plt.ylabel("x [m]")
            plt.title("Champ de vitesse")
            plt.axis("equal")
            plt.grid(True)
            plt.show()

            
            pressure = calculate_pressure(Xr,Yr,u_rot,v_rot)
            
            # Visualisation de la pression
            plt.figure(figsize=(14, 10))  # Augmenté de (12, 8)
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
    
    def canal_en_j(plot_solution=False, plot_velocity=False, plot_pressure=False,plot_streamlines=False):
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
        cl_2q5 = createCL2q5(dom_2,num_2, contour_obj_2)
        A, B = create_system(dom_2, num_2, cl_2)
        A3, B3 = create_system(dom_2, num_2, cl_2q3)
        A5, B5 = create_system(dom_2, num_2, cl_2q5)

        print(f"A shape: {A.shape}, B shape: {B.shape}\n")

        psi = solve_system(A, B)
        psi3 = solve_system(A3, B3)
        psi5 = solve_system(A5, B5)
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
                plt.figure(figsize=(14, 10))  # Dans canal_en_j(), augmenté de (12, 8)
                img = plt.imshow(sol_grid, cmap="PuOr", interpolation='nearest')
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=14) 
                cbar.set_label('Valeur de la solution', fontsize=16)
                plt.grid(True, which='both', linestyle='-', linewidth=0.5, color='black', alpha=0.3)
                plt.tight_layout(pad=3)
                plt.show()
            
            u, v = velocity(sol_grid, dom_2, h)
            u_filtered = np.full_like(u, np.nan, dtype=float)
            v_filtered = np.full_like(v, np.nan, dtype=float)
            u_filtered[dom_2 > 0] = u[dom_2 >0]
            v_filtered[dom_2 > 0] = v[dom_2 > 0]
            ny, nx = sol_grid.shape
            x = np.arange(nx) * h
            y = np.arange(ny) * h
            Xg, Yg = np.meshgrid(x, y)
            Xr = Yg.T
            Yr = Xg.T
            u_rot = v_filtered.T[::-1]
            v_rot = u_filtered.T[::-1]
            

            if plot_velocity:
                

                plt.figure(figsize=(13, 10))  # Augmenté de (10, 7)
                img = plt.imshow(u_rot, cmap=purple_gray_orange, extent=[0,1,0,1])
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04, ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse horizontale [m/s]', fontsize=16)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseHorizontale.pdf', dpi=300, format='pdf')
                plt.show()


                plt.figure(figsize=(12, 9))  # Augmenté de (10, 7)
                img = plt.imshow(v_rot, cmap=purple_gray_orange, extent=[0,1,0,1])
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04, ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse verticale [m/s]', fontsize=16)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseVerticale.pdf', dpi=300, format='pdf')
                plt.show()

            if plot_streamlines:
                
                plt.figure(figsize=(12, 9))  # Augmenté de (10, 7)
                speed = np.sqrt(u_filtered**2 + v_filtered**2)
                strm = plt.streamplot(Xg, Yg, u_filtered, v_filtered, 
                                    color=speed, cmap='plasma', 
                                    density=4, linewidth=0.5, 
                                    arrowsize=0.0)
                
                cbar = plt.colorbar(strm.lines,ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse [m/s]', fontsize=16, rotation=90, labelpad=15)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/Streamlines.pdf', dpi=300, format='pdf')
                plt.show()
                
                plt.figure(figsize=(10, 8))
                # Calculer les coordonnées physiques pour le zoom
                x_zoom = Xg[30:70, 10:50]  # Zone plus large
                y_zoom = Yg[30:70, 10:50]  # Zone plus large
                u_zoom = u_filtered[30:70, 10:50]  # Zone plus large
                v_zoom = v_filtered[30:70, 10:50]  # Zone plus large
                speed_zoom = np.sqrt(u_zoom**2 + v_zoom**2)

                strm_zoom = plt.streamplot(x_zoom, y_zoom, u_zoom, v_zoom,
                                         color=speed_zoom, cmap='plasma',
                                         density=1, linewidth=0.5,
                                         arrowsize=0.0, start_points= None,broken_streamlines= False)
                
                cbar = plt.colorbar(strm_zoom.lines)
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse [m/s]', fontsize=16, rotation=90, labelpad=15)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                plt.grid(True)
                plt.savefig('Figures/StreamlinesZoom.pdf', dpi=300, format='pdf')
                plt.show()


            # Calcul et affichage de la pression
            if plot_pressure:
                pressure = calculate_pressure(Xr, Yr, v.T, u.T, dom_2.T)
                pressure = pressure[::-1]

                fig, ax = plt.subplots(figsize=(14, 10))  # Augmenté de (12, 8)
                im = ax.imshow(pressure, cmap='jet', interpolation='bilinear', extent= [0,1,0,1])
                cbar = plt.colorbar(im,ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Pression [Pa]', fontsize=16)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                ax.tick_params(axis='both', which='major', labelsize=14)
                plt.grid(True)
                plt.savefig('Figures/Pressure.pdf', dpi=300, format='pdf')
                plt.show()
            
            # Calcul de la circulation autour de l'obstacle en T
            if len(contour_obj_2) > 0:
                # Extraire les coordonnées du contour
                x_contour = []
                y_contour = []
                u_contour = []
                v_contour = []
                
                # Pour chaque point du contour
                for i in range(contour_obj_2.shape[0]):
                    # Coordonnées du point dans la grille
                    row, col = contour_obj_2[i]
                    
                    # Convertir en coordonnées réelles
                    x_contour.append(col * h)
                    y_contour.append(row * h)
                    
                    # Obtenir les composantes de vitesse à ce point
                    # Vérifier que le point est dans le domaine fluide
                    if 0 <= row < u.shape[0] and 0 <= col < u.shape[1]:
                        u_contour.append(u[row, col])
                        v_contour.append(v[row, col])
                    else:
                        u_contour.append(0)
                        v_contour.append(0)
                
                # Fermer le contour en ajoutant le premier point à la fin
                if len(x_contour) > 0:
                    x_contour.append(x_contour[0])
                    y_contour.append(y_contour[0])
                    u_contour.append(u_contour[0])
                    v_contour.append(v_contour[0])
                
                # Calculer la circulation
                circulation = circu(u_contour, v_contour, x_contour, y_contour)
                print(f"\nCirculation autour de l'obstacle en T: {circulation:.6f} m²/s")
                
                # Visualiser le contour si demandé
                if plot_solution or plot_velocity:
                    plt.figure(figsize=(10, 7))
                    plt.plot(x_contour, y_contour, 'r-', linewidth=2)
                    plt.xlabel("x [m]")
                    plt.ylabel("y [m]")
                    plt.axis("equal")
                    plt.grid(True)
                    plt.savefig('Figures/Contour_Obstacle.pdf', dpi=300, format='pdf')
                    plt.show()

        if psi3.ndim == 1 and dom_2.shape:
            sol_grid3 = np.full_like(num_2, np.nan, dtype=float)
            sol_grid3[num_2 > 0] = psi3[num_2[num_2 > 0] - 1]

            if plot_solution:
                plt.figure(figsize=(14, 10))  # Augmenté de (12, 8)
                img = plt.imshow(sol_grid3, cmap="PuOr", interpolation='nearest', extent= [0,1,0,1])
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=16) 
                cbar.set_label('Valeur de la solution', fontsize=16)
                plt.grid(True, which='both', linestyle='-', linewidth=0.5, color='black', alpha=0.3)
                plt.tight_layout(pad=3)
                plt.xlabel("x")
                plt.ylabel("y")
                plt.show()
            
            u3, v3 = velocity(sol_grid3, dom_2, h)
            u3_filtered = np.full_like(u3, np.nan, dtype=float)
            v3_filtered = np.full_like(v3, np.nan, dtype=float)
            u3_filtered[dom_2 >0] = u3[dom_2 >0]
            v3_filtered[dom_2 >0] = v3[dom_2 >0]
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

                plt.figure(figsize=(13, 10))  # Augmenté de (10, 7)
                img = plt.imshow(u_rot3, cmap=purple_gray_orange, interpolation='bilinear', extent=[0,1,0,1])
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04,ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse horizontale [m/s]', fontsize=16)
                plt.xlabel("x [m]", fontsize = 16)
                plt.ylabel("y [m]", fontsize = 16)
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseHorizontaleQ3.pdf', dpi = 300, format = 'pdf')
                plt.show()


                plt.figure(figsize=(12, 9))  # Augmenté de (10, 7)
                img = plt.imshow(v_rot3, cmap=purple_gray_orange, interpolation='bilinear', extent=[0,1,0,1])
                cbar = plt.colorbar(img, fraction=0.046, pad=0.04,ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse horizontale [m/s]', fontsize=16)
                plt.xlabel("x [m]", fontsize = 16)
                plt.ylabel("y [m]", fontsize = 16)
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/VitesseVerticaleQ3.pdf', dpi = 300, format = 'pdf')
                plt.show()

            # Calcul et affichage de la pression
            if plot_pressure:
                pressure = calculate_pressure(Xr3, Yr3, v3.T,u3.T, dom_2.T)
                pressure = pressure[::-1]

                fig, ax = plt.subplots(figsize=(14, 10))  # Augmenté de (12, 8)
                im = ax.imshow(pressure, cmap='jet', interpolation='bilinear',extent= [0,1,0,1])
                cbar = plt.colorbar(im,ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=14)
                cbar.set_label('Pression [Pa]', fontsize=16)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                plt.title("Distribution de pression", fontsize=18, pad=15)
                ax.tick_params(axis='both', which='major', labelsize=14)
                plt.grid(True)
                plt.savefig('Figures/PressureQ3.pdf', dpi = 300, format = 'pdf')
                plt.show()
                
            # Calcul de la circulation autour de l'obstacle en T pour la solution Q3
            if len(contour_obj_2) > 0:
                # Extraire les coordonnées du contour
                x_contour3 = []
                y_contour3 = []
                u_contour3 = []
                v_contour3 = []
                
                # Pour chaque point du contour
                for i in range(contour_obj_2.shape[0]):
                    # Coordonnées du point dans la grille
                    row, col = contour_obj_2[i]
                    
                    # Convertir en coordonnées réelles
                    x_contour3.append(col * h)
                    y_contour3.append(row * h)
                    
                    # Obtenir les composantes de vitesse à ce point
                    # Vérifier que le point est dans le domaine fluide
                    if 0 <= row < u3.shape[0] and 0 <= col < u3.shape[1]:
                        u_contour3.append(u3[row, col])
                        v_contour3.append(v3[row, col])
                    else:
                        u_contour3.append(0)
                        v_contour3.append(0)
                
                # Fermer le contour en ajoutant le premier point à la fin
                if len(x_contour3) > 0:
                    x_contour3.append(x_contour3[0])
                    y_contour3.append(y_contour3[0])
                    u_contour3.append(u_contour3[0])
                    v_contour3.append(v_contour3[0])
                
                # Calculer la circulation
                circulation3 = circu(u_contour3, v_contour3, x_contour3, y_contour3)
                print(f"\nCirculation autour de l'obstacle en T (Q3): {circulation3:.6f} m²/s")
                
                # Visualiser le contour si demandé
                if plot_solution or plot_velocity:
                    plt.figure(figsize=(10, 7))
                    plt.plot(x_contour3, y_contour3, 'r-', linewidth=2)
                    plt.xlabel("x [m]")
                    plt.ylabel("y [m]")
                    plt.title("Contour de l'obstacle (Q3)")
                    plt.axis("equal")
                    plt.grid(True)
                    plt.savefig('Figures/Contour_Obstacle_Q3.pdf', dpi=300, format='pdf')
                    plt.show()
        if psi5.ndim == 1 and dom_2.shape:
            sol_grid5 = np.full_like(num_2, np.nan, dtype=float)
            sol_grid5[num_2 > 0] = psi5[num_2[num_2 > 0] - 1]
            u5, v5 = velocity(sol_grid5, dom_2, h)
            u5_filtered = np.full_like(u5, np.nan, dtype=float)
            v5_filtered = np.full_like(v5, np.nan, dtype=float)
            u5_filtered[dom_2 > 0] = u5[dom_2 > 0]
            v5_filtered[dom_2 > 0] = v5[dom_2 > 0]
            ny, nx = sol_grid.shape
            x = np.arange(nx) * h
            y = np.arange(ny) * h
            Xg, Yg = np.meshgrid(x, y)
            Xr = Yg.T
            Yr = Xg.T
            u5_rot = v_filtered.T[::-1]
            v5_rot = u_filtered.T[::-1]

            if plot_streamlines:
                """
                plt.figure(figsize=(12, 9))  # Augmenté de (10, 7)
                speed = np.sqrt(u5_filtered**2 + v5_filtered**2)
                strm = plt.streamplot(Xg, Yg, u5_filtered, v5_filtered, 
                                    color=speed, cmap='plasma', 
                                    density=3, linewidth=0.5, 
                                    arrowsize=0.0)
                
                cbar = plt.colorbar(strm.lines,ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse [m/s]', fontsize=16, rotation=90, labelpad=15)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                plt.axis("equal")
                plt.grid(True)
                plt.savefig('Figures/StreamlinesCorrect.pdf', dpi=300, format='pdf')
                plt.show()
                """

                plt.figure(figsize=(10, 8))
                # Calculer les coordonnées physiques pour le zoom
                x_zoom = Xg[30:70, 10:50]  # Zone plus large
                y_zoom = Yg[30:70, 10:50]  # Zone plus large
                u_zoom = u5_filtered[30:70, 10:50]  # Zone plus large
                v_zoom = v5_filtered[30:70, 10:50]  # Zone plus large
                speed_zoom = np.sqrt(u_zoom**2 + v_zoom**2)

                strm_zoom = plt.streamplot(x_zoom, y_zoom, u_zoom, v_zoom,
                                         color=speed_zoom, cmap='plasma',
                                         density=1, linewidth=0.5,
                                         arrowsize=0.0, start_points= None,broken_streamlines= False)
                
                cbar = plt.colorbar(strm_zoom.lines)
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label('Vitesse [m/s]', fontsize=16, rotation=90, labelpad=15)
                plt.xlabel("x [m]", fontsize=16)
                plt.ylabel("y [m]", fontsize=16)
                plt.grid(True)
                plt.savefig('Figures/StreamlinesZoomQ5.pdf', dpi=300, format='pdf')
                plt.show()

                    
            return psi
    def plot_initial_conditions(cl_2):
        dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
        num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
    
    # Créer la figure avec la même taille que l'exemple
        plt.figure(figsize=(10, 8))

    
    # Définir les couleurs personnalisées
        
        cl_2[dom_2 == 2]  = cl_2[dom_2 == 2]   
        cl_2[dom_2 == 1] = np.nan 
        cl_2[dom_2 == 0] = np.nan
        cl_2 = cl_2.T[::-1]
    # Afficher avec imshow
        plt.imshow(cl_2, cmap='coolwarm', extent=[0, 1, 0, 1])
    
    # Ajouter une barre de couleur
        cbar = plt.colorbar(ticks=[-1, 0, 1])
        cbar.set_label('Valeurs des conditions limite [m²/s]')
    
    # Paramètres du graphique
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.grid(True)
    
        plt.savefig('Figures/ConditionsInitiales.pdf', dpi=300, format='pdf', bbox_inches='tight')
        plt.show()

    #X_rect = canal_rectiligne()
    X_j  = canal_en_j(plot_velocity=True)
    dom_2 = np.loadtxt('CL/2-dom.txt', dtype=int)
    num_2 = np.loadtxt('CL/2-num.txt', dtype=int)
    contour_obj_2 = np.loadtxt('CL/2-contourObj.txt', dtype=int)
    cl2 = createCL2q5(dom_2,num_2,contour_obj_2)

    return  X_j, dom_2, num_2, contour_obj_2

if __name__ == "__main__":  # Ajouter cet appel avant les autres
    X_j, dom_2, num_2, contour_obj_2 = main()
