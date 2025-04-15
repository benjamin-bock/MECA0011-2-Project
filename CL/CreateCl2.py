import numpy as np

# Load the domain matrix
dom_matrix = np.loadtxt('2-dom.txt', dtype=int)
contour = np.loadtxt('2-contourObj.txt', dtype=int)

# Create the cl1 matrix based on the domain
CL = np.zeros(dom_matrix.shape, dtype=float)

X = 0   # X is the tens digit of the group number
Y = 1  # Y is the units digit of the group number
h = 0.01
Q_out = (10 * X + 5 * Y) * 10**(-3)  
Q_in = Q_out / 2 

L = 21*h
u_intlet = Q_in/L
u_outlet = Q_out/L

# 2 entrées du canal 
for i in range(0,21):
    CL[1][79 + i] = -u_intlet * h * i  
    CL[99][79 + i] = u_intlet * h * i

#2 canal d'entré interieu 
CL[1:41,79] = 0.000
CL[60:100,79] = 0.00

#canal d'entré exterieur
CL[1:41,99] = -u_intlet*h*21
CL[41:100,99] = u_intlet*h*21

#sortie du canal
for i in range(0,20): 
    CL[1][i+1] = u_outlet*h*i
    
CL[1:61, 1] = 0.0000

# Apply fully developed flow conditions
for index in contour:
    if index[0] >=50: 
        CL[tuple(index)] = Q_in/2
    if index[0]  < 50: 
        if index[1] >= 32: 
            CL[tuple(index)] = - Q_in/2
        if index[1] < 32 : 
            CL[tuple(index)] = Q_in/2
    
CL[1:41, 21] = u_outlet*h*21

CL[40,21:79] = u_intlet*h*21
CL[11:40,86] = -u_intlet*h*21

CL[39:42,85] = -u_intlet*h*21

CL[41:43,84] = -u_intlet*h*21
CL[42:44,83] = -u_intlet*h*21

CL[43:45,82] = -u_intlet*h*21






# Save the boundary conditions to a file
with open('2-cl.txt', 'w') as f:
    for row in CL: 
        formatted_row = '   '.join(f'{val:.8f}' for val in row)
        f.write(f'{formatted_row}\n')

print("Boundary conditions have been set successfully.")