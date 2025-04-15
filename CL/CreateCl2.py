import numpy as np

# Load the domain matrix
dom_matrix = np.loadtxt('2-dom.txt', dtype=int)
contour = np.loadtxt('2-contourObj.txt', dtype = int)
# Create the cl1 matrix based on the domain
cl1_matrix = np.zeros(dom_matrix.shape, dtype=float)


cl1_matrix[dom_matrix == 2] = 0.00000000
cl1_matrix[dom_matrix ==  1] = 0.00000000


X = 0   # X is the tens digit of the group number
Y = 1  # Y is the units digit of the group number
h = 0.01
Q_out = (10 * X + 5 * Y) * 10**(-3)  
Q_in = Q_out / 2 

y = np.arange(0,60 *h,60 )

cl1_matrix[1, 1:22] = -Q_out
cl1_matrix[1, 79:100] = Q_in
cl1_matrix[99, 79:100] = Q_in

for index in contour:
    cl1_matrix[tuple(index)] = Q_in
    


with open('2-cl.txt', 'w') as f:
    for row in cl1_matrix:
        formatted_row = '   '.join(f'{val:.8f}' for val in row)
        f.write(f'{formatted_row}\n')

print("Boundary conditions have been set successfully.")
print(contour)