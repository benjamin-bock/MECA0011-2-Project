import numpy as np

# Load the domain matrix
dom_matrix = np.loadtxt('2-dom.txt', dtype=int)

# Create the cl1 matrix based on the domain
cl1_matrix = np.zeros(dom_matrix.shape, dtype=float)

# Set nodes where dom_matrix is 2 to 0.00000000 (Dirichlet condition)
cl1_matrix[dom_matrix == 2] = 0.00000000

# Calculate Q_out based on group number (example: group 58)
X = 0   # X is the tens digit of the group number
Y = 1  # Y is the units digit of the group number
Q_out = (10 * X + 5 * Y) * 10**(-3)  # Q_out in m³/s
Q_in = Q_out / 2  # Q_in = Q_out / 2
q_spe = Q_in  # Uniform distribution over 18 nodes
cl1_matrix[dom_matrix == 2] = Q_out
# Apply Dirichlet conditions at the inlet and outlet
# Example: Set q_spe in the first row (row index 0) and columns 2 to 20
cl1_matrix[1, 1:21] = Q_out
cl1_matrix[1, 79:100] = q_spe  # Inlet
cl1_matrix[99, 79:100] = q_spe  # Outlet

# Apply Neumann conditions around the obstacle
# Example: Set derivative conditions around the obstacle
# This requires identifying the nodes around the obstacle and applying the derivative

# Save the cl1_matrix to a new text file
with open('2-cl.txt', 'w') as f:
    for row in cl1_matrix:
        formatted_row = '   '.join(f'{val:.8f}' for val in row)
        f.write(f'{formatted_row}\n')

print("Boundary conditions have been set successfully.")