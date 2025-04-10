import numpy as np

# Load the 2-dom.txt file
dom_matrix = np.loadtxt('CL/2-dom.txt', dtype=int)

# Create the new cl1 matrix based on 2-dom.txt
cl1_matrix = np.zeros(dom_matrix.shape, dtype=float)

# Set nodes where dom_matrix is 2 to 0.00000000
cl1_matrix[dom_matrix == 2] = 0.00000000

# Calculate Q_out, Q_in, and q_spe
Q_out = (10 * 0 + 5 * 1) * 10**(-3)  
# Load the 2-dom.txt file
dom_matrix = np.loadtxt('CL/2-dom.txt', dtype=int)

# Create the new cl1 matrix based on 2-dom.txt
cl1_matrix = np.zeros(dom_matrix.shape, dtype=float)
cl1_matrix[dom_matrix == 2] = 0.00000000 


# Save the cl1_matrix to a new text file
with open('CL/2-cl.txt', 'w') as f:
    for row in cl1_matrix:
        # Format each row to match the structure of 1-cl.txt
        formatted_row = '   '.join(f'{val:.8f}' for val in row)
        f.write(f'{formatted_row}\n')
    for i in range (2,21):
        Q_out = (10*0 + 5*1)*10**(-3)
        Q_in = Q_out/2
        q_spe = Q_in/ 18

print("cl2.txt has been created successfully.")# Q_out = 5e-3
Q_in = Q_out / 2                     # Q_in = 2.5e-3
q_spe = Q_in / 18                    # q_spe = (2.5e-3) / 18 ≈ 0.00013888888888888889

# Set q_spe in the first row (row index 0) and columns 2 to 20 (column indices 1 to 19)
cl1_matrix[1, 79:99] = q_spe
cl1_matrix[99, 79:99] = q_spe

# Save the cl1_matrix to a new text file
with open('CL/2-cl.txt', 'w') as f:
    for row in cl1_matrix:
        # Format each row to match the structure of 1-cl.txt
        formatted_row = '   '.join(f'{val:.8f}' for val in row)
        f.write(f'{formatted_row}\n')

print("cl2.txt has been created successfully.")