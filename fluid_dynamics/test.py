import time
import numpy as np
import matplotlib.pyplot as plt

# Import the functions to test
from laplacian import create_system, min_exclude_zero, solve_system
from laplacian import create_sparse_system, solve_sparse_system

def generate_test_data(rows, cols):
    """
    Generates test data for the system creation functions.
    
    Parameters
    ----------
    rows : int
        Number of rows in the grid.
    cols : int
        Number of columns in the grid.
        
    Returns
    -------
    dom_1 : matrix
        Domain matrix where all nodes are interior points (1).
    num_1 : matrix
        Node numbering matrix.
    cl_1 : matrix
        Boundary conditions matrix (all zeros for this test).
    """
    num_1 = np.arange(1, rows * cols + 1).reshape(rows, cols)
    dom_1 = np.ones_like(num_1)
    cl_1 = np.zeros_like(dom_1, dtype=float)
    
    # Set some boundary conditions for testing
    dom_1[0, :] = 2  # Top boundary
    dom_1[-1, :] = 2  # Bottom boundary
    dom_1[:, 0] = 2  # Left boundary
    dom_1[:, -1] = 2  # Right boundary
    
    # Set Dirichlet conditions on boundaries
    cl_1[0, :] = 1.0  # Top boundary
    cl_1[-1, :] = 0.0  # Bottom boundary
    cl_1[:, 0] = 0.0  # Left boundary
    cl_1[:, -1] = 1.0  # Right boundary
    
    return dom_1, num_1, cl_1

def compare_systems(A_dense, B_dense, A_sparse, B_sparse):
    """
    Compares the dense and sparse systems to ensure they are equivalent.
    
    Parameters
    ----------
    A_dense : array
        Dense coefficient matrix.
    B_dense : array
        Dense right-hand side vector.
    A_sparse : csr_matrix
        Sparse coefficient matrix.
    B_sparse : array
        Sparse right-hand side vector.
        
    Returns
    -------
    bool
        True if systems are equivalent, False otherwise.
    """
    # Compare right-hand side vectors
    if not np.allclose(B_dense, B_sparse):
        print("Right-hand side vectors differ!")
        return False
    
    # Compare coefficient matrices
    A_sparse_dense = A_sparse.toarray()
    if not np.allclose(A_dense, A_sparse_dense):
        print("Coefficient matrices differ!")
        return False
    
    print("Systems are equivalent.")
    return True

def main():
    # Grid sizes to test
    grid_sizes = [10, 20, 50, 100]
    
    non_sparse_times = []
    sparse_times = []
    
    for size in grid_sizes:
        print(f"\nTesting grid size: {size}x{size}")
        
        # Generate test data
        dom_1, num_1, cl_1 = generate_test_data(size, size)
        
        # Measure non-sparse function time
        start_time = time.time()
        A_dense, B_dense = create_system(dom_1, num_1, cl_1)
        non_sparse_time = time.time() - start_time
        non_sparse_times.append(non_sparse_time)
        print(f"Non-sparse time: {non_sparse_time:.6f} seconds")
        
        # Measure sparse function time
        start_time = time.time()
        A_sparse, B_sparse = create_sparse_system(dom_1, num_1, cl_1)
        sparse_time = time.time() - start_time
        sparse_times.append(sparse_time)
        print(f"Sparse time: {sparse_time:.6f} seconds")
        
        # Compare systems
        compare_systems(A_dense, B_dense, A_sparse, B_sparse)
        
        # Solve systems and compare solutions
        X_dense = solve_system(A_dense, B_dense)
        X_sparse = solve_sparse_system(A_sparse, B_sparse)
        
        if np.allclose(X_dense, X_sparse):
            print("Solutions are equivalent.")
        else:
            print("Solutions differ!")
    
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(grid_sizes, non_sparse_times, label='Non-Sparse Function', marker='o')
    plt.plot(grid_sizes, sparse_times, label='Sparse Function', marker='x')
    plt.xlabel('Grid Size')
    plt.ylabel('Time (seconds)')
    plt.title('Time Comparison Between Non-Sparse and Sparse Functions')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()