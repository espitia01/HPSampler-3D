import numpy as np
import matplotlib.pyplot as plt
from HPDistance import HPDistance

def HPEnergy(P, S, J, D_cache=None):
    N = len(P)
    
    # Compute the distance matrix using the optimized HPDistance function
    D = HPDistance(S, cache=D_cache)
    
    # Set the diagonal and sub-diagonal elements to 0
    np.fill_diagonal(D, 0)
    np.fill_diagonal(D[1:], 0)
    np.fill_diagonal(D[:, 1:], 0)
    
    # Create a mask for elements equal to 1
    mask = (D == 1)
    
    # Create a matrix with the product of the protein array
    P2 = np.outer(P, P)
    
    # Calculate energy using the mask and protein matrix
    E = np.sum(mask * P2) * 0.5 * J
    
    return E