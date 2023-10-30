import numpy as np
import matplotlib.pyplot as plt
from HPDistance import HPDistance

def HPEnergy(P, S, J):
    N = P.size
    D = HPDistance(S)

    # Set diagonal and sub-diagonal to a value other than 1
    for i in range(N):
        D[i, i] = 0
        if i < N - 1:
            D[i, i + 1] = 0
            D[i + 1, i] = 0

    # Create matrix with the product of the protein array
    P2 = np.outer(P, P)  # This replaces the loop and np.ones

    # Calculate energy
    E = np.sum((D == 1) * P2)

    E *= 0.5 * J

    return E

