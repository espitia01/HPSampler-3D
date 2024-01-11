import numpy as np

def HPDistance(S):
    # S is a 3D array with the movements in the x-, y-, and z-directions, respectively
    positions = np.cumsum(S, axis=1)
    N = S.shape[1]
    Distances = np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            if i == j:
                Distances[i, j] = 100
            else:
                Distances[i, j] = np.sqrt(np.sum((positions[:, i] - positions[:, j])**2))
    
    return Distances