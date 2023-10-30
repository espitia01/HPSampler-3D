import numpy as np

def HPDistance(S):
#S is a 2D array with the movements in the x- and y-directions, respectively
    positions = np.cumsum(S, axis=1)
    N = S.shape[1]
    Distances = np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            if i == j:
                Distances[i, j] = 100
            else:
                Distances[i, j] = np.sum((positions[:, i] - positions[:, j])**2)
    
    return Distances