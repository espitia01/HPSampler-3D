import numpy as np

def HPDistance(S, cache=None):
    if cache is None:
        cache = {}

    # Check if the result is already cached
    hash_key = tuple(S.flatten())
    if hash_key in cache:
        return cache[hash_key]

    # S is a 3D array with the movements in the x-, y-, and z-directions, respectively
    positions = np.cumsum(S, axis=1)
    N = S.shape[1]
    
    # Create a 3D array of pairwise position differences
    diff = positions[:, :, np.newaxis] - positions[:, np.newaxis, :]
    
    # Compute the squared Euclidean distances
    squared_distances = np.sum(diff**2, axis=0)
    
    # Set the diagonal elements to a large value (e.g., 100) to avoid self-interactions
    np.fill_diagonal(squared_distances, 100)
    
    # Take the square root of the squared distances
    distances = np.sqrt(squared_distances)
    
    # Cache the result
    cache[hash_key] = distances
    
    return distances