import matplotlib.pyplot as plt
import numpy as np
from HPShow import HPShow

S = np.array([[ 0,  1,  0, -1],
 [ 0,  0,  1,  0],
 [ 0,  0,  0,  0]])

P = np.array([1, 0, 1, 1])

HPShow(P, S)