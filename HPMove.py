import numpy as np

def HPMove(Si, Simin):
    dp = np.array([
        [-1, 1, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0],
        [0, 0, 0, 0, -1, 1]
    ])

    nd = dp.shape[1]

    # This loop picks a move at random until one is picked that is unequal
    # to the last two moves.
    d = -Simin

    while np.all(d == -Simin):
        idx = np.random.randint(nd)
        d = dp[:, idx]

    return d