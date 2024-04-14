import numpy as np
import matplotlib.pyplot as plt

def HPShow(P, S, Energy=0.0, Temperature=0.0, Time=0.0):
    Path = np.cumsum(S, axis=1)
    Path = (Path.T - np.mean(Path, 1)).T  # center the protein

    N = len(P)
    xy = np.sqrt(2 * N) / 1.0

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the backbone
    ax.plot(Path[0, :], Path[1, :], Path[2, :], '-k')

    # Set axis limits
    ax.set_xlim([-xy, +xy])
    ax.set_ylim([-xy, +xy])
    ax.set_zlim([-xy, +xy])

    # Scatter plot for amino acids
    colors = np.where(P == 0, 'b', 'r')
    ax.scatter(Path[0, :], Path[1, :], Path[2, :], c=colors, s=12)

    # Mark the first amino acid in green
    ax.scatter(Path[0, 0], Path[1, 0], Path[2, 0], c='g', s=4)

    # Add text annotations
    ax.text(xy * (0.05 - 1), xy * 0.9, xy * 0.8, f"T = {Temperature:.2f}")
    ax.text(xy * (0.05 - 1), xy * 0.8, xy * 0.7, f" t = {Time:.0f}")
    ax.text(xy * (0.05 - 1), xy * 0.7, xy * 0.6, f"E = {Energy:.2f}")

    plt.show(block=False)
    plt.pause(10)
    plt.close(fig)  # Close the figure

    return