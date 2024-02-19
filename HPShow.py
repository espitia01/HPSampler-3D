import numpy as np
import matplotlib.pyplot as plt

def HPShow(P,S,Energy=0.0,Temperature=0.0,Time=0.0):
    Path = np.cumsum(S, axis=1)
    Path = (Path.T-np.mean(Path,1)).T #center the protein
    N = len(P)
    xy = np.sqrt(2*N)/1.0
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot(Path[0,:], Path[1,:], Path[2,:],'-k') #the backbone
    ax.set_xlim([-xy, +xy])
    ax.set_ylim([-xy, +xy])
    ax.set_zlim([-xy, +xy])
    
    for i in range(N):
        if P[i]==0:
            ax.scatter(Path[0,i], Path[1,i], Path[2,i], c='b', s=12) #H(0) aa in blue
        elif P[i]==1:
            ax.scatter(Path[0,i], Path[1,i], Path[2,i], c='r', s=12) #P(1) aa in red
        ax.scatter(Path[0,0], Path[1,0], Path[2,0], c='g', s=4) #mark 1st aa in green
        
    ax.text(xy*(0.05-1),xy*0.9,xy*0.8,"T = %.2f"%Temperature)
    ax.text(xy*(0.05-1),xy*0.8,xy*0.7," t = %.0f"%Time)
    ax.text(xy*(0.05-1),xy*0.7,xy*0.6,"E = %.2f"%Energy)
    
    plt.show(block=False)
    plt.pause(10)
    plt.close(fig)  # Close the figure
    return