import numpy as np
import matplotlib.pyplot as plt

def HPShow(P,S,Energy=0.0,Temperature=0.0,Time=0.0):
	
	Path = np.cumsum(S, axis=1)
	Path = (Path.T-np.mean(Path,1)).T #center the protein
	N = len(P)
	xy = np.sqrt(2*N)/1.0
	
	plt.clf()
	plt.plot(Path[0,:], Path[1,:],'-k') #the backbone
	plt.axis([-xy, +xy, -xy, +xy])
	
	for i in range(N):
		if P[i]==0:
			plt.plot(Path[0,i], Path[1,i],'bo',markersize=12) #H(0) aa in blue
		elif P[i]==1:
			plt.plot(Path[0,i], Path[1,i],'ro',markersize=12) #P(1) aa in red
		plt.plot(Path[0,0], Path[1,0],'go',markersize=4) #mark 1st aa in green
		
	plt.text(xy*(0.05-1),xy*0.9,"T = %.2f"%Temperature)
	plt.text(xy*(0.05-1),xy*0.8," t = %.0f"%Time)
	plt.text(xy*(0.05-1),xy*0.7,"E = %.2f"%Energy)
	
	plt.show(block=False)
	plt.pause(0.1)
	return