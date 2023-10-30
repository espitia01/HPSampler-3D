import numpy as np

def HPMove(Si, Simin):
	
	dp = np.concatenate(([[-1,1,0,0]],[[0,0,-1,1]]))
	nd = np.size(dp,1)

	#This loop picks a move at random until one is picked that is unequal
	#to the last two moves.
	d = -Simin
	while all(d==-Simin): #or all(d==Si):
		d = dp[:,np.random.randint(nd)]

	#print(d, Si, -Simin)
	return d

#S = np.array([[0,-1,0,0,0,1,-1,1,0,1,0],[0,0,1,1,-1,0,0,0,-1,0,1]])
#P = np.array([0,1,0,1,1,0,0,0,1,0,0,1])
#print(HPEnergy(P,S,1))
