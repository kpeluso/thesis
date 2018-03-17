from ind2mass import *
from math import pi as PI

def GRS(N, n, pi, maxIters=n*10000, tol=TOL):
	'''
	INPUT:
		N :: Integer
			# population size
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Greedy Random Search Method
	'''
	# initialization and normalization
	if maxIters == None:
		maxIters = n*10000
	output = np.random.random((n,n))
	for col in xrange(n):
		output[:,col] /= sum(output[:,col])
	out = np.linalg.eig(output)
	while not (abs(out[0][0] - 1.0) < tol and listMatch(out[1][:,0], pi)) and maxIters != 0: # Step 2
		a = np.random.randint(0, high=n) # Step 3
		v = np.random.random((n,)) # Step 4
		c = 2*int(out[1][a,0] > pi[a]) - 1 # Step 5
		output[a,:] = np.random.dirichlet(output[a,:] + c*v) # Step 6
		for col in xrange(n):
			output[:,col] /= sum(output[:,col]) # Step 7
		out = np.linalg.eig(output)
		maxIters -= 1
	return output

