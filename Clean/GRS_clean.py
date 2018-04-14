from HELPERS_clean import np, TOL, listMatch, corrEv

def GRS(n, pi, maxIters=float('inf'), tol=TOL):
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
		maxIters = n*100000
	output = np.random.random((n,n))
	for col in xrange(n):
		output[:,col] /= sum(output[:,col])
	out = np.linalg.eig(output)
	ev = corrEv(out[1])
	while not listMatch(np.dot(output, pi), pi) and maxIters != 0: # Step 2
		a = np.random.randint(0, high=n) # Step 3
		c = 0.45*float((2*int(ev[a] < pi[a]) - 1)) # Step 4
		# Step 5
		output[a,a] = output[a,a]*(int((1.0+c)*output[a,a] > 1.0)+int(abs((1.0+c)*output[a,a]) < tol)) \
			+ (1.0+c)*output[a,a]*int((1.0+c)*output[a,a] <= 1.0)*int((1.0+c)*output[a,a] > 0.0)
		output[:,a] /= sum(output[:,a]) # Step 6
		out = np.linalg.eig(output) # preparation for Step 2
		ev = corrEv(out[1]) # preparation for Step 2
		maxIters -= 1 # Step 7
	return output # Step 8

