from HELPERS_clean import TOL, np, listMatch

def NRS(n, pi, iterCols=float('inf'), iterInCol=10000, tol=TOL):
	'''
	INPUT:
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Naive Random Search Method
	'''
	# initialization and normalization
	if iterCols == None:
		iterCols = n*10000
	output = np.random.random((n,n))
	for col in xrange(n):
		output[:,col] /= sum(output[:,col])
	out = np.linalg.eig(output)
	n0s = np.zeros(n)
	totalIts = iterCols; totalItsIn = iterInCol
	while (not abs(out[0][0] - 1.0) < tol or not listMatch(out[1][:,0], pi)) and iterCols!=0:
		for col in xrange(n):
			p_old = n0s
			iterInCol = totalItsIn
			temp = output
			while not listMatch(output[:,col], p_old) and iterInCol!=0:
				next_try = np.random.random((n,))
				next_try /= sum(next_try) # Line 5 in LaTeX
				temp[:,col] = next_try
				out = np.linalg.eig(temp)
				if np.linalg.norm(out[1][:,0] - pi)**2 < np.linalg.norm(p_old - pi)**2:
					output = temp
					p_old = output[:,col]
				iterInCol -= 1
		out = np.linalg.eig(output)
		iterCols -= 1
	return output

