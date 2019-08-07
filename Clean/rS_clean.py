from HELPERS_clean import TOL, np, resMat, corrEv, p22p1, genStat, listMatch

def rS(n, pi, iters=float('inf'), tol=TOL):
	'''
	INPUT:
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Random Search Method
			# only works in 2x2 case
			# BECAUSE IT'S 2x2 CASE I CAN USE 2x2 P^(M) SERIES FORMULA
	'''
	output = np.zeros([2,2]) # P^(I)
	for col in xrange(2):
		output[:,col] = np.transpose(genStat(2))
	ev = corrEv(output)
	b1 = 0.0 if pi[0] >= pi[1] else 1.0-(pi[0]/pi[1])
	b2 = 1.0
	output[1][1] = np.average([b1,b2])
	output = resMat(output)
	while not listMatch(np.dot(output, ev), pi, tol=tol) and iters > 0: # s1, loop
		output[1][1] = (b2-b1)*np.random.random()+b1
		output[0][0] = p22p1(output[1][1], pi) # calculate p_11
		output = resMat(output) # calculate p_12, p_21
		ev = corrEv(output)
		iters -= 1
	return output

