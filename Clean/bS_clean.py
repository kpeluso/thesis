from HELPERS_clean import TOL, np, listMatch, resMat, corrEv, p22p1

def bS(n, pi, iters=float('inf'), tol=TOL):
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
			# P^(I) via Binary Search Method
			# only works in 2x2 case
			# BECAUSE IT'S 2x2 CASE I CAN USE 2x2 P^(M) SERIES FORMULA
	'''
	output = np.ones([2,2])*0.5 # P^(I)
	ev = corrEv(output)
	b1 = 0.0 if pi[0] >= pi[1] else 1.0-(pi[0]/pi[1])
	b2 = 1.0
	output[1][1] = np.average([b1,b2])
	output = resMat(output)
	while not listMatch(np.dot(output, ev), pi, tol=tol) and iters > 0: # s1, loop
		if ev[1] < pi[1]:
			b1 = np.average([b1,b2])
			output[1][1] += np.average([b1,b2])-b1
		else:
			b2 = np.average([b1,b2])
			output[1][1] -= np.average([b1,b2])-b1
		output[0][0] = p22p1(output[1][1], pi)
		output = resMat(output)
		ev = corrEv(output)
		iters -= 1
	return output

