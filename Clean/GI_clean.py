from HELPERS_clean import TOL, np, corrEv, genStat, listMatch
from brS_clean import brS

def GI(n, pi, iters=float('inf'), tol=TOL):
	'''
	INPUT:
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
			# CAN'T HAVE ANY ZERO ENTRIES!
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Gibbs Sampling-inspired Method
	'''
	output = np.zeros([n,n]) # P^(I)
	for col in xrange(n):
		output[:,col] = np.transpose(genStat(n))
	ev = corrEv(output)
	indices = range(n)
	while not listMatch(np.dot(output, ev), pi) and iters > 0: # s1, loop
		# s2, isolate
		alter = np.random.choice(indices, size=[2], replace=False).astype(int)
		alter = np.array([min(alter), max(alter)]) # sort in order of lowest to highest
		subpi = np.zeros(2)
		subpi[0] = pi[alter[0]]
		subpi[1] = pi[alter[1]]
		# s3b, note how much space was formerly taken up
		resMass_mat = (output[alter[0]][alter[0]] + output[alter[1]][alter[0]], \
			output[alter[0]][alter[1]] + output[alter[1]][alter[1]])
		resMass_pi = sum(subpi)
		# s3, normalize
		subpi /= sum(subpi)
		# s4, optimize extracted 2-equation system
		submat = brS(n, subpi) # !!! Use bS, rS, brS methods. !!!
		# s5a, denormalize
		submat[:,0] *= resMass_mat[0]
		submat[:,1] *= resMass_mat[1]
		subpi *= resMass_pi
		# s5, substitute in new values renormalized to Q
		output[alter[0]][alter[0]] = submat[0][0]
		output[alter[1]][alter[0]] = submat[1][0]
		output[alter[0]][alter[1]] = submat[0][1]
		output[alter[1]][alter[1]] = submat[1][1]
		ev = corrEv(output)
		iters -= 1
	return output

