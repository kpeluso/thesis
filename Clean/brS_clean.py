from HELPER-MODULE import *

def brS(N, n, pi, iters=20000, tol=TOL):
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
			# P^(I) via Binary-Random Search Method
			# only works in 2x2 case
			# BECAUSE IT'S 2x2 CASE I CAN USE 2x2 P^(M) SERIES FORMULA
	'''
	output = np.ones([2,2])*.5 # P^(I)
	ev = corrEv(output)
	b1 = 0.0 if pi[0] >= pi[1] else 1.0-(pi[0]/pi[1])
	b2 = 1.0
	output[1][1] = np.average([b1,b2])
	output = resMat(output)
	while not listMatch(np.dot(output, ev), pi) and iters > 0: # s1, loop
		# look for random number within updated b1, b2
			if ev[1] < pi[1]:
				b1 = np.average([b1,b2])
				output[1][1] += next_term(b1,b2)
			else:
				b2 = np.average([b1,b2])
				output[1][1] -= next_term(b1,b2)
			output[0][0] = p22p1(output[1][1], pi)
			output = resMat(output)
			ev = corrEv(output)
			iters -= 1
	return output
