from ind2mass import *
from math import pi as PI

def NRS(N, n, pi, iterCols=None, iterInCol=10000, tol=TOL):
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
		#print 'iterCols', iterCols
		for col in xrange(n):
			p_old = n0s
			iterInCol = totalItsIn
			temp = output
			while not listMatch(output[:,col], p_old) and iterInCol!=0:
				#print 'iterInCol', iterInCol
				# testPt = truncnorm.pdf(np.random.random(), 0.0, 1.0)
				# next_try = truncnorm.pdf(testPt, 0.0, 1.0) # Line 5 in LaTeX
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


algo = False
if algo:
	print '\nNRS Tests:'
	g = NRS
else:
	print '\nGRS Tests:'
	g = GRS

testpi = np.array([12, 1])/13.0
N = 3
n = 2
rho = 0.4 # correlation
print '\n'
ans = g(N, n, testpi)
print 'eigenvalue:', np.linalg.eig(ans)[0][0]
ev = np.linalg.eig(ans)[1][:,0]
print 'eigenvector:', ev/sum(ev)
#print np.linalg.eig(ans)
print 'pi', testpi

testpi = np.array([2, 10, 1])/13.0
N = 13
n = 3
rho = 0.4 # correlation
print '\n'
ans = g(N, n, testpi)
print 'eigenvalue:', np.linalg.eig(ans)[0][0]
ev = np.linalg.eig(ans)[1][:,0]
print 'eigenvector:', ev/sum(ev)
#print np.linalg.eig(ans)
print 'pi', testpi

testpi = np.array([2, 7, 0, 3, 1])/13.0
N = 13
n = 5
rho = 0.4 # correlation
print '\n'
ans = g(N, n, testpi)
print 'eigenvalue:', np.linalg.eig(ans)[0][0]
ev = np.linalg.eig(ans)[1][:,0]
print 'eigenvector:', ev/sum(ev)
#print np.linalg.eig(ans)
print 'pi', testpi

