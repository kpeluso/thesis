import numpy as np
from VSEA_clean import stochTest
import cma

N = 3; n = 2 # people, bins
pi = [2.0/3.0, 1.0/3.0] # ideal stationary distribution
tol = 0.01
fitnessFunction = None
sigma0 = .25**n # ~1/4 of search domain width => try with .25 and (.25)**n
# build x0
x0 = np.zeros(n**2)
for i in range(n):
	nxt = np.random.uniform(0,1,n)
	x0[i*n:i*n+n] = nxt/sum(nxt)


def devectorize(vm, bins):
	'''
	INPUT:
		vm :: NPArray<Float>
			# bins^2-length vectorized P^I
		bins :: Integer
			# the number of bins
	OUTPUT:
		NPArray<NPArray<Float>>
			# devectorized vm as bins-x-bins-matrix
	'''
	output = np.zeros([bins, bins])
	col = 0
	for i in xrange(0, len(vm), bins):
		output[:,col] = vm[i:i+bins]
		col += 1
	return output


def anyNeg(l):
	for el in l:
		if el < 0.0:
			return True
	return False


def corrEig(vm, corr, bins, tol):
	'''
	INPUT:
		vm :: NPArray<Float>
			# bins^2-length vectorized P^I
		corr :: NPArray<Float>
			# bins-length correct stationary distribution
		bins :: Integer
			# the number of bins
		tol :: Float
			# tolerance
	OUTPUT:
		Float
			# The error between corr and the first eigenvector of vm
			# If first eigenvalue is not within tol of 1.0, then np.NaN is returned
	'''
	mat = devectorize(vm, bins)
	out = np.linalg.eig(mat)
	eigenvalue1 = out[0][0]
	# if first eigenvalue is not 1.0 or columns not probability distributions (between 0-1 and sum to 1)
	if not abs(eigenvalue1 - 1.0) < tol or anyNeg(vm) or not stochTest(mat, tol=tol):
		return np.NaN
	else:
		eigenvector1 = out[1][:,0]
		return sum(abs(corr - eigenvector1)) # 1-D loss
		#return corr - eigenvector1 # [bins]-D loss


f = cma.fitness_functions.FitnessFunctions()
preFit = lambda x: corrEig(x, pi, n, tol)
fitnessFunction = lambda x: f.fun_as_arg(x, preFit)
# result = cma.fmin(fitnessFunction, x0, sigma0)
# print result
def CMAES(fitFun, start_x, start_sig):
	return cma.fmin(fitFun, start_x, start_sig)

