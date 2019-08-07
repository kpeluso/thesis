from HELPERS_clean import TOL, np, corrEv, genStat, listMatch, normize
from brS_clean import brS

#
# MAYBE MAKE A PATCH TYPE WITH el(i,j) METHODS
#
# class Patch:
# 	def __init__(self, mat):
# 		self.mat = mat

def el(p, i, j):
	'''
	* Returns index of matrix element i,j according to patch p
	* i,j in {0,1}
	'''
	return (p[0][i], p[1][j])

def patchR(n):
	'''
	Generates random 2x2 patch
	'''
	output = np.zeros([2,2]) # [row_indices, col_indices]
	for i in range(2):
		output[i] = np.random.choice(range(n), size=[2], replace=False).astype(int)
	return output

def patch(n):
	'''
	Deterministic patch generator
	'''
	for i in range(n-1):
		for j in range(n-1):
			yield np.array([[i, i+1], [j, j+1]]).astype(int)

def AVG(n, pi, maxIters=float('inf'), tol=TOL):
	'''
	INPUT:
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
			# CAN'T HAVE ANY ZERO ENTRIES!
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via AVG Method
	'''
	L = lambda mat: np.linalg.norm(np.dot(mat, pi) - pi)
	M = lambda mat: np.linalg.norm(corrEv(mat) - pi)
	if n < 2:
		print "n must be >= 2"
	elif n == 2:
		output = brS(2, pi)
		return output, L(output), M(output)
	N = lambda x: normize(np.random.random([x,x]))
	output = N(n) # P^(I)
	p = patch(n)
	Z = (n-1)**2
	A = np.zeros([n,n])
	for _ in range((n-1)**2):
		cp = next(p) # current patch
		Ai = np.zeros([n,n])
		# Ai = np.identity(n)
		s = sum(pi[cp[0]])
		Ai[cp[0][0]:(cp[0][1]+1), cp[1][0]:(cp[1][1]+1)] = brS(2, pi[cp[0]]/s)
		A += Ai
	output = A/Z
	return output, L(output), M(output)

print '\nRUNNING AVG:\n'

errors = []
errors2 = []
b1 = 2
b2 = 11
for i in range(b1,b2):
	errors.append(AVG(i, genStat(i))[1])
	errors2.append(AVG(i, genStat(i))[2])
e = np.array(errors)
e2 = np.array(errors2)
print e
print e2
print sum(e)/ float(b2-b1)
print sum(e2)/ float(b2-b1)

# errors = []
# errors2 = []
# b = 10
# L = lambda mat, pi: np.linalg.norm(np.dot(mat, pi) - pi)
# M = lambda mat, pi: np.linalg.norm(corrEv(mat) - pi)
# for _ in range(b):
# 	pi = genStat(2)
# 	errors.append(L(brS(2, pi), pi))
# 	errors2.append(M(brS(2, pi), pi))
# e = np.array(errors)
# e2 = np.array(errors2)
# # print e
# print sum(e)/ float(b)
# print sum(e2)/ float(b)

n = 100
print AVG(n, genStat(n))[1]
print AVG(n, genStat(n))[2]

print '\nAVG DONE!\n'
