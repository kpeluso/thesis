import numpy as np
from scipy.special import comb

TOL = 0.0001

p22p1 = lambda p2,pi: 1.0 - (1.0 - p2)*(pi[1]/pi[0]) # proof of formula given in rough draft of thesis

vect = lambda mat: [mat[0][0], mat[1][0], mat[0][1], mat[1][1]]

next_term = lambda a,b: (b-a)*np.random.random() # b >= a

wc_count = lambda n,k: int(comb(n+k-1,k-1)) # count total number of wcs, n=balls, k=boxes


def resMat(mat):
	mat[1][0] = 1.0-mat[0][0]
	mat[0][1] = 1.0-mat[1][1]
	return mat


def normize(mat):
	for i in xrange(len(mat)):
		mat[:,i] /= sum(mat[:,i])
	return mat


def corrEv(mat, tol=TOL):
	'''
	INPUT:
		mat :: List<List<Float>>
			# 2x2
	OUTPUT:
		NPArray<Float>
			# returns eigenvector corresponding to eigenvalue of 1
	'''
	info = np.linalg.eig(mat)
	ev = info[1][:,0] if abs(info[0][0] - 1.0) < tol else info[1][:,1]
	return ev/sum(ev)


def corrEv_m(mat, tol=TOL):
	'''
	INPUT:
		mat :: List<List<Float>>
			# square
	OUTPUT:
		NPArray<Float>
			# returns eigenvector corresponding to eigenvalue of 1
	'''
	info = np.linalg.eig(mat)
	for i in xrange(len(info[1][:,0])):
		if abs(info[0][i] - 1.0) < tol:
			ev = info[1][:,i]
			return ev/sum(ev)


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


#
# Source: (2nd answer in link)
#  https://stackoverflow.com/questions/4647120/next-composition-of-n-into-k-parts-does-anyone-have-a-working-algorithm
#
def weak_compositions(balls, boxes, parent=tuple()):
	if boxes > 1:
		for i in xrange(balls + 1):
			for x in weak_compositions(i, boxes - 1, parent + (balls - i,)):
				yield x
	else:
		yield parent + (balls,)


#
# Source: LeetCode
#
def permute(nums):
	perms = [[]]
	for n in nums:
		new_perms = []
		for perm in perms:
			for i in range(len(perm)+1):
				new_perms.append(perm[:i] + [n] + perm[i:])
		perms = new_perms
	return perms


def perm(n, r):
	return factorial(n)/factorial(n-r)


def genStat(n):
	'''
	INPUT:
		n :: Integer
			# size of stationary distribution
	OUTPUT:
		List<Float>
			# generate a random stationary distribution
	'''
	samp = np.random.uniform(0,1,n)
	return samp/sum(samp)


def listMatch(l1, l2, tol=TOL):
	'''
	INPUT:
		l1, l2 :: List<numeric>
			# lists must be the same length
		tol :: Float
	OUTPUT:
		True if all elements in l1, l2 are equal within tol, else False
	'''
	ll1 = len(l1)
	if ll1 != len(l2):
		print '\n ERROR - listMatch() - Lists have unequal length!\n'
		quit()
	for i in xrange(ll1):
		if abs(l1[i] - l2[i]) >= tol:
			return False
	return True


def matMatch(m1, m2, tol=TOL):
	'''
	INPUT:
		m1, m2 :: List<List<numeric>>
			# matrices must have the same shapes
		tol :: Float
	OUTPUT:
		True if all elements in m1, m2 are equal within tol, else False
	'''
	lm1 = np.shape(m1); lm2 = np.shape(m2)
	if lm1[0] != lm2[0] or lm1[1] != lm2[1]:
		print '\n ERROR - matMatch() - Matrices have unequal shapes!\n'
		return False
	for i in xrange(lm1[0]):
		for j in xrange(lm1[1]):
			if abs(m1[i][j] - m2[i][j]) >= tol:
				return False
	return True


def stochTest(mat, tol=TOL):
	'''
	INPUT:
		mat :: List<List<Float>>
			# square matrix
		tol :: Float
	OUTPUT:
		True if columns sum to 1 within tol, else False
	'''
	s = np.shape(mat)
	for i in xrange(len(mat)):
		if abs(sum(mat[:,i]) - 1.0) >= tol:
			return False
	return True

