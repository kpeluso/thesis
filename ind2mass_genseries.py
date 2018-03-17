
import numpy as np
from scipy.special import comb
from copy import deepcopy
from ind2mass import listMatch, stochTest

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
# for x in weak_compositions(3, 5): print x


# def bounded_wcs(balls, boxes, minBalls, maxBalls, currBox=0, parent=tuple()):
# 	'''
# 	Each bin i in a given generated output must have at least minBalls[i] balls
# 		but no more than maxBalls[i] balls.
# 	sum(minBalls) <= balls <= sum(maxBalls)
# 	'''
# 	if boxes > 1:
# 		i=max(0, balls-maxBalls[currBox]) # balls leftover after currBox
# 		while i <= balls-minBalls[currBox]:
# 			for x in bounded_wcs(i, boxes-1, minBalls, maxBalls, currBox+1, parent+(balls-i,)):
# 				yield x
# 			i+=1
# 	else:
# 		yield parent + (balls,)

# def bounded_wcs(balls, boxes, minBalls, maxBalls, currBox=0, parent=tuple(), first=True, iters=0):
# 	'''
# 	Each bin i in a given generated output must have at least minBalls[i] balls
# 		but no more than maxBalls[i] balls.
# 	sum(minBalls) <= balls <= sum(maxBalls)
# 	currBox indexed at 0
# 	'''
# 	num_funct_calls.append(iters)
# 	if first and sum(minBalls) > balls:
# 		print '\nError! :: Too few balls given the constraint minBalls!\n'; quit()
# 	if first and sum(maxBalls) < balls:
# 		print '\nError! :: Too many balls given the constraint maxBalls!\n'; quit()
# 	# check if current number of balls exceeds any max constraints yet to be parsed
# 	# 		This is the part that can be sped-up - I shouldn't have to implement
# 	#		this conditional in the first place.
# 	if balls > sum(maxBalls[currBox:]):
# 		pass
# 	else:
# 		if boxes > 1:
# 			i = max(0, balls-maxBalls[currBox]) # balls leftover after currBox
# 			while i <= balls-minBalls[currBox]:
# 				for x in bounded_wcs(i, boxes-1, minBalls, maxBalls, currBox+1, parent+(balls-i,), False, iters+1):
# 					yield x
# 				i+=1
# 		else:
# 			yield parent + (balls,)

def bounded_wcs(balls, boxes, minBalls, maxBalls, currBox=0, parent=tuple(), first=True, iters=0):
	'''
	Each bin i in a given generated output must have at least minBalls[i] balls
		but no more than maxBalls[i] balls.

	sum(minBalls) <= balls <= sum(maxBalls)

	currBox indexed at 0
	'''
	#num_funct_calls.append(iters)
	if first and sum(minBalls) > balls:
		print '\nError! :: Too few balls given the constraint minBalls!\n'; quit()
	if first and sum(maxBalls) < balls:
		print '\nError! :: Too many balls given the constraint maxBalls!\n'; quit()
	# calculate degrees of freedom and only vary within them
	#	(^unless this is already optimal)
	if balls > sum(maxBalls[currBox:]):
		pass
	else:
		if boxes > 1:
			i = max(0, balls-maxBalls[currBox]) # balls leftover after currBox
			while i <= balls-minBalls[currBox]:
				for x in bounded_wcs(i, boxes-1, minBalls, maxBalls, currBox+1, parent+(balls-i,), False, iters+1):
					yield x
				i+=1
		else:
			yield parent + (balls,)

# num_funct_calls = []
# m = [0,0,1,1,0]
# M = [0,2,3,3,3]
# c = 0
# for x in bounded_wcs(3, 5, m, M): print x; c += 1
# print 'count, m, M: ', c, m, M
# print num_funct_calls
# quit()
# import pdb; pdb.set_trace() # L6

# Would love to know:
# 1. How many generated outputs I ought to have?
# 2. Why my last output for >=6,5,m,M is not permissible?
# 	- It always adds more to very last bin
#	*- I should be sure to vary only balls per bins given a certain
#		number of degrees of freedom e.g. 'wiggle room' between the
#		number of balls and the maxBalls constraints
#	- PUT AN EXTRA CONDITION ON i IN/ON WHILE-LOOP ???


def mn(n, k):
	'''
	INPUT:
		n :: Integer
		k :: List<Integer>
			# sum(k) <= n
	OUTPUT:
		Integer
			# multinomial choose calculator i.e. multichoose
	'''
	result = 1
	less = 0
	for i in k:
		result *= comb(n-less,i)
		less += i
	return result
# print 'mn test output:'
# print mn(0,[0])
# print mn(1,[1])
# print mn(2,[1,1])
# print mn(2,[2,0])
# print mn(2,[0,2])
# print mn(2,[0,0,2,0])
# print mn(2,[1,0,1,0])
# print mn(2,[0,1,0,1])
# print mn(13,[3,4,5,1])
# quit()


# def p_wc(ls, initial_maxs=np.zeros(len(ls)).astype(int), initial_output=[], n=len(ls)):
# 	'''
# 	INPUT:
# 		ls :: NPArray<Integer>
# 	OUTPUT:
# 		List<NPArray<Integer>>
# 			# every possible successively bounded wc list as generator
# 	'''
# 	n = len(ls); n0s = np.zeros(n)
# 	for donor in ls:
# 		next_wc = np.array(bounded_wcs(donor, n, n0s, initial_maxs)).astype(int)
# 		initial_output.append(next_wc)
# 		initial_maxs -= next_wc


# def p_wc(ls, n, initial_maxs, initial_output=[], currIdx=0):
# 	'''
# 	INPUT:
# 		ls :: NPArray<Integer>
# 		n :: Integer
# 			# len(ls)
# 		initial_maxs :: NPArray<Integer>
# 			# should have length n
# 	OUTPUT:
# 		generated List<NPArray<Integer>>
# 			# a successively bounded wc list as generated
# 	'''
# 	# print ' '
# 	# print initial_output
# 	# print initial_maxs
# 	# print currIdx
# 	if currIdx == n or ls == []:
# 		yield initial_output
# 	else:
# 		n0s = np.zeros(n)
# 		for donor in xrange(n-currIdx):
# 			for bwc in bounded_wcs(ls[donor], n, n0s, initial_maxs):
# 				next_wc = np.array(bwc).astype(int)
# 				i_o = deepcopy(initial_output)
# 				i_o.append(next_wc)
# 				if currIdx == n-1:
# 					new_ls = []
# 				else:
# 					new_ls = ls[currIdx+1:]
# 				yield p_wc(new_ls, n, initial_maxs-next_wc, i_o, currIdx+1)

def p_wc(ls, n, initial_maxs, initial_output=[], currIdx=0):
	'''
	INPUT:
		ls :: NPArray<Integer>
		n :: Integer
			# len(ls)
		initial_maxs :: NPArray<Integer>
			# should have length len(ls)
	OUTPUT:
		generated List<NPArray<Integer>>
			# a successively bounded wc list as generated
	'''
	n0s = np.zeros(n)
	if currIdx == n or ls == [] or listMatch(initial_maxs, n0s):
		yield initial_output
	else:
		for bwc in bounded_wcs(ls[currIdx], n, n0s, initial_maxs):
			next_wc = np.array(bwc).astype(int)
			i_o = deepcopy(initial_output)
			i_o.append(next_wc)
			for x in p_wc(ls, n, initial_maxs-next_wc, i_o, currIdx+1):
				yield x
# print '\np_wc test output:'
# print '\ntest 1'
# ls = np.array([2,0,1]).astype(int)
# n = 3
# maxs = np.array([2,2,2]).astype(int)
# for i in p_wc(ls, n, maxs):
# 	print i
# print '\ntest 2'
# ls = np.array([2,1]).astype(int)
# n = 2
# maxs = np.array([2,2]).astype(int)
# for i in p_wc(ls, n, maxs):
# 	print i
# print '\ntest 3'
# ls = np.array([4,2,1]).astype(int)
# n = 3
# maxs = np.array([4,7,7]).astype(int)
# for i in p_wc(ls, n, maxs):
# 	print i
# print '\ntest 4'
# ls = np.array([4,2,1]).astype(int)
# n = 3
# maxs = np.array([2,4,1]).astype(int)
# for i in p_wc(ls, n, maxs):
# 	print i
# print '\ntest 5'
# ls = np.array([4,2,1]).astype(int)
# n = 3
# maxs = np.array([0,10,3]).astype(int)
# for i in p_wc(ls, n, maxs):
# 	print i
# quit()


wc_count = lambda n,k: int(comb(n+k-1,k-1)) # count total number of wcs, n=balls, k=boxes


def ind2mass_genseries(N, n, indMat):
	'''
	INPUT:
		N :: Integer
			# number of people
		n :: Integer
			# number of bins
		indMat :: List<List<Float>>
			# individual matrix
	OUTPUT:
		List<List<Float>>
			# mass matrix according to my generalized series formula (without Java speedup)
	'''
	sl = wc_count(N,n) # side length
	print 'Progress:'
	output = np.zeros([sl, sl])
	n0s = np.zeros([n, n])
	i = 0 # current row index
	for Ti in weak_compositions(N, n): # rows of mass matrix (s^(2))
		j = 0 # current column index
		for Tj in weak_compositions(N, n): # columns of mass matrix (s^(1))

			b1 = 0 # block 1
			for scwc in p_wc(np.array(Tj), n, np.array(Ti)): # loop over successively-constrained wcs
				b2 = 1 # block 2
				for d_idx, donor in enumerate(Tj):
					if donor > 0:
						b2 *= mn(donor, scwc[d_idx])
						for t in xrange(n): # block 3
							b2 *= indMat[t][d_idx]**scwc[d_idx][t]
				b1 += b2

			output[i][j] = b1
			j += 1
		i += 1
		print 100.0*float(i)/float(sl), '%' # progress output
	return output

# print 'ind2mass_genseries output:'
# print 'test 1:'
# N = 3; n = 2 # people, bins
# mat1 = np.array([
# 	[7.0/8.0, 3.0/4.0],
# 	[1.0/8.0, 1.0/4.0]
# 	])
# print 'side length:', wc_count(N,n)
# print stochTest(ind2mass_genseries(N, n, mat1))
# print 'test 2:',
# N = 13; n = 2
# mat1 = np.array([
# 	[2.0/8.0, 2.0/4.0],
# 	[6.0/8.0, 2.0/4.0]
# 	])
# print 'side length:', wc_count(N,n)
# print stochTest(ind2mass_genseries(N, n, mat1))
# print 'test 3',
# N = 200; n = 2
# mat1 = np.array([
# 	[2.0/8.0, 2.0/4.0],
# 	[6.0/8.0, 2.0/4.0]
# 	])
# print 'side length:', wc_count(N,n)
# print stochTest(ind2mass_genseries(N, n, mat1))

