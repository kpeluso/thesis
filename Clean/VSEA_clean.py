import numpy as np
from scipy.special import comb
from sklearn.utils.extmath import cartesian
from sklearn.preprocessing import normalize
from copy import deepcopy

TOL = 0.001

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
# Probably from LeetCode (ask The Dude, again)
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


from scipy.misc import factorial
def perm(n, r):
  return factorial(n)/factorial(n-r)


def timeRev(mat, stat, tol=TOL):
	'''
	INPUT:
		mat :: List<List<Float>>
			# square transition matrix
		stat :: List<Float>
			# stationary distribution of mat
		tol :: Float
	OUTPUT:
		Boolean
			# True if mat is time-reversible (mat converges to mat of stationary dists), else False
	'''
	sl = np.shape(mat)
	for i in xrange(sl[0]): # loop through rows
		for j in xrange(sl[1]): # loop through columns
			if abs(mat[i][j] - mat[j][i]*stat[j]/stat[i]) >= tol:
				return False
	return True


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


def massEntry_2x2(s1, s2, indMat):
	'''
	INPUT:
		s1, s2 :: Tuple< Integer, Integer>
			# start state, end state
		indMat :: List<List<Float>>
			# individual matrix
			# ASSUMES: 2x2
	OUTPUT:
		Float
			# mass matrix entry corresponding to s1, s2
			# calculated according to my series formula
	'''
	N = sum(s1)
	m = min(s1)
	if s1[0] <= s1[1]:
		m2 = s2[0]
		probs = lambda x: (indMat[0][0]**x) \
							*(indMat[1][0]**(m-x)) \
							*(indMat[0][1]**(s2[0]-x)) \
							*(indMat[1][1]**(s2[1]-m+x))
	else:
		m2 = s2[1]
		probs = lambda x: (indMat[0][0]**(s2[0]-m+x)) \
							*(indMat[1][0]**(s2[1]-x)) \
							*(indMat[0][1]**(m-x)) \
							*(indMat[1][1]**x)
	total = 0.0
	for i in xrange(m+1):
		total += comb(N-m,m2-i)*probs(i)
		if s1 == (3, 0) and s2 == (1, 2):
			print ' '
			print m
			print N
			print m2
			print probs(i)
			print comb(N-m,m2-i)
			print comb(N-m,m2-i)*probs(i)
	return total


def ind2mass_series(N, n, indMat):
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
			# mass matrix according to my series formula
			# ice cream example only!
	'''
	sl = int(comb(N+n-1, n-1)) # side length
	output = np.zeros([sl, sl])
	i = 0 # current row index
	for Ti in weak_compositions(N, n): # rows of mass matrix
		j = 0 # current column index
		for Tj in weak_compositions(N, n): # columns of mass matrix
			output[i][j] = massEntry_2x2(Tj,Ti,iceCream)
			j+=1
		i += 1
	return output


def buildVectList(size, mat):
	'''
	INPUT:
		size :: Integer
			# n
		mat :: List<List<Integer>>
			# individual matrix
	OUTPUT:
		List< Tuple< List<Integer>, Float > >
			# first index is permutation of (n-2)*0,1,-1, next is P[source -> sink]
	'''
	output = []
	for perm in iter(permute([1,-1]+[0]*(size-2))):
		oneIdx = perm.index(1); negOneIdx = perm.index(-1)
		output.append((np.array(perm), mat[oneIdx][negOneIdx]))
	return output


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


def allotCombos(st, leavs, ents):
	'''
	INPUTS:
		st :: List<Integer>
			# number of people in each bin in starting state (Tj)
		leavs :: List<Ineger>
			# N-length list of people leaving each bin
		ents :: dict<Integer:List<Integer>>
			# key=source bin : value=N-length list of numbers of unit sinks per bin
	OUTPUT:
		Integer
			# total number of ways of sorting people into bins given variable, known space in each bin
	'''
	output = 1
	for idx, bi in enumerate(st):
		output *= comb(bi, leavs[idx])
		if bi > 0 and output > 0 and leavs[idx] > 0:
			sinked = 0 # number of people already allotted
			for bi2 in ents[idx]:
				output *= comb(leavs[idx]-sinked, bi2)
				sinked += bi2
	return output


def ind2mass_nU(N, n, indMC): # non-unique, identity-less individuals
	'''
	INPUT:
		N :: Integer
			# number of people
		n :: Integer
			# number of bins
		indMC :: NPArray<NPArray<Float>>
			# individual matrix
	OUTPUT:
		List<List<Float>>
			# mass matrix
	'''
	#
	# Block A
	#
	sl = int(comb(N+n-1, n-1)) # side length
	output = np.zeros([sl, sl])
	TESTLIST = []
	vectList = buildVectList(n, indMC) # lacks 0 vectors
	zeroVect = np.zeros(n)
	#
	# Block B
	#
	i = 0 # current row index
	for Ti in weak_compositions(N, n): # rows of mass matrix
		j = 0 # current column index
		for Tj in weak_compositions(N, n): # columns of mass matrix
			#
			# Block C
			#
			diff = np.array(Tj) - np.array(Ti) # Prospective state transition: Tj -> Ti
			vectList_copy = deepcopy(vectList)
			# add 0 vectors to vectList
			max_0s_per_bin = []
			for bn in xrange(n):
				max_0s_per_bin.append(min(Tj[bn], Ti[bn]))
				if Tj[bn]>0 and Ti[bn]>0:
					vectList_copy.append((np.zeros(n), indMC[bn][bn], bn))
			#
			# Block D
			#
			# get every way of allotting N tokens across all vectors in vectList
			lvlc = len(vectList_copy)
			transition_prob_terms = []
			for wc in weak_compositions(N, lvlc): # N = number of tokens, lvlc = number of vectors available
				sum_vect = np.zeros(n)
				breakIt = False
				static_in_bin = np.zeros(n) # number of people staying in their bin
				possibly_correct_vects = [] # :: List<Tuple< vector, prob, number_of_times_used, idx0 >>  # idx0 if 0 vect, gives corresponding diagonal index of particular 0 vect
				for vectTup_idx in xrange(lvlc):
					if wc[vectTup_idx] > 0: #check if vector is used
						# check if too many of the same zero vector has been invoked
						if listMatch(vectList_copy[vectTup_idx][0], zeroVect): # when we encounter a zero vector
							static_in_bin[vectList_copy[vectTup_idx][2]] += wc[vectTup_idx]
							if max_0s_per_bin[vectList_copy[vectTup_idx][2]] < wc[vectTup_idx]:
								breakIt = True # too many of the same zero vector has been invoked
								break
						sum_vect += float(wc[vectTup_idx])*vectList_copy[vectTup_idx][0]
						possibly_correct_vects.append((vectList_copy[vectTup_idx][0], vectList_copy[vectTup_idx][1], wc[vectTup_idx]))
				if breakIt: # skip to next wc
					continue
				#
				# Block E
				#
				if listMatch(diff, sum_vect): # check if current weak composition wc of vectors equals diff
					# coefficient determination
					leaving = np.zeros(n); entering = {}; prod = 1
					for pre_pvs in possibly_correct_vects:
						pvs = np.ndarray.tolist(pre_pvs[0])
						if not listMatch(pvs, zeroVect): # when we DON'T encounter a zero vector
							for bi in xrange(n):
								if pvs[bi] == -1:
									leaving[bi] += pre_pvs[2]
									if bi not in entering:
										entering[bi] = np.zeros(n)
									entering[bi][pvs.index(1)] += pre_pvs[2]
									break # jump to next pre_pvs
					# check if too many people are leaving a bin <-- REALLY HACKY! >:(
					for bi in xrange(n):
						if leaving[bi] + static_in_bin[bi] > Ti[bi]:
							breakIt = True
							break
					if breakIt: # skip to next wc
						continue
					#
					# Block F
					#
					prod *= allotCombos(Ti, leaving, entering) # get coefficient
					# get P^I probabilities that correspond to selected basis vectors
					#   get exponents for all probabilities that correspond to the
					#   quantities of each of their respective basis vector
					pvstEST = 1
					for pvs in possibly_correct_vects:
						prod *= pvs[1]**pvs[2]
						pvstEST *= pvs[1]**pvs[2]
					transition_prob_terms.append(prod) # this is just one term that contributes to one entry in the mass matrix
			output[j][i] = sum(transition_prob_terms); j += 1 # add element to mass matrix
		i += 1
	return output

