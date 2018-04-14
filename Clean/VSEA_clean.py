from HELPERS_clean import TOL, np, weak_compositions, permute, perm, listMatch
from scipy.special import comb
from scipy.misc import factorial
from copy import deepcopy


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


def VSEA(N, n, indMC): # non-unique, identity-less individuals
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
		print 'Progress:', 100.0*float(sl*i)/float(sl**2), '%'
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

