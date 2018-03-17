from ind2mass import *

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
	sl = int(comb(N+n-1, n-1)) # side length
	output = np.zeros([sl, sl])
	vectList = buildVectList(n, indMC) # lacks 0 vectors
	zeroVect = np.zeros(n)
	i = 0 # current row index
	for Ti in weak_compositions(N, n): # rows of mass matrix
		j = 0 # current column index
		for Tj in weak_compositions(N, n): # columns of mass matrix
			diff = np.array(Ti) - np.array(Tj) # Prospective state transition: Tj -> Ti
			vectList_copy = deepcopy(vectList)
			# add 0 vectors to vectList
			max_0s_per_bin = []
			for bn in xrange(n):
				max_0s_per_bin.append(min(Tj[bn], Ti[bn]))
				if Tj[bn]>0 and Ti[bn]>0:
					vectList_copy.append((np.zeros(n), indMC[bn][bn], bn))
			# get every way of allotting N tokens across all vectors in vectList
			lvlc = len(vectList_copy)
			transition_prob_terms = []
			for wc in weak_compositions(N, lvlc): # N = number of tokens, lvlc = number of vectors available
				sum_vect = np.zeros(n)
				breakIt = False
				static_in_bin = np.zeros(n) # number of people staying in their bin
				possibly_correct_vects = [] # :: List<Tuple< vector, prob, number_of_times_used, idx0 >>  # idx0 if 0 vect, gives corresponding diagonal index of particular 0 vect
				for vectTup_idx in xrange(lvlc):
					if wc[vectTup_idx] > 0: # check if vector is used
						# check if too many of the same zero vector has been invoked
						if listMatch(vectList_copy[vectTup_idx][0], zeroVect): # when we encounter a zero vector
							static_in_bin[vectList_copy[vectTup_idx][2]] += 1
							if max_0s_per_bin[vectList_copy[vectTup_idx][2]] < wc[vectTup_idx]:
								breakIt = True # too many of the same zero vector has been invoked
								break
						sum_vect += float(wc[vectTup_idx])*vectList_copy[vectTup_idx][0]
						# 	possibly_correct_vects.append(( \
						# 			vectList_copy[vectTup_idx][0], vectList_copy[vectTup_idx][1], \
						# 			wc[vectTup_idx], vectList_copy[vectTup_idx][2]))
						# else:
						possibly_correct_vects.append((vectList_copy[vectTup_idx][0], vectList_copy[vectTup_idx][1], wc[vectTup_idx]))
				if breakIt:
					continue
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
						# else: # when we encounter a zero vector
						# 	entering[pre_pvs[3]] = np.zeros(n)
					# check if too many people are supposed to be leaving a bin <-- REALLY HACKY! >:(
					if Ti == Tj:
						for bi in xrange(n):
							if leaving[bi] + static_in_bin[bi] > Tj[bi]:
								breakIt = True
								continue
						if breakIt:
							continue
					prod *= allotCombos(Tj, leaving, entering)
					# get probability factor
					for pvs in possibly_correct_vects:
						prod *= pvs[1]**pvs[2]
					transition_prob_terms.append(prod) # this is just one term that contributes to one entry in the mass matrix
			output[i][j] = sum(transition_prob_terms); j += 1 # add element to mass matrix
		i += 1
	return output