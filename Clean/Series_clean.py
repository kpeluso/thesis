from HELPERS_clean import np, listMatch, stochTest, wc_count, weak_compositions, comb
from copy import deepcopy

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


def Series(N, n, indMat):
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

