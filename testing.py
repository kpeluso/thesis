# files to test
from ind2mass import *

# ind2mass.py dependencies:
# import numpy as np
# from scipy.special import comb
# from sklearn.utils.extmath import cartesian
# from sklearn.preprocessing import normalize
# from copy import deepcopy

print '\nTESTING:\n'

# COPY AND PASTE THIS FOR QUICK, IN-FILE TESTING:

# ans = GRS(3,np.array([.3,.3,.4]))
# print 'ans', ans
# ans2 = np.linalg.eig(ans)
# print 'ans2[0]', ans2[0]
# print 'eig vects:'
# print ans2[1][:,0]/sum(ans2[1][:,0])
# print ans2[1][:,1]/sum(ans2[1][:,1])
# print ans2[1][:,2]/sum(ans2[1][:,2])
#
# print '\n\n'
#
# ans = GRS(4,np.array([.2,.5,.2,.1]))
# print 'ans', ans
# ans2 = np.linalg.eig(ans)
# print 'ans2[0]', ans2[0]
# print 'eig vects:'
# print ans2[1][:,0]/sum(ans2[1][:,0])
# print ans2[1][:,1]/sum(ans2[1][:,1])
# print ans2[1][:,2]/sum(ans2[1][:,2])
# print ans2[1][:,3]/sum(ans2[1][:,3])


# print 'timeRev:'
# for i in xrange(5):
# 	out = genTimeRev(i+1)
# 	print out
# 	print stochTest(out)
# quit()

#timeRev(mat, stat, tol=TOL)

import difflib
from itertools import izip

def similarity(l1,l2):
	'''
	INPUT:
		l1,l2 :: List<Integer>
	OUTPUT:
		Float
			# similarity between l1,l2 as percentage
	'''
	# Source:
	#   https://stackoverflow.com/questions/6709693/calculating-the-similarity-of-two-lists
	sm = difflib.SequenceMatcher(None, l1, l2)
	return sm.ratio()

def edit_error(l1, l2):
	'''
	INPUT:
		l1,l2 :: List<Float>
			# pi and the first eigenvector of the calculated P^(I)
	OUTPUT:
		Float
			# similarity between the arrays of indices of l1,l2 sorted by the values of l1,l2
	'''
	# Source:
	#   https://stackoverflow.com/questions/4576115/convert-a-list-to-a-dictionary-in-python
	d1 = {}; c = 0
	for i in iter(l1):
		d1[c] = i
		c += 1
	d2 = {}; c = 0
	for i in iter(l2):
		d2[c] = i
		c += 1
	# print 'd1', d1
	# print 'd2', d2
	sorted_l1 = []; sorted_l2 = []
	# Source:
	#   https://www.saltycrane.com/blog/2007/09/how-to-sort-python-dictionary-by-keys/
	#for key in sorted(d1.iterkeys()):
	for key, value in sorted(d1.iteritems(), key=lambda (k,v): (v,k)):
		#sorted_l1.append((key, d1[key]))
		sorted_l1.append(key)
	#for key in sorted(d2.iterkeys()):
	for key, value in sorted(d2.iteritems(), key=lambda (k,v): (v,k)):
		#sorted_l2.append((key, d2[key]))
		sorted_l2.append(key)
	# print 'sorted_l1', sorted_l1
	# print 'sorted_l2', sorted_l2
	return similarity(sorted_l1, sorted_l2)


def sum_error(l1, l2, norm_type):
	'''
	INPUT:
		l1,l2 :: NPArray<Float>
			# pi and the first eigenvector of the calculated P^(I)
	OUTPUT:
		Float
			# 'norm_type'-norm of l1,l2
	'''
	return np.linalg.norm(l2-l1, norm_type)

print '\nError Function testing'
a = [1,1,2]; b = [1,1,3]; c = [1,2,1]
print edit_error(a,a)
print edit_error(b,b)
print edit_error(a,b)
print edit_error(b,a)
print edit_error(a,c)
print edit_error(c,a)
a = np.array(a); b = np.array(b)
print sum_error(a,a,1)
print sum_error(b,b,2)
print sum_error(a,b,1)
print sum_error(a,b,2)
quit()


print 'stochTest:'
mat1 = np.array([
	[.5, .6, .72],
	[.4, .0, .01],
	[.1, .4, .27],
	])
print 'Test 1:', stochTest(mat1) == True
mat2 = np.array([
	[.5, .6, .72],
	[.4, .1, .01],
	[.1, .4, .27],
	])
print 'Test 2:', stochTest(mat2) == False
mat1_ans = np.array( \
	[[ 0.66992188,  0.57421875,  0.4921875,   0.421875  ], \
	 [ 0.28710938,  0.35546875,  0.3984375,   0.421875  ], \
	 [ 0.04101562,  0.06640625,  0.1015625,   0.140625  ], \
	 [ 0.00195312,  0.00390625,  0.0078125,   0.015625  ]])
print 'Test 3:', stochTest(mat1_ans) == True # This matrix is the Mass MC for the Ice Cream Example


print 'matMatch:'
print 'Test 1:', matMatch(mat1, mat1) == True
print 'Test 2:', matMatch(mat1, mat2) == False
print 'Test 3:', matMatch(mat1_ans, mat1_ans) == True
print 'Test 4:', matMatch(mat1, mat1_ans) == False


print 'allotCombos:'
st = [3,1,0] # to [1,2,1]
leavs = [2, 0, 0]
ents = {0:[0, 1, 1]} # key in ents is 1 less than the corresponding bin's number
print 'Test 1:', allotCombos(st, leavs, ents) == 6
st = [3,1,0] # to [1,2,1]
leavs = [2, 1, 0]
ents = {0:[0, 2, 0], 1:[0,0,1]}
print 'Test 2:', allotCombos(st, leavs, ents) == 3
st = [2,1] # to [3,0]
leavs = [0, 1]
ents = {1:[1,0]}
print 'Test 3:', allotCombos(st, leavs, ents) == 1
st = [3,0] # to [2,1]
leavs = [1, 0]
ents = {0:[0,1]}
print 'Test 4:', allotCombos(st, leavs, ents) == 3


print 'ind2mass_nU:'
# ice cream example
N = 3; n = 2 # people, bins
mat1 = np.array([
	[7.0/8.0, 3.0/4.0],
	[1.0/8.0, 1.0/4.0]
	])
out = ind2mass_nU(N, n, mat1)
print 'Test 1a:', matMatch(out, mat1_ans) == True
print 'Test 1b:', stochTest(out) == True
# vary the individual matrix's values
N = 3; n = 2
mat1 = np.array([
	[2.0/8.0, 2.0/4.0],
	[6.0/8.0, 2.0/4.0]
	])
out = ind2mass_nU(N, n, mat1)
#print 'Test 2a:', out#matMatch(out, mat1_ans) == True
print 'Test 2b:', stochTest(out) == True
# increase the number of people in the ice cream example
N = 13; n = 2
mat1 = np.array([
	[2.0/8.0, 2.0/4.0],
	[6.0/8.0, 2.0/4.0]
	])
out = ind2mass_nU(N, n, mat1)
#print 'Test 3a:', out#matMatch(out, mat1_ans) == True
print 'Test 3b:', stochTest(out) == True
# large increase in number of people
N = 100; n = 2
mat1 = np.array([
	[2.0/8.0, 2.0/4.0],
	[6.0/8.0, 2.0/4.0]
	])
out = ind2mass_nU(N, n, mat1)
#print 'Test 4a:', out#matMatch(out, mat1_ans) == True
print 'Test 4b:', stochTest(out) == True

# increase the individual matrix's size
N = 13; n = 4
mat1 = np.array([
	[.5, .6, .25, .125],
	[.3, .0, .25, .3],
	[.1, .2, .25, .075],
	[.1, .2, .25, .5]
	])
out = ind2mass_nU(N, n, mat1)
#print 'Test 5a:', out#matMatch(out, mat1_ans) == True
print 'Test 5b:', stochTest(out) == True
print np.shape(out)
# m = np.shape(out)
# for i in xrange(m[1]):
# 	print out[:,i]
# 	print sum(out[:,i])
