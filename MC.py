#
# Author: Kenny Peluso
# Elevator Description: Sampling the final states of individuals moving about their own Markov Chains.
#

import numpy as np
from numpy.random import randint
from numpy import linalg as la
from sklearn.preprocessing import normalize
# http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.normalize.html
# https://stackoverflow.com/questions/29661574/normalize-numpy-array-columns-in-python

n = 50000 # number of individuals
iters = 5000 # number of timesteps a single individual takes within their Markov Chain

# ice cream example
P_I_iceCream = np.array([
	[7.0/8.0, 3.0/4.0],
	[1.0/8.0, 1.0/4.0]
]) # ice cream individual matrix
p11 = P_I_iceCream[0][0]; p21 = P_I_iceCream[0][1]
p12 = P_I_iceCream[1][0]; p22 = P_I_iceCream[1][1]
print 'Individual MC:'
print P_I_iceCream
eigenvalues, eigenvectors = la.eig(P_I_iceCream)
eigenvectors = normalize(eigenvectors, axis=0, norm='l1')
print eigenvalues; print eigenvectors

# P_I_iceCream's eigenvalues and eigenvectors:
# [ 1.     0.125]
# [[ 0.85714286 -0.5       ]
#  [ 0.14285714  0.5       ]]

# ice cream mass matrix
P_M_iceCream = np.transpose(np.array([
	[p11**3, 		3*p11**2*p12, 				3*p11*p12**2, 				p12**3],
	[p11**2*p21, 	p11**2*p22+2*p11*p12*p21, 	2*p12*p22*p11+p21*p12**2, 	p22*p12**2],
	[p11*p21**2, 	2*p11*p21*p22+p12*p21**2, 	2*p22*p12*p21+p11*p22**2, 	p22**2*p12],
	[p21**3, 		3*p21**2*p22, 				3*p21*p22**2, 				p22**3]
]))
print 'Mass MC:'
print P_M_iceCream
eigenvalues, eigenvectors = la.eig(P_M_iceCream)
eigenvectors = normalize(eigenvectors, axis=0, norm='l1')
print eigenvalues; print eigenvectors

# using the first Mass Eigenvector to weight all possible states
print "Weighted sum Mass States weighted by first Mass Eigenvector:"
result = 0.63*np.array([3,0]) + 0.315*np.array([2,1]) + 0.052*np.array([1,2]) + 0.003*np.array([0,3])
print result/sum(result)

quit() # ! ! !

# P_M_iceCream:
# [[ 0.66992188  0.57421875  0.4921875   0.421875  ]
#  [ 0.28710938  0.35546875  0.3984375   0.421875  ]
#  [ 0.04101562  0.06640625  0.1015625   0.140625  ]
#  [ 0.00195312  0.00390625  0.0078125   0.015625  ]]
# P_M_iceCream's eigenvalues and eigenvectors:
# [ 1.          0.125       0.015625    0.00195312]
# [[ 0.62973761  0.5         0.27272727 -0.125     ]
#  [ 0.3148688  -0.33333333 -0.5         0.375     ]
#  [ 0.05247813 -0.15277778  0.18181818 -0.375     ]
#  [ 0.00291545 -0.01388889  0.04545455  0.125     ]]
# ^first column sums to 1.177. Normalizing first column with 1.177 creates a probability distribution.

def sampleState(n_inds, t_steps, trans_mat, init_state=None):
	'''
	OUTPUT: List<Integer, Integer> - All sampled states any individual can possibly end up in as indexes
									with values as the number of individual that ended up there after t_steps.
	'''
	num_states = len(trans_mat)
	if not init_state: # use 1 init_state for all individuals
		init_state = randint(1, num_states)
	final_states = [0]*num_states # total numper of individuals ending up in states 1 and 2
	tr = np.transpose(trans_mat)
	for p in range(n_inds):
		curr_state = init_state-1
		for _ in range(t_steps):
			curr_state = np.random.choice(num_states, 1, p=tr[curr_state])[0]
		final_states[curr_state] += 1
		if not p%(n_inds/10): print 100*p/n_inds, '% of individuals done.' # print progress
	print '100 % Done.\n', final_states

# function calls / tests
sampleState(n, iters, P_I_iceCream) # [856, 144] when n=1000, [42752, 7248] when n=50000
# Therefore, seems to consistently be 85% and 15% in Individual State 1 and 2, respectively.

