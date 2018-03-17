# Source:
#   https://en.wikipedia.org/wiki/Power_iteration
# See also (for second derivative calculation):
#   https://math.stackexchange.com/questions/1114777/approximate-the-second-largest-eigenvalue-and-corresponding-eigenvector-given

# POWER ITERATION ALGORITHM

import numpy as np

NUM_SIMULATIONS = 50
test_mat = np.array([[0.5, 0.5], [0.2, 0.8]])

def firstEig(mat, num_simulations):
	# Ideally choose a random vector to decrease the chance
	#   that our vector is orthogonal to the eigenvector
	b_k = np.random.rand(mat.shape[0])
	for _ in range(num_simulations):
		b_k1 = np.dot(mat, b_k) # calculate the matrix-by-vector product Ab
		b_k1_norm = np.linalg.norm(b_k1) # calculate the norm
		b_k = b_k1 / b_k1_norm # re normalize the vector
	result = np.dot(mat, b_k)
	eigVa = 0.0
	l = len(result)
	for i in range(l):
		if result[i] or b_k[i]:
			eigVa = result[i]/b_k[i]
			if not eigVa:
				break
	return (b_k, eigVa)

def secondEig(mat, num_simulations):
	res1 = firstEig(mat, num_simulations)
	new_mat = np.abs(mat - res1[1]*np.outer(res1[0], res1[0]))
	return firstEig(new_mat, num_simulations)

# execute firstEig
result1 = firstEig(test_mat, NUM_SIMULATIONS)
print 'Eigenvector 1:', result1[0]
print 'Eigenvalue 1:', result1[1]

# execute secondEig
result2 = secondEig(test_mat, NUM_SIMULATIONS)
print 'Eigenvector 2:', result2[0]
print 'Eigenvalue 2:', result2[1]

# double-check
print 'Test Matrix:', test_mat
print result1[1]*np.outer(result1[0], result1[0]) \
	+ result2[1]*np.outer(result2[0], result2[0])
print result1[1]*np.outer(result1[0], result1[0])
print result2[1]*np.outer(result2[0], result2[0])

output = np.linalg.eig(test_mat)
print 'Eigenvector 2:', np.array([ \
	output[1][0][0]/output[1][1][0], \
	output[1][1][0]/output[1][1][0] \
	])

