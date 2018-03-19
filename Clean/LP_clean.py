import numpy as np
from scipy.optimize import linprog

def LP(n, pi):
    '''
	INPUT:
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via LP Method, solved with a black box Simplex Method
	'''
    # build objective function - trace(P^(I))
    c = []
    for i in xrange(n):
        c += [0]*i + [1] + [0]*(n-1-i)
    c = -1*np.array(c) # default is to minimize, so we multiply -1 to maximize
    # constraint 1 - InnerProduct(P^(I)[j,:], pi) == pi[j] for all j=1:n
    M = np.zeros([2*n,n**2])
    pi_l = list(pi)
    for i in xrange(n):
        M[i,:] += np.array([0]*(i*n) + pi_l + [0]*(n**2-i*n-n))
    # constraint 2 - Sum(P^(I)[:,j]) == 1 for all j=1:n
    for i in xrange(n,2*n):
        M[i,:] += np.array([0]*((i-n)*n) + [1]*n + [0]*(n**2-(i-n)*n-n))
    # constraint 3 - P^(I)[i,j] >= 0 for all i,j=1:n (nonnegativity constraints)
    bounds = [(0.0,1.0)]*n**2
    # Mx = b
    b = pi_l + [1]*n
    return linprog(c, A_ub=M, b_ub=b, bounds=bounds)

