from HELPER-MODULE import *
from scipy.optimize import linprog

# subproblem 6
c = -1.0*np.array([5, 2])
A_ub_val = np.array([
	[3, 1],
	[1, 1],
	])
b_ub_val = np.array([12, 5])
bounds_val = [(None, 3), (None, 3)]
print '\nSubproblem 6:'
print linprog(c, A_ub=A_ub_val, b_ub=b_ub_val, bounds=bounds_val)

def LP(n,pi):
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
    M = np.zeros([n,n**2])
    pi_l = list(pi)
    for i in xrange(n):
        M[i,:] += np.array([0]*(i*n) + pi_l + [0]*(n**2-i*n-n))
    # constraint 2 - Sum(P^(I)[:,j]) == 1 for all j=1:n
	for i in xrange(n,2*n):
        M[i,:] += np.array([0]*(i*n) + [1]*n + [0]*(n**2-i*n-n))
    # constraint 3 - P^(I)[i,j] >= 0 for all i,j=1:n
    bounds = [(0.0,1.0)]*n**2
    # print np.dot(devectorize(c,n), pi)
    return linprog(c, A_ub=M, b_ub=pi, bounds=bounds)
# print LP(5,np.array([.1,.2,.3,.4,.0]))

