#
# Author: Flavio C. Coelho
# Source:
#   http://code.activestate.com/recipes/413086-gibbs-sampler/
#
from math import sqrt
import numpy as np
from matplotlib.pylab import * 
# initialize parameters for test cases
my_n = 10000 # number of iterations
my_k = 10 # number of components (parameters)
my_rho = .99 # correlation

#
# 2-COMPONENT GIBBS SAMPLER
#
def oneCompGibbs(n, k, rho):
	# initialization
	m1 = 10; m2 = 20 # means
	s1 = 1; s2 = 1 # standard deviations
	x = np.zeros(n); y = np.zeros(n)
	sd = sqrt(1-rho**2)
	# Sample recursively from two normal distributions.
	#   The mean for the current sample is updated in each step.
	for i in range(1,n):
		x[i] = np.random.normal(loc=m1+rho*(y[i-1]-m2)/s2, scale=s1*sd)
		y[i] = np.random.normal(loc=m2+rho*(x[i-1]-m1)/s1, scale=s2*sd)
	# graphing
	scatter(x,y,marker='d',c='r')
	title('2-Component Gibbs Sampler')
	xlabel('x')
	ylabel('y')
	grid()
	show()

#
# k-COMPONENT GIBBS SAMPLER
#
# (adapted from code above)
#
def kCompGibbs(n, k, rho):
	# initialization
	means = np.random.random_integers(10, high=20, size=k)
	stddevs = np.ones(k)
	params = np.zeros((k, n)) # numParams x numberOfIterations
	sd = sqrt(1-rho**2)
	# Sample recursively from two normal distributions.
	#   The mean for the current sample is updated in each step.
	for i in range(1,n):
		for j in range(k):
			nV = -1 if j+1 > k-1 else j+1 # index of next variable
			params[j][i] = np.random.normal( \
				loc=means[j]+rho*(params[nV][i-1]-means[nV])/stddevs[nV], \
				scale=stddevs[j]*sd \
				)
			#
			# ^THIS CONDITIONAL DISTRIBUTION MAY NOT BE RIGHT!
			#   (only conditioned on previous mean,
			#   not (all previous ones updated and all ones nonprevious ones unupdated))
			#
			# x[i] = np.random.normal(loc=m1+rho*(y[i-1]-m2)/s2, scale=s1*sd)
			# y[i] = np.random.normal(loc=m2+rho*(x[i-1]-m1)/s1, scale=s2*sd)
	# graphing 2 random parameters
	g = np.random.random_integers(0, high=k-1, size=2)
	scatter(params[g[0]], params[g[1]], marker='d',c='r')
	title(str(k)+'-Component Gibbs Sampler - showing parameters: '+str(g[0])+', '+str(g[1]))
	xlabel('x')
	ylabel('y')
	grid()
	show()

NUM_TRIALS = 10 # number of graphs to generate and show
for i in xrange(NUM_TRIALS):
	kCompGibbs(my_n, my_k, my_rho)

