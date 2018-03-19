# Packages
# numpy as np already imported

# Algorithms
from LP_clean import *
from rS_clean import *
from bS_clean import *
from brS_clean import *
from NRS_clean import *
from GRS_clean import *
from GI_clean import *
from CMAES_clean import *
from VSEA_clean import *
from Series_clean import *

# test cases to generate P^I - (n, pi)
TESTS_PI = [ \
    #(1, np.array([1])), \
    (2, np.array([.5,.5])), \
    (2, np.array([.25,.75])), \
    (2, np.array([.75,.25])), \
    (5, np.array([.2,.1,.3,.0,.4])) \
]

# test cases to generate P^M - (N, n, P^I)
TESTS_PM = [ \
    (1, 2, np.array([[.25, .25],[.75, .75]])), \
    (14, 2, np.array([[.25, .65],[.75, .35]])), \
    (4, 4, np.array([[.25, .60, .75, .1], \
                    [.1, .1, .05, .2], \
                    [.3, .25, .05, .3], \
                    [.1, .1, .05, .2]])), \
    (14, 4, np.array([[.25, .60, .75, .1], \
                    [.1, .1, .05, .2], \
                    [.3, .25, .05, .3], \
                    [.1, .1, .05, .2]])) \
]

# test suite
c = 0
for T in TESTS_PI:
    print '\nTEST #' + str(c); c += 1
    print 'LP:', LP(T[0], T[1])
    print 'rS:', rS(T[0], T[1])
    print 'bS:', bS(T[0], T[1])
    print 'brS:', brS(T[0], T[1])
    print 'NRS:', NRS(T[0], T[1])
    print 'GRS:', GRS(T[0], T[1])
    print 'GI:', GI(T[0], T[1])
    print 'CMAES:', CMAES(T[0], T[1])

c = 0
for T in TESTS_PM:
    print '\nTEST #' + str(c); c += 1
    ans = ind2mass_nU(T[0], T[1], T[2])
    print 'VSEA:', ans
    print '\nVSEA verification:', len(ans), len(ans[0]), '\n'
    ans = ind2mass_genseries(T[0], T[1], T[2])
    print 'Series:', ans
    print '\nSeries verification:', len(ans), len(ans[0]), '\n'

# !!!
# print out realtime runtimes of each algorithm
# !!!

