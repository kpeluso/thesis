# Packages
# numpy already imported as np
from time import clock

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
    (2, np.array([.5,.5])), \
    (2, np.array([.25,.75])), \
    (2, np.array([.75,.25])), \
    (5, np.array([.2,.1,.3,.01,.39])), \
    (10, np.array([.2/2,.1/2,.3/2,.01/2,.39/2,.2/2,.1/2,.3/2,.01/2,.39/2])), \
    (15, np.array([.2/3,.1/3,.3/3,.01/3,.39/3,.2/3,.1/3,.3/3,.01/3,.39/3,.2/3,.1/3,.3/3,.01/3,.39/3])) \
]

# test cases to generate P^M - (N, n, P^I)
TESTS_PM = [ \
    (2, 2, np.array([[.25, .25],[.75, .75]])), \
    (14, 2, np.array([[.25, .65],[.75, .35]])), \
    (4, 4, np.array([[.25, .60, .75, .1], \
                    [.1, .1, .05, .2], \
                    [.3, .25, .05, .3], \
                    [.1, .1, .05, .2]])), \
    (14, 4, np.array([[.25, .60, .75, .1], \
                    [.1, .1, .05, .2], \
                    [.3, .25, .05, .3], \
                    [.1, .1, .05, .2]])), \
    (30, 4, np.array([[.25, .60, .75, .1], \
                    [.1, .1, .05, .2], \
                    [.3, .25, .05, .3], \
                    [.1, .1, .05, .2]])) \
]

# test suite
c = 0
PI_TIMES = [[],[],[],[],[],[],[],[]]
for T in TESTS_PI:
    print '\nTEST #' + str(c); c += 1
    start = clock(); print 'LP:', LP(T[0], T[1]); PI_TIMES[0].append(clock() - start)
    if c <= 2:
        start = clock(); print 'rS:', rS(T[0], T[1]); PI_TIMES[1].append(clock() - start)
        start = clock(); print 'bS:', bS(T[0], T[1]); PI_TIMES[2].append(clock() - start)
        start = clock(); print 'brS:', brS(T[0], T[1]); PI_TIMES[3].append(clock() - start)
    # start = clock(); print 'NRS:', NRS(T[0], T[1]); PI_TIMES[4].append(clock() - start)
    # start = clock(); print 'GRS:', GRS(T[0], T[1]); PI_TIMES[5].append(clock() - start)
    start = clock(); print 'GI:', GI(T[0], T[1]); PI_TIMES[6].append(clock() - start)
    # start = clock(); print 'CMAES:', CMAES(T[0], T[1]); PI_TIMES[7].append(clock() - start)
# the CMA-ES doesn't always converge in a pragmatic amount of time

c = 0
PM_TIMES = [[],[]]
for T in TESTS_PM:
    print '\nTEST #' + str(c); c += 1
    # start = clock(); ans = VSEA(T[0], T[1], T[2]); PM_TIMES[0].append(clock() - start)
    # print 'VSEA:', ans; print '\nVSEA verification:', len(ans), len(ans[0]), '\n'
    start = clock(); ans = Series(T[0], T[1], T[2]); PM_TIMES[1].append(clock() - start)
    # print 'Series:', ans; print '\nSeries verification:', len(ans), len(ans[0]), '\n'
# the VSEA doesn't always converge in a pragmatic amount of time

# print out realtime runtimes of each algorithm (in seconds)
print 'PI_TIMES:', PI_TIMES
print 'PM_TIMES:', PM_TIMES

