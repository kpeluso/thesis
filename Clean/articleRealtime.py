# Packages
# numpy already imported as np
import sys
import pprint as pp
import matplotlib.pyplot as plt
from scipy.spatial.distance import cosine
from time import time
from HELPERS_clean import genStat, normize, listMatch, corrEv_m
from TIMEOUT_clean import timeout

# Algorithms
from GRS_clean import *
from GI_clean import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def graph(nfuns, domain, rng, labels, xl=None, yl=None, title=None, logx=False, logy=False):
    '''
    INPUT:
        nfuns :: Integer
        domain, rng :: List<Integer>
        labels :: List<String>
    '''
    lines = []
    if not (logx or logy):
        f = plt.plot
    elif logx and not logy:
        f = plt.semilogx
    elif not logx and logy:
        f = plt.semilogy
    elif logx and logy:
        f = plt.loglog
    for i in xrange(nfuns):
        aLine, = f(domain, rng[i], label=labels[i])
        lines.append(aLine)
    plt.legend(handles=lines)
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.title(title)
    plt.show()

def exactDist(v,w,tol=TOL):
    '''
    INPUT:
        v,w :: List<Integer>
            # equal lengths
            # v,w are randomly generated so likely no repeated values within either
    OUTPUT:
        Float
            # Exact match distance divided by n between v and a permutation of v, aka w
    '''
    total = 0
    n = len(v)
    for i in xrange(n):
        if abs(v[i]-w[i]) < tol:
            total += 1
    return float(n-total)/float(n)

def vals2rank(v):
    '''
    INPUT:
        v :: List<Float>
    OUTPUT:
        List<Integer>
            # Rank of each element in v, Rank 0 is smallest element
    '''
    big_val = np.max(v)
    big_idx = np.argmax(v)
    n = len(v)
    output = [n-1]*n
    for i in xrange(n):
        output[np.argmin(v)] = i
        v[np.argmin(v)] = big_val+1
    return output

def devDist(v,w):
    '''
    INPUT:
        v,w :: NPArray<Integer>
            # equal lengths
            # v,w are randomly generated so likely no repeated values within either
    OUTPUT:
        Float
            # Deviation distance between v and a permutation of v, aka w
    '''
    n = len(v)
    v_ranks = vals2rank(v)
    w_ranks = vals2rank(w)
    total = 0
    for i in xrange(n):
        total += abs(v_ranks.index(i) - w_ranks.index(i))
    return float(total)/float(n-1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Note that the digit-precision is constant between and within all test suites
# system arguments: 2, n2, ur, m
NUM_TIMES = 10 # number of trials = times we run each function (per each epsilon) per each test case
START_CASE = 2 # minimum value of n
MAX_SECS = 60 # number of seconds the algorithm has to complete per unit n
MAX_PIS = 7 # maximum value of n+1; a value of 4 will test 2 cases: n=2,3
ITERS = 10000
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
print 'REAL-RUNTIME and PROPORTION COMPLETED test suite for GRS, GI'

TEST_F = [GRS, GI]
D = ['GRS', 'GI']
results_1 = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F))] # time data
results_2 = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F))] # proportion of trials completed data

for f in xrange(len(TEST_F)):
    print '\nFUNCTION:', f, ',', D[f]
    for case in xrange(START_CASE, MAX_PIS):
        print '\t\tn_VALUE:', case
        TO_AVG_1 = [] # we will average the values in this array to get the average real-time runtime per n
        TO_AVG_2 = [] # we will average the values in this array to get the proportion of successful completions in 60n seconds
        for x in xrange(NUM_TIMES): # we'll average time for each test here
            test_stat = genStat(case)
            start = time()
            print '\t\t\tTRIAL:', x, ', pi_VALUE:', test_stat
            output = timeout(MAX_SECS*case, TEST_F[f], case, test_stat)
            TO_AVG_1.append(time() - start)
            print 'output', output
            TO_AVG_2.append(output)
        results_1[f][case-START_CASE] = np.average(TO_AVG_1)
        results_2[f][case-START_CASE] = np.average(TO_AVG_2)
print '\n', D
res1 = np.array(results_1)
res2 = np.array(results_2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
print 'LIMITED ACCURACY test suite for GRS, GI'

TEST_F = [GRS, GI]
D = ['GRS', 'GI']
results_ex = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F))] # accuracy data with exactDist
results_dev = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F))] # accuracy data with devDist
results_cos = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F))] # accuracy data with cosine
for f in xrange(len(TEST_F)):
    print '\nFUNCTION:', f
    for case in xrange(START_CASE, MAX_PIS):
        test_stat = genStat(case)
        print '\tn_VALUE:', case, ', pi_VALUE:', test_stat
        # we will average the values in these arrays to get the average accuracy per n
        TO_AVG_ex = []
        TO_AVG_dev = []
        TO_AVG_cos = []
        for x in xrange(NUM_TIMES): # we'll average time for each test here
            print '\t\tTRIAL:', x
            v = corrEv_m(TEST_F[f](case, test_stat, maxIters=ITERS*case))
            TO_AVG_ex.append(exactDist(v, test_stat))
            TO_AVG_dev.append(devDist(v, test_stat))
            TO_AVG_cos.append(cosine(v, test_stat))
        results_ex[f][case-START_CASE] = np.average(TO_AVG_ex)
        results_dev[f][case-START_CASE] = np.average(TO_AVG_dev)
        results_cos[f][case-START_CASE] = np.average(TO_AVG_cos)
print '\n', D
res_ex = np.array(results_ex)
res_dev = np.array(results_dev)
res_cos = np.array(results_cos)

for f in xrange(len(TEST_F)):
    print '\nRESULTS for FUNCTION', f
    print res1[f]
    print res2[f]
graph(len(TEST_F), range(START_CASE, MAX_PIS), res1, D, \
    xl='n Value', yl='Runtime (seconds)', title='\"n>=2\"-Case Direct Methods \n Real-Time Runtime')
graph(len(TEST_F), range(START_CASE, MAX_PIS), res2, D, \
    xl='n Value', yl='Proportion of '+str(NUM_TIMES)+' Trials Completed', title='\"n>=2\"-Case Direct Methods \n Proportion of Completed Trials')

for f in xrange(len(TEST_F)):
    print '\nRESULTS for FUNCTION', f
    print res_ex[f]
    print res_dev[f]
    print res_cos[f]
graph(len(TEST_F), range(START_CASE, MAX_PIS), res_ex, D, \
    xl='n Value', yl='Exact Distance Accuracy', title='\"n>=2\"-Case Direct Methods \n Exact Distance with '+str(ITERS)+'n Max Iterations')
graph(len(TEST_F), range(START_CASE, MAX_PIS), res_dev, D, \
    xl='n Value', yl='Deviation Distance Accuracy', title='\"n>=2\"-Case Direct Methods \n Deviation Distance with '+str(ITERS)+'n Max Iterations')
graph(len(TEST_F), range(START_CASE, MAX_PIS), res_cos, D, \
    xl='n Value', yl='Cosine Distance Accuracy', title='\"n>=2\"-Case Direct Methods \n Cosine Similarity with '+str(ITERS)+'n Max Iterations')

