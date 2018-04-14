# Packages
# numpy already imported as np
import sys
import pprint as pp
import matplotlib.pyplot as plt
from time import time
from HELPERS_clean import genStat, normize, listMatch
from TIMEOUT_clean import timeout

# Algorithms
from bS_clean import *
from rS_clean import *
from brS_clean import *
from LP_clean import *
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
    (2, np.array([.6943,.3057])), \
    (5, np.array([.2,.1,.3,.01,.39])), \
    (10, np.array([.2/2,.1/2,.3/2,.01/2,.39/2,.2/2,.1/2,.3/2,.01/2,.39/2])), \
    (15, np.array([.2/3,.1/3,.3/3,.01/3,.39/3,.2/3,.1/3,.3/3,.01/3,.39/3,.2/3,.1/3,.3/3,.01/3,.39/3])) \
]

# test cases to generate P^M - (N, n, P^I)
TESTS_PM = [ \
    (2, 2, np.array([[.25, .25],[.75, .75]])), \
    (10, 2, np.array([[.25, .65],[.75, .35]])), \
    (14, 2, np.array([[.25, .65],[.75, .35]])), \
    (21, 2, np.array([[.25, .65],[.75, .35]])), \
    (30, 2, np.array([[.25, .65],[.75, .35]])), \
    (32, 2, np.array([[.25, .65],[.75, .35]])), \

    (2, 3, np.array([[.2, .05, .01],[.7, .35, .9],[.1, .6, .09]])), \
    (2, 4, np.array([[.25, .60, .75, .1],[.1, .1, .05, .2],[.3, .25, .05, .3],[.1, .1, .05, .2]])), \
    (2, 5, np.array([[.3, .65, .05, .2, .4],[.1, .1, .05, .2, .05],[.4, .05, .8, .2, .1],[.1, .1, .05, .2, .4],[.1, .1, .05, .2, .05]])), \

    (4, 4, np.array([[.25, .60, .75, .1],[.1, .1, .05, .2],[.3, .25, .05, .3],[.1, .1, .05, .2]])), \
    (14, 4, np.array([[.25, .60, .75, .1],[.1, .1, .05, .2],[.3, .25, .05, .3],[.1, .1, .05, .2]]))
]

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

def expRT(n):
    '''
    Expected runtime of GI given n, the number of bins
    '''
    from scipy.special import digamma
    from mpmath import euler
    from scipy.misc import comb
    A = comb(n,2)
    return A*(digamma(A+1.0) + float(euler))

def divide(a,b):
    return float(a)/float(b)

vecexp = np.vectorize(expRT)
vecdiv = np.vectorize(divide)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Note that the digit-precision is constant between and within all test suites
# system arguments: 2, n2, ur, m
NUM_TIMES = 30 # number of trials = times we run each function (per each epsilon) per each test case
START_CASE = 2 # minimum value of n
if sys.argv[1] == '2':
    print '"n=2"-case direct methods test suite'
    #   vary epsilon
    #   plot log(TIME) x log(n)
    TEST_F = [bS, rS, brS]
    D = ['bS', 'rS', 'brS']
    NUM_EPSILONS = 10
    TEST_EPSILONS = [1.0/(10.0**i) for i in xrange(1,NUM_EPSILONS+1)]
    NUM_PIS = 4 # number of test cases (test pi values)
    results = [[[0]*NUM_PIS for x in xrange(NUM_EPSILONS)] for y in xrange(len(TEST_F))]
    for f in xrange(len(TEST_F)):
        print '\nFUNCTION:', f
        for T in xrange(len(TEST_EPSILONS)):
            print '\tEPSILON:', T, ', VALUE:', TEST_EPSILONS[T]
            for case in xrange(NUM_PIS):
                print '\t\tTESTCASE:', case, ', VALUE:', TESTS_PI[case][1]
                TO_AVG = [] # we will average the values in this array to get the per function per epsilon per test-case average runtime
                for x in xrange(NUM_TIMES): # we'll average time for each test here
                    start = time()
                    print '\t\t\tTRIAL:', x
                    TEST_F[f](2, TESTS_PI[case][1], tol=TEST_EPSILONS[T])
                    TO_AVG.append(time() - start)
                results[f][T][case] = np.average(TO_AVG)
    print D
    res = np.array(results)
    for f in xrange(len(TEST_F)):
        print '\nRESULTS for FUNCTION', f
        print res[f]
    plot_res = [[0]*NUM_EPSILONS for T in xrange(len(TEST_F))] # each column is an epsilon, each rows is a function
    for f in xrange(len(TEST_F)):
        for T in xrange(NUM_EPSILONS):
            plot_res[f][T] = np.average(res[f,T,:]) # average over all test cases
    for f in xrange(len(TEST_F)):
        print '\nGRAPHING DATA for FUNCTION', f
        print plot_res[f]
    graph(len(TEST_F), TEST_EPSILONS, plot_res, D, \
        xl='Epsilon Value', yl='Runtime (seconds)', title='\"n=2\"-Case Direct Methods', logx=True)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
elif sys.argv[1] == 'n2':
    print '"n>=2"-case direct methods test suite'
    #   vary n (cases from TESTS_PI)
    #   plot f(TIME) x f(n), where f is worst-case complexity of respective method
    TEST_F = [LP, GI, expRT]
    D = ['LP', 'GI', 'Expected GI']
    MAX_PIS = 7 # maximum value of n+1; a value of 4 will test 2 cases: n=2,3
    results = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F)-1)]
    varis = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F)-1)]
    for f in xrange(len(TEST_F)-1):
        print '\nFUNCTION:', f
        for case in xrange(START_CASE, MAX_PIS):
            print '\t\tn_VALUE:', case
            TO_AVG = [] # we will average the values in this array to get the per function per epsilon per test-case average runtime
            for x in xrange(NUM_TIMES): # we'll average time for each test here
                test_stat = genStat(case)
                start = time()
                print '\t\t\tTRIAL:', x, ', pi_VALUE:', test_stat
                print results
                r = TEST_F[f](case, test_stat)
                print 'RESULT:', r
                print 'Worked?:', listMatch(np.dot(r, test_stat), test_stat)
                TO_AVG.append(time() - start)
            results[f][case-START_CASE] = np.average(TO_AVG)
            varis[f][case-START_CASE] = np.var(TO_AVG)
    results.append(vecexp(range(START_CASE, MAX_PIS))) # get expected runtime of GI
    varis.append([None])
    print '\n', D
    res = np.array(results)
    v = np.array(varis)
    for f in xrange(len(TEST_F)):
        print '\nRESULTS for FUNCTION', f
        print res[f]
        print 'VARIANCES for FUNCTION', f
        print v[f]
    # only LP
    graph(1, range(START_CASE, MAX_PIS), [res[0]], [D[0]], \
        xl='n Value', yl='Runtime (seconds)', title='\"n>=2\"-Case Reliable Direct Methods')
    # graph GI and its expected runtime, 4 different ways
    graph(2, range(START_CASE, MAX_PIS), res[1:], D[1:], \
        xl='n Value', yl='Runtime (seconds)', title='\"n>=2\"-Case Reliable Direct Methods', logx=True, logy=True)
    graph(2, range(START_CASE, MAX_PIS), res[1:], D[1:], \
        xl='n Value', yl='Runtime (seconds)', title='\"n>=2\"-Case Reliable Direct Methods', logx=True, logy=False)
    graph(2, range(START_CASE, MAX_PIS), res[1:], D[1:], \
        xl='n Value', yl='Runtime (seconds)', title='\"n>=2\"-Case Reliable Direct Methods', logx=False, logy=True)
    graph(2, range(START_CASE, MAX_PIS), res[1:], D[1:], \
        xl='n Value', yl='Runtime (seconds)', title='\"n>=2\"-Case Reliable Direct Methods')
    # see difference between expected GI and actual GI real-time runtimes
    print '\nRESULTS for GI/E[GI]'
    print vecdiv(results[1],results[2])
    graph(1, range(START_CASE, MAX_PIS), [vecdiv(results[1],results[2])], ['Ratio GI/E[GI]'], \
        xl='n Value', yl='Ratio (GI)/(Expected GI) Realtime Runtime', title='\"n>=2\"-Case Direct Methods')
    # graph everything
    graph(2, range(START_CASE, MAX_PIS), res, D, \
        xl='n Value', yl='Runtime (seconds)', title='\"n>=2\"-Case Direct Methods', logx=False, logy=True)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
elif sys.argv[1] == 'ur':
    print 'un-reliable (NRS, GRS, CMAES) direct methods test suite'
    #   vary n (cases randomly generated from genStat())
    #   plot proportion_of_correct_answers_within_an_hour x test_case_index
    MAX_SECS = 60
    TEST_F = [NRS, GRS, CMAES]
    D = ['NRS', 'GRS', 'CMA-ES']
    MAX_PIS = 6 # maximum value of n+1; a value of 4 will test 2 cases: n=2,3
    results = [[0]*(MAX_PIS-START_CASE) for y in xrange(len(TEST_F))]
    for f in xrange(len(TEST_F)):
        print '\nFUNCTION:', f
        for case in xrange(START_CASE, MAX_PIS):
            test_stat = genStat(case)
            print '\tn_VALUE:', case, ', pi_VALUE:', test_stat
            TO_AVG = [] # we will average the values in this array to get the per function per epsilon per test-case average runtime
            for x in xrange(NUM_TIMES): # we'll average time for each test here
                print '\t\tTRIAL:', x
                TO_AVG.append(timeout(MAX_SECS*case, TEST_F[f], case, test_stat))
            results[f][case-START_CASE] = np.average(TO_AVG)
    print D
    res = np.array(results)
    for f in xrange(len(TEST_F)):
        print '\nRESULTS for FUNCTION', f
        print res[f]
    graph(len(TEST_F), range(START_CASE, MAX_PIS), res, D, \
        xl='n Value', yl='Proportion of '+str(NUM_TIMES)+' Successful Trials Completed in an Hour', title='Unreliable Direct Methods')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
elif sys.argv[1] == 'm':
    print 'P^(M) generating methods test suite'
    #   vary n and N (cases form TESTS_PM)
    #   plot TIME x test_case_index
    TEST_F = [VSEA, Series]
    D = ['VSEA', 'Series']
    NUM_TESTS = 10 # the number of tests to run: 0 to len(TESTS_PM)
    NUM_TESTS = min(max(NUM_TESTS,0), len(TESTS_PM)) # quality control

    results = [[0]*NUM_TESTS for y in xrange(len(TEST_F))]
    for f in xrange(len(TEST_F)):
        print '\nresults', results
        print '\nFUNCTION:', f
        for case in xrange(NUM_TESTS):
            print '\tTESTCASE:', case, ', CASE_VALUE:\n\t', TESTS_PM[case]
            start = time()
            TEST_F[f](TESTS_PM[case][0], TESTS_PM[case][1], TESTS_PM[case][2])
            results[f][case] = time() - start
    print D
    res = np.array(results)
    for f in xrange(len(TEST_F)):
        print '\nRESULTS for FUNCTION', f
        print res[f]
    graph(1,len(TEST_F), range(NUM_TESTS), res, D, \
        xl='Test Case', yl='Runtime (seconds)', title='Mass Matrix Generation Methods')
    graph(1, range(NUM_TESTS), [res[0]], [D[0]], \
        xl='Test Case', yl='Runtime (seconds)', title='Mass Matrix Generation Methods - VSEA')
    graph(1, range(NUM_TESTS), [res[1]], [D[1]], \
        xl='Test Case', yl='Runtime (seconds)', title='Mass Matrix Generation Methods - Series')

