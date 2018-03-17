from ind2mass import listMatch, TOL, genStat, stochTest
from scipy.stats import truncnorm
import numpy as np


p22p1 = lambda p2,pi: 1.0 - (1.0 - p2)*(pi[1]/pi[0]) # proof of formula given in rough draft of thesis

vect = lambda mat: [mat[0][0], mat[1][0], mat[0][1], mat[1][1]]

next_term = lambda a,b: (b-a)*np.random.random() # b >= a


def resMat(mat):
	mat[1][0] = 1.0-mat[0][0]
	mat[0][1] = 1.0-mat[1][1]
	return mat


def normize(mat):
	for i in xrange(len(mat)):
		mat[:,i] /= sum(mat[:,i])
	return mat


def corrEv(mat, tol=TOL):
	'''
	INPUT:
		mat :: List<List<Float>>
			# 2x2
	OUTPUT:
		NPArray<Float>
			# returns eigenvector corresponding to eigenvalue of 1
	'''
	info = np.linalg.eig(mat)
	ev = info[1][:,0] if abs(info[0][0] - 1.0) < tol else info[1][:,1]
	return ev/sum(ev)


def corrEv_m(mat, tol=TOL):
	'''
	INPUT:
		mat :: List<List<Float>>
			# square
	OUTPUT:
		NPArray<Float>
			# returns eigenvector corresponding to eigenvalue of 1
	'''
	info = np.linalg.eig(mat)
	for i in xrange(len(info[1][:,0])):
		if abs(info[0][i] - 1.0) < tol:
			ev = info[1][:,i]
			return ev/sum(ev)


def brS(N, n, pi, iters=20000, tol=TOL):
	'''
	INPUT:
		N :: Integer
			# population size
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Binary-Random Search Method
			# only works in 2x2 case
			# BECAUSE IT'S 2x2 CASE I CAN USE 2x2 P^(M) SERIES FORMULA
	'''
	output = np.ones([2,2])*.5 # P^(I)
	ev = corrEv(output)
	b1 = 0.0 if pi[0] >= pi[1] else 1.0-(pi[0]/pi[1])
	b2 = 1.0
	output[1][1] = np.average([b1,b2])
	output = resMat(output)
	while not listMatch(np.dot(output, ev), pi) and iters > 0: # s1, loop
		# look for random number within updated b1, b2
			if ev[1] < pi[1]:
				b1 = np.average([b1,b2])
				output[1][1] += next_term(b1,b2)
			else:
				b2 = np.average([b1,b2])
				output[1][1] -= next_term(b1,b2)
			output[0][0] = p22p1(output[1][1], pi)
			output = resMat(output)
			ev = corrEv(output)
			iters -= 1
	return output

# print 'brS RESULTS:'
# testpi = np.array([2, 1])/3.0
# ans = brS(3,2,testpi)
# print ans
# ev = np.linalg.eig(ans)[1][:,0]; ev /= sum(ev)
# print ev
# print testpi
# quit()


def rS(N, n, pi, iters=20000, tol=TOL):
	'''
	INPUT:
		N :: Integer
			# population size
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Random Search Method
			# only works in 2x2 case
			# BECAUSE IT'S 2x2 CASE I CAN USE 2x2 P^(M) SERIES FORMULA
	'''
	output = np.zeros([2,2]) # P^(I)
	for col in xrange(2):
		output[:,col] = np.transpose(genStat(2))
	ev = corrEv(output)
	b1 = 0.0 if pi[0] >= pi[1] else 1.0-(pi[0]/pi[1])
	b2 = 1.0
	output[1][1] = np.average([b1,b2])
	output = resMat(output)
	while not listMatch(np.dot(output, ev), pi) and iters > 0: # s1, loop
		output[1][1] = (b2-b1)*np.random.random()+b1
		output[0][0] = p22p1(output[1][1], pi) # calculate p_11
		output = resMat(output) # calculate p_12, p_21
		ev = corrEv(output)
		iters -= 1
	return output

# print 'rS RESULTS:'
# testpi = np.array([2, 1])/3.0
# print rS(3,2,testpi)
# ev = np.linalg.eig(rS(3,2,testpi))[1][:,0]; ev /= sum(ev)
# print ev
# print testpi
# for i in xrange(10000):
# 	output = rS(3,2,testpi)
# 	for j in vect(output):
# 		if j<0:
# 			print output
# 			print 'Error!'
# 			quit()
# print '\nall good!\n'
# quit()


def bS(N, n, pi, iters=2000, tol=TOL):
	'''
	INPUT:
		N :: Integer
			# population size
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Binary Search Method
			# only works in 2x2 case
			# BECAUSE IT'S 2x2 CASE I CAN USE 2x2 P^(M) SERIES FORMULA
	'''
	output = np.ones([2,2])*0.5 # P^(I)
	ev = corrEv(output)
	b1 = 0.0 if pi[0] >= pi[1] else 1.0-(pi[0]/pi[1])
	b2 = 1.0
	output[1][1] = np.average([b1,b2])
	output = resMat(output)
	while not listMatch(np.dot(output, ev), pi) and iters > 0: # s1, loop
		if ev[1] < pi[1]:
			b1 = np.average([b1,b2])
			output[1][1] += np.average([b1,b2])-b1
		else:
			b2 = np.average([b1,b2])
			output[1][1] -= np.average([b1,b2])-b1
		output[0][0] = p22p1(output[1][1], pi)
		output = resMat(output)
		ev = corrEv(output)
		iters -= 1
	return output

print '\nbS RESULTS:'
# testpi = np.array([2, 1])/3.0
# print bS(3,2,testpi)
# ev = np.linalg.eig(bS(3,2,testpi))[1][:,0]; ev /= sum(ev)
# print ev
# print testpi
# quit()

# testpi = np.array([ 0.16666667,  0.83333333])
# out = bS(3,2,testpi)
# print out
# ev = corrEv(out)
# print ev
# print testpi
# quit()


def gibbs(N, n, pi, iters=20000, tol=TOL):
	'''
	INPUT:
		N :: Integer
			# population size
		n :: Integer
			# number of bins
		pi :: NPArray<Float>
			# optimal stationary distribution
			# CAN'T HAVE ANY ZERO ENTRIES!
	OUTPUT:
		NPArray<NPArray<Float>>
			# P^(I) via Gibbs Sampling-inspired Method
	'''
	output = np.zeros([n,n]) # P^(I)
	for col in xrange(n):
		output[:,col] = np.transpose(genStat(n))
	ev = corrEv(output)
	indices = range(n)

	# print '\n'
	# print iters
	# print 'np.dot(output, ev)', np.dot(output, ev)
	# print 'pi', pi
	# print 'len(np.dot(output, ev)), len(pi)', len(np.dot(output, ev)), len(pi)
	# print listMatch(np.dot(output, ev), pi)

	while not listMatch(np.dot(output, ev), pi) and iters > 0: # s1, loop
		# s2, isolate
		alter = np.random.choice(indices, size=[2], replace=False).astype(int)
		alter = np.array([min(alter), max(alter)]) # sort in order of lowest to highest
		subpi = np.zeros(2)
		subpi[0] = pi[alter[0]]
		subpi[1] = pi[alter[1]]
		# s3b, note how much space was just (formerly) taken up (call it Q)
		resMass_mat = (output[alter[0]][alter[0]] + output[alter[1]][alter[0]], \
			output[alter[0]][alter[1]] + output[alter[1]][alter[1]])
		resMass_pi = sum(subpi)
		# s3, normalize
		subpi /= sum(subpi)
		# s4, optimize (perhaps try again with Binary Search on 2x2)
		submat = brS(N, n, subpi) # !!! Use bS, rS, brS methods. !!!
		# s5a, denormalize
		submat[:,0] *= resMass_mat[0]
		submat[:,1] *= resMass_mat[1]
		subpi *= resMass_pi
		# s5, substitute in new values renormalized to Q
		output[alter[0]][alter[0]] = submat[0][0]
		output[alter[1]][alter[0]] = submat[1][0]
		output[alter[0]][alter[1]] = submat[0][1]
		output[alter[1]][alter[1]] = submat[1][1]
		ev = corrEv(output)

		iters -= 1
	return output

# print '\nGibbs output:'
# testpi = np.array([2, 1])/3.0
# N = 3
# n = 2
# rho = 0.4 # correlation
# print '\n'
# print 'test 1 results:'
# ans = gibbs(N, n, testpi, 10000*n)
# print 'ans', ans
# ev = corrEv_m(ans)
# print 'eigenvector:', ev
# print 'pi', testpi
# print 'mat is stoch', stochTest(ans)
# print 'ev is stoch', sum(ev)
#
# testpi = np.array([2, 10, 1])/13.0
# N = 13
# n = 3
# print '\n'
# ans = gibbs(N, n, testpi, 10000*n)
# print 'ans', ans
# ev = corrEv_m(ans)
# print 'test 2 results:'
# print 'eigenvector:', ev
# print 'pi', testpi
# print 'mat is stoch', stochTest(ans)
# print 'ev is stoch', sum(ev)
#
# testpi = np.array([23, 50, 7])/80.0
# N = 80
# n = 3
# print '\n'
# print 'test 3 results:'
# ans = gibbs(N, n, testpi, 10000*n)
# print 'ans', ans
# ev = corrEv_m(ans)
# print 'eigenvector:', ev
# print 'pi', testpi
# print 'mat is stoch', stochTest(ans)
# print 'ev is stoch', sum(ev)
#
# testpi = np.array([2, 5, 1, 1, 1, 1, 1, 1])/13.0
# N = 13
# n = 8
# print '\n'
# print 'test 4 results:'
# ans = gibbs(N, n, testpi, 10000*n)
# print 'ans', ans
# ev = corrEv_m(ans)
# print 'eigenvector:', ev
# print 'pi', testpi
# print 'mat is stoch', stochTest(ans)
# print 'ev is stoch', sum(ev)
#
