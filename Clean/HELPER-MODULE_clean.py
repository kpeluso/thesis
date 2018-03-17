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

