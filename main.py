#
# Author: Kenny Peluso
# Elevator Description: Create Python class objects to hold all necessary info and methods.
#

import numpy as np
from ind2mass import *

class Institution:
	'''
	DESCRIPTION:
		metastate of the entire system
	INPUT:
		MC_dict :: dict< String : Tuple< Integer, NPArray<NPArray<Float>> > >
		# String is from np.tostring(NPArray), Integer is the number of times that MC occurs
	'''

	def __init__(self, MCdict):
		self.MCdict = MCdict

	def addMC(MC):
		'''
		INPUT:
			MC :: NPArray<NPArray<Float>>
		'''
		theStr = np.tostring(MC)
		if theStr in self.MCdict:
			self.MCdict[theStr] += 1
		else:
			self.MCdict[theStr] = 1

	def delMC_once(MC):
		'''
		INPUT:
			MC :: NPArray<NPArray<Float>>
		'''
		theStr = np.tostring(MC)
		if theStr in self.MCdict and self.MCdict[theStr] > 1:
			self.MCdict[theStr] -= 1
		elif theStr in self.MCdict and self.MCdict[theStr] <= 1:
			del self.MCdict[theStr]
		else:
			print "\nError! ~~ delMC_once() ~~ Your input MC is not in the class dictionary!\n"

	def delMC_all(MC):
		'''
		INPUT:
			MC :: NPArray<NPArray<Float>>
		'''
		theStr = np.tostring(MC)
		if theStr in self.MCdict:
			del self.MCdict[theStr]
		else:
			print "\nError! ~~ delMC_once() ~~ Your input MC is not in the class dictionary!\n"

