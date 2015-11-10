import matplotlib.pyplot
#from openalea.stat_tool.distribution import DistributionIdentifierType as f
from _dmcdistributions import _DiscreteMultivariateConditionalDistribution
#from _dmcdistributions import DUCPLinks as Links
from _dmcdistributions import DMCDistributions as Ids
from _dmcdistributions import _DiscreteMultivariateConditionalParametric
from _dmcdistributions import _DiscreteMultivariateConditionalData
from openalea.tree_statistic.ducdistributions import DiscreteUnivariateConditionalParametric

class DiscreteMultivariateConditionalDistribution(_DiscreteMultivariateConditionalDistribution):
	"""
		A class for Discrete Multivariate Conditional Distributions in python
	"""

	def __init__(self, *args):
		"""
			Initialize a DiscreteMultivariateConditionalDistribution

				:parameter:
					`nb_covariables` (int) - Number of covariables for DiscreteMultivariateConditionalDistribution
					`nb_variables` (int) - Number of variables for DiscreteMultivariateConditionalDistribution
					`mass` (list) : A list of covariates and response mass
				
				:parameter:
					`ducd` (DiscreteMultivariateConditionalDistribution) - A DiscreteMultivariateConditionalDistribution to copy

				:usage:
					mass = [[[0], [[0,0],[0,1],[1,0],[1,1]],[0.25,0.25,0.25,0.25]]] 
					d = DiscreteMultivariateConditionalDistribution(1, 2, mass) # initialize a DiscreteMultivariateConditionalDistribution from masses
					dbis = DiscreteMultivariateConditionalDistribution(d) # copy a DiscreteMultivariateCondtionalDistribution
		"""
		try:
			if len(args) == 3:
				if not(isinstance(args[0], int)):
					raise Exception('Error: Invalid `nb_covariables` parameter !')
				elif not(args[0] > 0):
					raise Exception('Error: Invalid `nb_covariables` parameter !')
				elif not(isinstance(args[1], int)):
					raise Exception('Error: Invalid `nb_variables` parameter !')
				elif not(args[1] > 0):
					raise Exception('Error: Invalid `nb_variables` parameter !')
				elif not(isinstance(args[2], list)):
					raise Exception('Error: Invalid	`mass` parameter !')
				else:
					for i in args[2]:
						if not(isinstance(i, list)):
							raise Exception('Error: Invalid `mass` parameter !')
						elif not(len(i) == 3):
							raise Exception('Error: Invalid `mass` parameter !')
						elif not(isinstance(i[0], list)):
							raise Exception('Error: Invalid covariates values (list) !')
						elif not(len(i[0]) == args[0]):
							raise Exception('Error: Wrong number of covariates !')
						else:
							for j in i[0]:
								if not(isinstance(j, int)):
									raise Exception('Error: Invalid covariates value (list of int)')
							if not(len(i[1]) == len(i[2])):
								raise Exception('Error: variables events number and probabilies number not equal !')
							for j in i[1]:
								if not(isinstance(j, list)):
									raise Exception('Error: Invalid variable event !')
								elif not(len(j) == args[1]):
									raise Exception('Error: Wrong number of variables !')
								else:
									for k in j:
										if not(isinstance(k ,int)):
											raise Exception('Error: Invalid variable event !')
							for j in i[2]:
								if not(isinstance(j, float)):
									raise Exception('Error: Invalid probability (float) !')
			elif len(args) == 1:
				if not(_DiscreteMultivariateConditionalDistribution.__instancecheck__(args[0]) or DiscreteMultivariateConditionalDistribution.__instancecheck__(args[0])):
					raise Exception('Error: `dmcd` parameter not valid !')
			else:
				raise Exception('Error: Invalid parameters for DiscreteMultivariateConditionalDistribution initialization !')
		except Exception as e:
			print e
		else:
			_DiscreteMultivariateConditionalDistribution.__init__(self, *args)

	def GetNbCovariables(self):
		"""
			Return number of covariables for DiscreteMultivariateConditionalDistribution

				:usage:
					d.GetNbCovariables()
		"""
		return _DiscreteMultivariateConditionalDistribution._GetNbCovariables(self)

	def GetNbVariables(self):
		"""
			Return number of variables for DiscreteMultivariateConditionalDistribution

				:usage:
					d.GetNbVariables()
		"""
		return _DiscreteMultivariateConditionalDistribution._GetNbVariables(self)

	def GetNbParameters(self):
		"""
			Return number of parameters for DiscreteMultivariateConditionalDistribution

				:usage:
					d.GetNbParameters()
		"""
		return _DiscreteMultivariateConditionalDistribution._GetNbParameters(self)

	def GetMass(self, covariates, responses):
		"""
			Return mass of event for DiscreteMultivariateConditionalDistribution

				:parameter:
					`covariates` (list) : Covariates for event
					`responses` (list) : Response variable event

				:usage:
					d.GetMass([0],[1,1])
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			elif not(len(covariates) == self.GetNbCovariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbCovariables()) + ' covariables but got ' +str(len(covariates))+ ' !')
			elif not(isinstance(responses, list)):
				raise Exception('Error: Invalid `responses` parameter !')
			elif not(len(responses) == self.GetNbVariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbVariables()) + ' covariables but got ' +str(len(responses))+ ' !')
			else:
				for i in covariates:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
				for i in responses:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `responses` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalDistribution._GetMass(self, covariates, responses)

	def Simulate(self, covariates):
		"""
			TODO
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			for i in covariates:
				if not(isinstance(i, list)):
					raise Exception('Error: Invalid `covariates` parameter !')
				if not(len(i) == self.GetNbCovariables()):
					raise Exception('Error: Invalid `covariates` parameter !')
				for j in i:
					if not(isinstance(j, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
		except Exception as e:
			print e
		else:
			return DiscreteMultivariateConditionalData(_DiscreteMultivariateConditionalDistribution._Simulate(self, covariates))

	def GetMeans(self, covariates):
		"""
			Return mean with covariates for DiscreteMultivariateConditionalDistribution

				:parameter:
					`covariates` (list) : Covariates

				:usage:
					d.GetMeans([0])
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			elif not(len(covariates) == self.GetNbCovariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbCovariables()) + ' covariables but got ' +str(len(covariates))+ ' !')
			else:
				for i in covariates:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalDistribution._GetMeans(self, covariates)

	def GetVarianceCovarianceMatrix(self, covariates):
		"""
			Return variance with covariates for DiscreteMultivariateConditionalDistribution

				:parameter:
					`covariates` (list) : Covariates

				:usage:
					d.GetVarianceCovarianceMatrix([0])
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			elif not(len(covariates) == self.GetNbCovariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbCovariables()) + ' covariables but got ' +str(len(covariates))+ ' !')
			else:
				for i in covariates:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalDistribution._GetVariancesCovariances(self, covariates)

class DiscreteMultivariateConditionalParametric(_DiscreteMultivariateConditionalParametric):
	"""
		A class for Discrete Multivariate Conditional Parametric Distributions in python
	"""

	def __init__(self, *args):
		"""
			Initialize a DiscreteMultivariateConditionalParametric

				:parameter:
					`nb_covariables` (int) : Number of covariables.
					`nb_variables` (int) : Number of variables.
					`ident` : an identifier DiscreteMultivariateConditionalParametric distribution
					`ducp` (DiscreteUnivariateConditionalParametric) : A compound distribution
					`parameters` (list) : A list of parameters

				:parameter:
					`dmcp` (DiscreteMultivariateConditionalParametric) - A DiscreteMultivariateConditionalParametric to copy

				:usage:
					p = DiscreteMultivariateConditionalParametric(2, 2, d, Ids.MultinomialSumCompound, [0.2,0.8]) # initialize a DiscreteUnivariateConditionalParametric
					pbis = DiscreteMultivariateConditionalParametric(p) # initialize a DiscreteMultivariateCondtionalParametric by copy
		"""
		try:
			if len(args) == 5:
				if not(isinstance(args[0], int)):
					raise Exception('Error: Invalid `nb_covariables` parameter !')
				elif not(isinstance(args[1], int)):
					raise Exception('Error: Invalid `nb_variables` parameter !')
				elif not(DiscreteUnivariateConditionalParametric.__instancecheck__(args[3])):
					raise Exception('Error: Invalid `ducp` parameter !')
				elif not(args[3].GetNbCovariables() == args[0]):
					raise Exception('Error: Expecting '+ str(args[0])+ ' covariables but `ducd` has '+str(args[3].GetNbCovariables()) + ' covariables !')
				elif not(Ids.__instancecheck__(args[2])):
					raise Exception('Error: Invalid `ident` parameter !')
				elif not(isinstance(args[4], list)):
					raise Exception('Error: Invalid `parameters` parameter !')
				else:
					for i in args[4]:
						if not(isinstance(i, float)):
							raise Exception('Error: Invalid `parameters` parameter')
			elif len(args) == 1:
				if not(_DiscreteMultivariateConditionalParametric.__instancecheck__(args[0]) or DiscreteMultivariateConditionalParametric.__instancecheck__(args[0])):
					raise Exception('Error: `dmcp` parameter not valid !')
			else:
				raise Exception('Error: Invalid parameters for DiscreteUnivariateConditionalParametric initialization !')
		except Exception as e:
			print e
		else:
			_DiscreteMultivariateConditionalParametric.__init__(self, *args)

	def GetCompound(self):
		"""
			Return Compound distribution of DiscreteMultivariateConditionalParametric

				:usage:
					p.GetCompound()
		"""
		return DiscreteUnivariateConditionalParametric(_DiscreteMultivariateConditionalParametric._GetCompound(self))

	def GetNbCovariables(self):
		"""
			Return number of covariables for DiscreteMultivariateConditionalParametric
		"""
		return _DiscreteMultivariateConditionalParametric._GetNbCovariables(self)

	def GetNbVariables(self):
		"""
			Return number of variables for DiscreteMultivariateConditionalParametric
		"""
		return _DiscreteMultivariateConditionalParametric._GetNbVariables(self)

	def GetNbParameters(self):
		"""
			Return number of parameters for DiscreteUnivariateConditionalParametric
		"""
		return _DiscreteMultivariateConditionalParametric._GetNbParameters(self)

	def GetParameters(self):
		"""
			Return parameters of DiscreteUnivariateConditionalParametric
		"""
		return _DiscreteMultivariateConditionalParametric._GetParameters(self)

	def GetMass(self, covariates, responses):
		"""
			Return mass of event for DiscreteUnivariateConditionalParametric

				:parameter:
					`covariates` (list) : Covariates for event
					`responses` (list) : Responses variables event

				:usage:
					d.GetMass([0,1],[0,1])
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			elif not(len(covariates) == self.GetNbCovariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbCovariables()) + ' covariables but got ' +str(len(covariates))+ ' !')
			if not(isinstance(responses, list)):
				raise Exception('Error: Invalid `responses` parameter !')
			elif not(len(responses) == self.GetVariables()):
				raise Exception('Error: Expecting '+ str(self.GetVariables()) + ' variables but got ' +str(len(responses))+ ' !')
			else:
				for i in covariates:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
				for i in responses:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `responses` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalParametric._GetMass(self, covariates, responses)

	def Simulate(self, covariates):
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			for i in covariates:
				if not(isinstance(i, list)):
					raise Exception('Error: Invalid `covariates` parameter !')
				if not(len(i) == self.GetNbCovariables()):
					raise Exception('Error: Invalid `covariates` parameter !')
				for j in i:
					if not(isinstance(j, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
		except Exception as e:
			print e
		else:
			return DiscreteMultivariateConditionalData(_DiscreteMultivariateConditionalParametric._Simulate(self, covariates))

	def GetMeans(self, covariates):
		"""
			Return means with covariates for DiscreteMultivariateConditionalParametric

				:parameter:
					`covariates` (list) : Covariates

				:usage:
					d.GetMeans([0,0])
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			elif not(len(covariates) == self.GetNbCovariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbCovariables()) + ' covariables but got ' +str(len(covariates))+ ' !')
			else:
				for i in covariates:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalParametric._GetMeans(self, covariates)

	def GetVarianceCovarianceMatrix(self, covariates):
		"""
			Return variance co-variance matrix with covariates for DiscreteMultivariateConditionalDistribution

				:parameter:
					`covariates` (list) : Covariates

				:usage:
					d.GetVarianceCovarianceMatrix([0,0])
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			elif not(len(covariates) == self.GetNbCovariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbCovariables()) + ' covariables but got ' +str(len(covariates))+ ' !')
			else:
				for i in covariates:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalParametric._GetVariancesCovariances(self, covariates)

class DiscreteMultivariateConditionalData(_DiscreteMultivariateConditionalData):
	"""
		A class for Discrete Multivariate Conditional Data estimation in Python
	"""

	def __init__(self, *args):
		"""
			Initialize a DiscreteMultivariateConditionalData

				:parameter:
					`nb_covariables` (int) - Number of covariables
					`nb_variables` (int) - Number of variables
					`covariates` (list) - A list of covariates values
					`responses` (list) - Corresponding responses values

				:parameter:
					`vectors` (Vectors) - 

				:usage:
					c = [[0],[0]]
					r = [[0,1],[1,0]]
					e = DiscreteMultivariateConditionalData(1,2,c,r)
		"""
		try:
			if len(args) == 4:
				if not(isinstance(args[0], int)):
					raise Exception('Error: Invalid `nb_covariables` parameter !')
				elif not(args[0] > 0):
					raise Exception('Error: Invalid `nb_covariables` parameter !')
				elif not(isinstance(args[1], int)):
					raise Exception('Error: Invalid `nb_variables` parameter !')
				elif not(args[1] > 0):
					raise Exception('Error: Invalid `nb_variables` parameter !')
				elif not(isinstance(args[2], list) and isinstance(args[3], list)):
					raise Exception('Error: Invalid `covariates` (nested list) or `response` (nested list) parameter !')
				elif not(len(args[2]) == len(args[3])):
					raise Exception('Error: Invalid `covariates` or `response` parameter (not same number of observations) !')
				else:
					for i in args[2]:
						if not(isinstance(i, list)):
							raise Exception('Error: Invalid `covariates` parameter (nested list)')
						elif not(len(i) == args[0]):
							raise Exception('Error: Expecting '+str(args[0])+ ' covariables but got '+str(len(i))+ ' !')
						else:
							for j in i:
								if not(isinstance(j, int)):
									raise Exception('Error: Invalid `covariates` parameter (nested list of int)!')
					for i in args[3]:
						if not(isinstance(i, list)):
							raise Exception('Error: Invalid `responses` parameter (nested list)')
						elif not(len(i) == args[1]):
							raise Exception('Error: Expecting '+str(args[1])+ ' responses but got '+str(len(i))+ ' !')
						else:
							for j in i:
								if not(isinstance(j, int)):
									raise Exception('Error: Invalid `responses` parameter (nested list of int)!')
			elif len(args) == 1:
				if not(DiscreteMultivariateConditionalData.__instancecheck__(args[0]) or _DiscreteMultivariateConditionalData.__instancecheck__(args[0])):
					raise Exception('Error: Invalid `dmcd` parameter !')
			else:
				raise Exception('Error: Invalid initialization !')
		except Exception as e:
			print e
		else:
			_DiscreteMultivariateConditionalData.__init__(self, *args)

	def Display(self):
		"""
			Print conditional histogram of DiscreteMultivariateConditionalData

				:usage:
					e.Display()
		"""
		print _DiscreteMultivariateConditionalData._Display(self)

	def DistributionEstimation(self):
		"""
			Estimate a DiscreteMultivariateConditionalDistribution

				:usage:
					d = e.DistributionEstimation()
		"""
		return DiscreteMultivariateConditionalDistribution(_DiscreteMultivariateConditionalData._DistributionEstimation(self))

	def ParametricEstimation(self, **kwargs):
		"""
			Estimate a DiscreteMultivariateConditionalDistribution

				:usage:
					d = e.ParametricEstimation()
		"""
		return DiscreteMultivariateConditionalParametric(_DiscreteMultivariateConditionalData._TypeParametricEstimation(self, int(10e2)))

	def KullbackDistance(self, dmcd):
		"""
			return the Kullback distance between DiscreteMultivariateConditionalData and DiscreteMultivariateConditionalDistribution or DiscreteMultivariateConditionalParametric

				:parameter:
					`dmcd` (DiscreteMultivariateConditionalDistribution or DiscreteMultivariateConditionalParametric) - Distribution to use
					
				:usage:
					e.KullbackDistance(d)
		"""
		try:
			if not(DiscreteMultivariateConditionalDistribution.__instancecheck__(dmcd)):
				raise Exception('Error: Invalid `dmcd` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalData._KullbackDistance(self, dmcd)

	def LogLikelihood(self, dmcd):
		"""
			Return log-likelihood of DiscreteMultivariateConditionalData with considered DiscreteMultivariateConditionalDistribution (`dmcd`)

				:parameter:
					`dmcd` (DiscreteMultivariateConditionalDistribution or DiscreteMultivariateConditionalParameteric) - Discrete Multivariate Conditional Distribution used to compute log-likelihood of data

				:usage:
					e.LogLikelihood(d)

				:seealso:
					GetMass
		"""
		try:
			if not(DiscreteMultivariateConditionalDistribution.__instancecheck__(dmcd)):
				raise Exception('Error: Invalid `dmcd` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteMultivariateConditionalData._Likelihood(self, dmcd, True)
