import matplotlib.pyplot
from openalea.stat_tool.distribution import DistributionIdentifierType as Family
from _ducdistributions import _DiscreteUnivariateConditionalDistribution
from _ducdistributions import DUCPLinks as Links
from _ducdistributions import _DiscreteUnivariateConditionalParametric
from _ducdistributions import _DiscreteUnivariateConditionalData


class DiscreteUnivariateConditionalDistribution(_DiscreteUnivariateConditionalDistribution):
	"""
		A class for Discrete Univariate Conditional Distributions in python
	"""

	def __init__(self, *args):
		"""
			Initialize a DiscreteUnivariateConditionalDistribution

				:parameter:
					`nb_covariables` (int) : Number of covariables for DiscreteUnivariateConditionalDistribution
					`mass` (list) : A list of covariates and response mass
				
				:parameter:
					`ducd` (DiscreteUnivariateConditionalDistribution) - A DiscreteUnivariateConditionalDistribution to copy

				:usage:
					mass = [[[0,0], [0, 1], [0.25, 0.75]], [[0,1], [0, 1], [0.45, 0.55]], [[1,0], [0, 1], [0.10, 0.90]], [[1,1], [0, 1], [0.05, 0.95]]]
					distribution = DiscreteUnivariateConditionalDistribution(2, mass) # initialize a DiscreteUnivariateConditionalDistribution from masses
					copy = DiscreteUnivariateConditionalDistribution(distribution) # initialize a DiscreteUnivariateCondtionalDistribution from an existing DiscreteUnivariateConditionalDistribution
		"""
		try:
			if len(args) == 2:
				if isinstance(args[0], int):
					if not(args[0] > 0):
						raise Exception('Error: invalid `nb_covariables` parameter !')
				else:
					raise Exception('Error: invalid `nb_covariables` parameter !')
				if not(isinstance(args[1], list)):
					raise Exception('Error: Invalid `mass` parameter !')
				else:
					for i in args[1]:
						if not(len(i[0]) == args[0]):
							raise Exception('Error: Exepecting '+ str(args[0]) + ' covariables but got '+ str(len(i[0])))
						for j in i[0]:
							if not(isinstance(j, int)):
								raise Exception('Error: Invalid `mass` parameter !')
						for j in i[1]:
							if not(isinstance(j, int)):
								raise Exception('Error: Invalid `mass` parameter !')
						if not(len(i[1]) == len(i[2])):
							raise Exception('Error: invalid `mass` parameter !')
						for j in i[2]:
							if not(isinstance(j, float)):
								raise Exception('Error: Invalid `mass` parameter !')
			elif len(args) == 1:
				if not(_DiscreteUnivariateConditionalDistribution.__instancecheck__(args[0]) or DiscreteUnivariateConditionalDistribution.__instancecheck__(args[0])):
					raise Exception('Error: `ducd` parameter not valid !')
			else:
				raise Exception('Error: Invalid parameters for DiscreteGraphicalDistribution initialization !')
		except Exception as e:
			print e
		else:
			_DiscreteUnivariateConditionalDistribution.__init__(self, *args)

	def GetNbCovariables(self):
		"""
			Return number of covariables for DiscreteUnivariateConditionalDistribution

				:usage:
					distribution.GetNbCovariables()
		"""
		return _DiscreteUnivariateConditionalDistribution._GetNbCovariables(self)

	def GetNbParameters(self):
		"""
			Return number of parameters for DiscreteUnivariateConditionalDistribution

				:usage:
					distribution.GetNbParameters()
		"""
		return _DiscreteUnivariateConditionalDistribution._GetNbParameters(self)

	def GetMass(self, covariates, response):
		"""
			Return mass of event for DiscreteUnivariateConditionalDistribution

				:parameter:
					`covariates` (list) : Covariates for event
					`response` (int) : Response variable event

				:usage:
					distribution.GetMass([0,1],1)
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			elif not(len(covariates) == self.GetNbCovariables()):
				raise Exception('Error: Expecting '+ str(self.GetNbCovariables()) + ' covariables but got ' +str(len(covariates))+ ' !')
			elif not(isinstance(response, int)):
				raise Exception('Error: Invalid `response` parameter !')
			else:
				for i in covariates:
					if not(isinstance(i, int)):
						raise Exception('Error: Invalid `covariates` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteUnivariateConditionalDistribution._GetMass(self, covariates, response)

	def Simulate(self, covariates):
		"""
			Simulate a DiscreteUnivariateConditionalDistribution with covariates

				:parameter:
					`covariates` (nested list) - covariates from which response will be simulated.

				:usage:
					distribution.Simulate([[0,1]])
		"""
		try:
			if not(isinstance(covariates, list)):
				raise Exception('Error: Invalid `covariates` parameter !')
			for i in covariates:
				if not(isinstance(i, list)):
					raise Exception('Error: Invalid `covariates` parameter !')
				if not(len(i) == self.GetNbCovariables()):
					raise Exception('Error: Wrong number of covariables !')
				for j in i:
					if not(isinstance(j, int)):
						raise Exception('Error: Invalid `covariates` parameter (nested list of int) !')
		except Exception as e:
			print e
		else:
			return DiscreteUnivariateConditionalData(_DiscreteUnivariateConditionalDistribution._Simulate(self, covariates))

	def GetMean(self, covariates):
		"""
			Return mean with covariates for DiscreteUnivariateConditionalDistribution

				:parameter:
					`covariates` (list) : Covariates

				:usage:
					distribution.GetMean([0,1])
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
			return _DiscreteUnivariateConditionalDistribution._GetMean(self, covariates)

	def GetVariance(self, covariates):
		"""
			Return variance with covariates for DiscreteUnivariateConditionalDistribution

				:parameter:
					`covariates` (list) : Covariates

				:usage:
					distribution.GetVariance([0,1])
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
			return _DiscreteUnivariateConditionalDistribution._GetVariance(self, covariates)


class DiscreteUnivariateConditionalParametric(_DiscreteUnivariateConditionalParametric, DiscreteUnivariateConditionalDistribution):
	"""
		A class for Discrete Univariate Conditional Parametric Distributions in python
	"""

	def __init__(self, *args):
		"""
			Initialize a DiscreteUnivariateConditionalParametric

				:parameter:
					`nb_covariables` (int) - Number of covariables.
					`family` () : An identifier for family distribution.
					`link` (Links) : An identifier for link function.
					`parameters` (list) : A list of parameters.

				:parameter:
					`ducp` (DiscreteUnivariateConditionalParametric) - A DiscreteUnivariateConditionalParametric to copy

				:usage:
					p = DiscreteUnivariateConditionalParametric(stat_tool.Binomial, Links.Logit) # initialize a DiscreteUnivariateConditionalParametric
					pbis = DiscreteUnivariateConditionalParametric(d) # initialize a DiscreteUnivariateCondtionalParametric by copy
		"""
		try:
			if len(args) == 4:
				if not(isinstance(args[0], int)):
					raise Exception('Error: Invalid `nb_covariables` parameter !')
				elif not(Links.__instancecheck__(args[2])):
					raise Exception('Error: Invalid `link` parameter !')
				elif not(isinstance(args[3], list)):
					raise Exception('Error: Invalid `parameters` parameter')
				else:
					for i in args[3]:
						if not(isinstance(i, float)):
							raise Exception('Error: Invalid `parameters` parameter')
			elif len(args) == 1:
				if not(_DiscreteUnivariateConditionalParametric.__instancecheck__(args[0]) or DiscreteUnivariateConditionalParametric.__instancecheck__(args[0])):
					raise Exception('Error: `ducp` parameter not valid !')
			else:
				raise Exception('Error: Invalid parameters for DiscreteUnivariateConditionalParametric initialization !')
		except Exception as e:
			print e
		else:
			_DiscreteUnivariateConditionalParametric.__init__(self, *args)

	def GetFamily(self):
		"""
			Return family of DiscreteUnivariateConditionalParametric

				:usage:
					p.GetFamily()
		"""
		return Family.values[_DiscreteUnivariateConditionalParametric._GetFamily(self)]

	def Family(self):
		return Family.values[_DiscreteUnivariateConditionalParametric._GetFamily(self)].name

	def GetLink(self):
		"""
			Return link function of DiscreteUnivariateConditionalParametric

				:usage:
					p.GetLink()
		"""
		return Links.values[_DiscreteUnivariateConditionalParametric._GetLink(self)]

	def Link(self):
		return Links.values[_DiscreteUnivariateConditionalParametric._GetLink(self)].name
		


class DiscreteUnivariateConditionalData(_DiscreteUnivariateConditionalData):
	"""
		A class for Discrete Univariate Conditional Data estimation in Python
	"""

	def __init__(self, *args):
		"""
			Initialize a DiscreteUnivariateConditionalEstimation

				:parameter:
					`nb_covariables` (int) - Number of covariables
					`covariates` (list) - A list of covariates values
					`response` (list) - Corresponding response values

				:parameter:
					`ducd` (DiscreteUnivariateConditionalData) - A DiscreteUnivariateConditianalData to copy

				:usage:
					c = [[0,1],[0,1],[0,0],[0,1],[1,1]]
					r = [0,1,1,0,0]
					e = DiscreteUnivariateConditionalData(2,c,r)
		"""
		try:
			if len(args) == 3:
				if not(isinstance(args[0], int)):
					raise Exception('Error: Invalid `nb_covariables` parameter !')
				elif not(isinstance(args[1], list) and isinstance(args[2], list)):
					raise Exception('Error: Invalid `covariates` (nested list) or `response` (list) parameter !')
				elif not(len(args[1]) == len(args[2])):
					raise Exception('Error: Invalid `covariates` or `response` parameter (not same number of observations) !')
				else:
					for i in args[1]:
						if not(isinstance(i, list)):
							raise Exception('Error: Invalid `covariates` parameter (nested list)')
						elif not(len(i) == args[0]):
							raise Exception('Error: Expecting '+str(args[0])+ ' covariables but got '+str(len(i))+ ' !')
						else:
							for j in i:
								if not(isinstance(j, int)):
									raise Exception('Error: Invalid `covariates` parameter (nested list of int)!')
					for i in args[2]:
						if not(isinstance(i, int)):
							raise Exception('Error: Invalid `response` parameter (list of int)!')
			elif len(args) == 1:
				if not(DiscreteUnivariateConditionalData.__instancecheck__(args[0]) or _DiscreteUnivariateConditionalData.__instancecheck__(args[0])):
					raise Exception('Error: Invalid `ducd` paramater !')
			else:
				raise Exception('Error: Invalid initialization !')
		except Exception as e:
			print e
		else:
			_DiscreteUnivariateConditionalData.__init__(self, *args)

	def Display(self):
		"""
			Print conditional histogram of DiscreteUnivariateConditionalData

				:usage:
					e.Display()
		"""
		print _DiscreteUnivariateConditionalData._Display(self)

	def DistributionEstimation(self):
		"""
			Estimate a DiscreteUnivariateConditionalDistribution
		"""
		return DiscreteUnivariateConditionalDistribution(_DiscreteUnivariateConditionalData._DistributionEstimation(self))
	
	def ParametricEstimation(self, **kwargs):
		"""
			Estimate a DiscreteUnivariateConditionalParametric with family and link function

				:parameter:
					`family` () - Family of DiscreteUnivariateConditionalParametric.
					`link` (Links) - Link function from esperence of DiscreteUnivariateConditionalParametric and linear regressor.

				:usage:
					e.ParametricEstimation(family = stat_tool.Binomial)
					e.ParametricEstimation(family = stat_tool.Binomial, link = Links.LogLog)
					e.ParametricEstimation(family = stat_tool.Poisson, family = Links.Identity)
					e.ParametricEstimation()
		"""
		try:
			if 'link' in kwargs and 'family' in kwargs:
				if not(Links.__instancecheck__(kwargs['link'])):
					raise Exception('Error: Invalid `link` parameter !')
				elif not(_DiscreteUnivariateConditionalData._CheckCorrespondanceLinkFamily(self, kwargs['family'], kwargs['link'])):
					raise Exception('Error: Not valid `link` parameter for `family` parameter !')
			elif 'family' in kwargs :
				if(not(Family.__instancecheck__(kwargs['family']))):
					raise Exception('Error: Invalid `family` parameter')
				elif kwargs['family'] == Family.BINOMIAL:
					link = Links.Logit
				elif kwargs['family'] == Family.POISSON:
					link = Links.Log
				else:
					raise Exception('Error: Invalid `family` parameter !')
			if 'maxits' in kwargs:
				if not(isinstance(kwargs['maxits'], int)):
					raise Exception('Error: Invalid `maxits` parameter !')
				else:
					maxits = kwargs['maxits']
			else:
				maxits = int(10e2)
		except Exception as e:
			print e
		else:
			if 'link' in kwargs and 'family' in kwargs:
				return DiscreteUnivariateConditionalParametric(_DiscreteUnivariateConditionalData._ParametricEstimation(self, kwargs['family'], kwargs['link'], maxits))
			elif 'family' in kwargs:
				return DiscreteUnivariateConditionalParametric(_DiscreteUnivariateConditionalData._ParametricEstimation(self, kwargs['family'], link, maxits))
			else:
				return DiscreteUnivariateConditionalParametric(_DiscreteUnivariateConditionalData._TypeParametricEstimation(self, maxits))

	def LogLikelihood(self, ducd):
		"""
			Return log-likelihood od DiscreteUnivariateConditionalData with considered DiscreteUnivariateConditionalDistribution (`ducd`)

				:parameter:
					`ducd` (DiscreteUnivariateConditionalDistribution or DiscreteUnivariateConditionalParameteric) - Discrete Univariate Conditional Distribution used to compute log-likelihood of data

				:usage:
					d.LogLikelihood(p)

				:seealso:
					GetMass
		"""
		try:
			if not(DiscreteUnivariateConditionalDistribution.__instancecheck__(ducd) or DiscreteUnivariateConditionalParametric.__instancecheck__(ducd)):
				raise Exception('Error: Invalid `ducd` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteUnivariateConditionalData._Likelihood(self, ducd, True)
			
	def KullbackDistance(self, ducd):
		"""
			Return Kullback Distance between DiscreteUnivariateConditionalData and DiscreteUnivariateConditionalDistribution (`ducd`)

				:parameter:
					`ducd` (DiscreteUnivariateConditionalDistribution or DiscreteUnivariateConditionalParameteric) - Discrete Univariate Conditional Distribution used to compute the Kullback distance

				:usage:
					d.KullbackDistance(p)

				:seealso:
					GetMass
		"""
		try:
			if not(DiscreteUnivariateConditionalDistribution.__instancecheck__(ducd) or DiscreteUnivariateConditionalParametric.__instancecheck__(ducd)):
				raise Exception('Error: Invalid `ducd` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteUnivariateConditionalData._KullbackDistance(self, ducd)

	def GetPearsonResiduals(self, ducd):
		"""
			Return Pearson residuals when fitting DiscreteUnivariateConditionalData with DiscreteUnivariateConditionalDistribution

				:parameter:
					`ducd` (DiscreteUnivariateConditionalDistribution or DiscreteUnivariateConditionalParametric) - Distribution used for fitting DiscreteUnivariateConditionalData

				:usage:
					e.GetPearsonResiduals(d)

				:seealso:
					PearsonResiduals

		"""
		try:
			if not(DiscreteUnivariateConditionalDistribution.__instancecheck__(ducd) or DiscreteUnivariateConditionalParametric.__instancecheck__(ducd)):
				raise Exception('Error: Invalid `ducd` parameter !')
		except Exception as e:
			print e
		else:
			return _DiscreteUnivariateConditionalData._GetPearsonResiduals(self, ducd)

	def PearsonResiduals(self, ducd, **kwargs):
		"""
			Plot Pearson residuals when fitting DiscreteUnivariateConditionalData with DiscreteUnivariateConditionalDistribution

				:parameter:
					`ducd` (DiscreteUnivariateConditionalDistribution or DiscreteUnivariateConditionalParametric) - Distribution used for fitting DiscreteUnivariateConditionalData

				:usage:
					e.GetPearsonResiduals(d)
					plt.legend()
					plt.show()

				:seealso:
					GetPearsonResiduals

		"""
		try:
			if not(DiscreteUnivariateConditionalDistribution.__instancecheck__(ducd) or DiscreteUnivariateConditionalParametric.__instancecheck__(ducd)):
				raise Exception('Error: Invalid `ducd` parameter !')
			elif 'color' in kwargs:
				if not(isinstance(kwargs['color'], str)):
					raise Exception('Error: Invalid `color` parameter !')
				else:
					color = kwargs['color']
			else:
				color = 'b'
		except Exception as e:
			print e
		else:
			pearson_residuals = self.GetPearsonResiduals(ducd)
			if(DiscreteUnivariateConditionalParametric.__instancecheck__(ducd)):
				matplotlib.pyplot.scatter(pearson_residuals[0], pearson_residuals[1], s = pearson_residuals[2], c = color, label = ducd.Family() + ', '+ ducd.Link())
			else:
				matplotlib.pyplot.scatter(pearson_residuals[0], pearson_residuals[1], s = pearson_residuals[2], c = color, label = 'Non-parametric')
			matplotlib.pyplot.xlabel('Predicted mean')
			matplotlib.pyplot.ylabel('Pearson residual')
