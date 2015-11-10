import matplotlib.pyplot
#from openalea.stat_tool.distribution import DistributionIdentifierType as f
from _dgdistributions import _DiscreteGraphicalDistribution
from _dgdistributions import DGEAlgorithms as Algorithms
from _dgdistributions import DGECriterions as Criterions
from _dgdistributions import _DiscreteGraphicalParametric
#from _dmcdistributions import _DiscreteUnivariateConditionalParametric
from _dgdistributions import _DiscreteGraphicalData
import networkx
import tempfile
import math

class DiscreteGraphicalDistribution(_DiscreteGraphicalDistribution):
	"""
		A class for Discrete Graphical Distributions in python

	"""

	def __init__(self, *args):
		"""
			Initialize a DiscreteGraphicalDistribution object

				:parameter:
					`dgd` (DiscreteGraphicalDistribution) - A Discrete Graphical Distribution to copy

				:usage:
					distribution_bis = DiscreteGraphicalDistribution(distribution)
		"""
		try:
			if len(args) == 1:
				if not(DiscreteGraphicalDistribution.__instancecheck__(args[0]) or _DiscreteGraphicalDistribution.__instancecheck__(args[0])):
					raise Exception('Error: Invalid `dgd` parameter !')
			else:
				raise Exception('Error: Invalid initialization for DiscreteGraphicalDistribution !')
		except Exception as e:
			print e
		else:
			_DiscreteGraphicalDistribution.__init__(self, *args)

	def GetNbVariables(self):
		"""
			Return number of vertices in DiscreteGraphicalDistribution

				:usage:
					distribution.GetNbVariables()
		"""
		return _DiscreteGraphicalDistribution._GetNbVariables(self)

	def GetNbEdges(self):
		"""
			Return number of edges in DiscreteGraphicalDistribution graph

				:usage:
					distribution.GetNbEdges()
		"""
		return _DiscreteGraphicalDistribution._GetNbEdges(self)

	def GetNbParameters(self):
		"""
			Return number of parameters for DiscreteGraphicalDistribution

				:usage:
					distribution.GetNbParameters()
		"""
		return _DiscreteGraphicalDistribution._GetNbParameters(self)

	def GetMass(self, event):
		"""
			Return mass of event

				:parameter:
					event (list) - Event for which mass needs to be computed

				:usage:
					distribution.GetMass(event)
		"""
		try:
			if not(isinstance(event, list)):
				raise Exception('Error: Invalid `event` parameter !')
			elif not(len(event) == self.GetNbVariables()):
				raise exception('Error: `event` parameter should have length of '+ str(self.GerNbVariables()) + ' but got '+ str(len(event)))
			else:
				for i in event:
					if not(isinstance(i ,int)):
						raise Exception('Error: `'+str(i)+'` should be an integer !')
		except Exception as e:
			print e
		else:
			_DiscreteGraphicalDistribution._GetMass(self, event)

	def Compare(self, graph):
		"""
			Compare graph with DiscreteGraphicalDistribution graph

				:parameter:
					`ugraph` (networkx.Graph) : A undirected reference graph

				:usage:
					d.Compare(G)
		"""
		try:
			if isinstance(graph, networkx.Graph):
				if graph.is_directed():
					raise Exception('Error: `ugraph` should be an undirected graph !')
				else:
					f = tempfile.NamedTemporaryFile(delete=False)
					filename = f.name
					f.close()
					networkx.write_dot(graph, filename)
					ugraph = open(filename,'r').read()
					f.unlink;		
			else:
				raise Exception('Error: `ugraph` parameter not valid !')
		except Exception as e:
			print e
		else:
			return _DiscreteGraphicalDistribution._Compare(self, ugraph)

	def GetGraph(self):
		"""
			Return DiscreteGraphicalDistribution Graph

				:usage:
					G = distribution.GetGraph()
		"""
		f = tempfile.NamedTemporaryFile(delete=False)
		f.close()
		_DiscreteGraphicalDistribution._GetGraph(self, f.name)
		G = networkx.read_dot(f.name)
		f.unlink(f.name)
		return G

	def Graph(self):
		"""
			Plot DiscreteGraphicalDiscreteGraph

				:usage:
					distribution.Graph()
					matplotlib.pyplot.show()
		"""
		G = self.GetGraph()
		try:
			pos = networkx.graphviz_layout(G)
		except:
			pos = networkx.spring_layout(G,iterations=20)
		networkx.draw_networkx_edges(G, pos)
		networkx.draw_networkx_nodes(G, pos, node_color = 'w')
		networkx.draw_networkx_labels(G, pos)
		matplotlib.pyplot.title(_DiscreteGraphicalDistribution._GetFormula(self))
		matplotlib.pyplot.axis('off')

	def Formula(self):
		"""
			Print DiscreteGraphicalDistribution decomposition formula

				:usage:
					distribution.Decomposition()
		"""
		print(_DiscreteGraphicalDistribution._GetFormula(self))


class DiscreteGraphicalParametric(_DiscreteGraphicalParametric, DiscreteGraphicalDistribution):
	"""
		A class for Discrete Graphical Parametric Distributions in python

	"""
	def __init__(self, *args):
		try:
			if(len(args) != 1):
				raise Exception('Error: Invalid initialization')
			elif not(_DiscreteGraphicalParametric.__instancecheck__(args[0])):
				raise Exception('Error: Invalid `dgp` parameter !')
		except Exception as e:
			print e
		else:
			_DiscreteGraphicalParametric.__init__(self, *args)

class DiscreteGraphicalData(_DiscreteGraphicalData):
	"""
		A class for Discrete Graphical Data in python
	"""

	def __init__(self, *args):
		"""
			Initialize a Discrete Graphical Data object

				:parameter:
					`dgd` (DiscreteGraphicalData) - A DicreteGraphicalData to copy

				:parameter:
					`nb_vertices` (int) - Number of vertices for DiscreteGraphicalData
					`events` (nested list of int) - Events to consider in DiscreteGraphicalData

				:usage:
					events = [[0,0],[0,1],[0,1],[0,0],[1,1],[1,0],[1,0],[0,1],[1,0]]
					data = DiscreteGraphicalData(2, events)
					data_bis = DiscreteGraphicalData(data)
		"""
		try:
			if len(args) == 2:
				if not(isinstance(args[0], int)):
					raise Exception('Error: Invalid `nb_vertices` (int) parameter !')
				elif not(isinstance(args[1], list)):
					raise Exception('Error: Invalid `events` (nested list of int) parameter !')
				else:
					for i in args[1]:
						if not(len(i) == args[0]):
							raise Exception('Error: nested list in events should be of length '+ str(args[0]) +' but got '+ str(len(i))+ ' !')
						for j in i:
							if not(isinstance(j, int)):
								raise Exception('Error: `'+ str(j) +'` is not an integer !')
			elif len(args) == 1:
				if not(_DiscreteGraphicalData.__instancecheck__(args[0]) or DiscreteGraphicalData.__instancecheck__(args[0])):
					raise Exception('Error: Invalid `dgd` (DiscreteGraphicalData) parameter !')
			else:
				raise Exception('Error: Invalid initialization for DiscreteGraphicalData !')
		except Exception as e:
			print e
		else:
			_DiscreteGraphicalData.__init__(self, *args)

	def Display(self):
		"""
			Print histogram

		"""
		print(_DiscreteGraphicalData._Display(self))

	def DistributionEstimation(self, algorithm, **kwargs):
		"""
			DiscreteGraphicalDistribution estimation from DiscreteGraphicalData

				:parameter:
					`algorithm` (Algorithms) - Algorithm to use for estimation
					`threshold` (float) - Threshold to use with algorithm

				:parameter:
					`algorithm` (Algorithms) - Algorithm to use for estimation
					`number_of_edges` (int) - Number of edges for final graph choosen with algorithm

				:parameter:
					`algorithm` (Algorithms) - Algorithm to use for estimation
					`all` (bool) - Return all distribution infered using algorithm

				:usage:

					data.DistributionEstimation(Algorithms.CLT)
					data.DistributionEstimation(Algorithms.CLT, number_of_edges = 2)
					data.DistributionEstimation(Algorithms.CLT, all = True)
		"""
		try:
			if not(Algorithms.__instancecheck__(algorithm)):
				raise Exception('Error: Invalid `algorithm` (Algorithms) parameter !')
			elif 'all' in kwargs:
				if not(isinstance(kwargs['all'], bool)):
					raise Exception('Error: Invalid `all` (bool) parameter !')
			elif 'threshold' in kwargs:
				if not(isinstance(kwargs['threshold'], float)):
					raise Exception('Error: Invalid `threshold` (float) parameter !')
				else:
					threshold = kwargs['threshold']
			elif 'number_of_edges' in kwargs:
				if not(isinstance(kwargs['number_of_edges'], int)):
					raise Exception('Error: invalid `number_of_edges` (int) parameter !')
				else:
					number_of_edges = kwargs['number_of_edges']
		except Exception as e:
			print e
		else:
			if 'all' in kwargs:
				dgds = _DiscreteGraphicalData._StepwiseDistributionEstimation(self, algorithm)
				for i in range(0, len(dgds)):
					dgds[i] = DiscreteGraphicalDistribution(dgds[i])
				return dgds
			elif 'threshold' in kwargs:
				return DiscreteGraphicalDistribution(_DiscreteGraphicalData._DistributionEstimationThreshold(self, algorithm, threshold))
			elif 'number_of_edges' in kwargs:
				return DiscreteGraphicalDistribution(_DiscreteGraphicalData._DistributionEstimation(self, algorithm, number_of_edges))
			else:
				return DiscreteGraphicalDistribution(_DiscreteGraphicalData._DistributionEstimation(self, algorithm, -1))

	def ParametricEstimation(self, algorithm, **kwargs):
		"""
			DiscreteGraphicalParametric estimation from DiscreteGraphicalData

				:parameter:
					`algorithm` (Algorithms) - Algorithm to use for estimation
					`threshold` (float) - Threshold to use with algorithm

				:parameter:
					`algorithm` (Algorithms) - Algorithm to use for estimation
					`number_of_edges` (int) - Number of edges for final graph choosen with algorithm

				:parameter:
					`algorithm` (Algorithms) - Algorithm to use for estimation
					`all` (bool) - Return all distribution infered using algorithm

				:usage:

					data.ParametricEstimation(Algorithms.CLT)
					data.ParametricEstimation(Algorithms.CLT, number_of_edges = 2)
					data.ParametricEstimation(Algorithms.CLT, all = True)
		"""
		try:
			if not(Algorithms.__instancecheck__(algorithm)):
				raise Exception('Error: Invalid `algorithm` (Algorithms) parameter !')
			elif 'all' in kwargs:
				if not(isinstance(kwargs['all'], bool)):
					raise Exception('Error: Invalid `all` (bool) parameter !')
			elif 'threshold' in kwargs:
				if not(isinstance(kwargs['threshold'], float)):
					raise Exception('Error: Invalid `threshold` (float) parameter !')
				else:
					threshold = kwargs['threshold']
			elif 'number_of_edges' in kwargs:
				if not(isinstance(kwargs['number_of_edges'], int)):
					raise Exception('Error: invalid `number_of_edges` (int) parameter !')
				else:
					number_of_edges = kwargs['number_of_edges']
		except Exception as e:
			print e
		else:
			if 'all' in kwargs:
				dgds = list()#_DiscreteGraphicalData._StepwiseDistributionEstimation(self, algorithm)
				for i in range(0, len(dgds)):
					dgds[i] = DiscreteGraphicalParametric(dgds[i])
				return dgds
			elif 'threshold' in kwargs:
				return DiscreteGraphicalParametric(_DiscreteGraphicalData._ParametricEstimationThreshold(self, algorithm, threshold))
			elif 'number_of_edges' in kwargs:
				return DiscreteGraphicalParametric(_DiscreteGraphicalData._ParametricEstimation(self, algorithm, number_of_edges))
			else:
				return DiscreteGraphicalParametric(_DiscreteGraphicalData._ParametricEstimation(self, algorithm, -1))

	def GetLogLikelihood(self, dgd):
		"""
			get log-likelihood of DiscreteGraphicalDistribution(s)

				:parameter:
					`dgd` (DiscreteGraphicalDistribution or list) - DiscreteGraphicalDistribuion(s) from which Log-Likelihood of DiscreteGraphicalData is computed

				:usage:
					data.GetLogLikelihood(dgd)
		"""
		try:
			if not(DiscreteGraphicalDistribution.__instancecheck__(dgd) or isinstance(dgd, list)):
				raise Exception('Error: Invalid `dgd` parameter !')
			elif isinstance(dgd, list):
				for i in dgd:
					if not(DiscreteGraphicalDistribution.__instancecheck__(i)):
						raise Exception('Error: Invalid `dgd` parameter !')
		except Exception as e:
			print e
		else:
			if isinstance(dgd, list):
				loglikelihood = list()
				for i in dgd:
					loglikelihood.append(_DiscreteGraphicalData._LogLikelihood(self, i))
			else:
				loglikelihood = _DiscreteGraphicalData._LogLikelihood(self, dgd)
			return loglikelihood

	def GetPenalizedLogLikelihood(self, dgd, criterion):
		"""
			get penalized log-likelihood of DiscreteGraphicalDistribution(s)

				:parameter:
					`dgd` (DiscreteGraphicalDistribution or list) - DiscreteGraphicalDistribuion(s) from which Log-Likelihood of DiscreteGraphicalData is computed
					`criterion` (Criterions) - Criterion to use for penalized log-likelihood

				:usage:
					data.GetPenalizedLogLikelihood(dgd, Criterions.AIC)
		"""
		try:
			if not(DiscreteGraphicalDistribution.__instancecheck__(dgd) or isinstance(dgd, list)):
				raise Exception('Error: Invalid `dgd` parameter !')
			elif isinstance(dgd, list):
				for i in dgd:
					if not(DiscreteGraphicalDistribution.__instancecheck__(i)):
						raise Exception('Error: Invalid `dgd` parameter !')
			elif not(Criterions.__instancecheck__(criterion)):
				raise Exception('Error: Invalid `criterion` parameter !')
		except Exception as e:
			print e
		else:
			if isinstance(dgd, list):
				ploglikelihood = list()
				for i in dgd:
					ploglikelihood.append(_DiscreteGraphicalData._PenalizedLogLikelihood(self, i, criterion))
			else:
				ploglikelihood = _DiscreteGraphicalData._LogLikelihood(self, dgd, criterion)
			return ploglikelihood

	def LogLikelihood(self, dgd, **kwargs):
		"""
			Plot or get LogLikelihood of DiscreteGraphicalData using DiscreteGraphicalDistribution(s)

				:parameter:
					`dgd` (DiscreteGraphicalDistribution or list) - DiscreteGraphicalDistribuion(s) from which Log-Likelihood of DiscreteGraphicalData is computed

				:usage:
					data.LogLikelihood(dgd)
		"""
		try:
			if not(DiscreteGraphicalDistribution.__instancecheck__(dgd) or isinstance(dgd, list)):
				raise Exception('Error: Invalid `dgd` parameter !')
			elif isinstance(dgd, list):
				for i in dgd:
					if not(DiscreteGraphicalDistribution.__instancecheck__(i)):
						raise Exception('Error: Invalid `dgd` parameter !')
		except Exception as e:
			print e
		else:
			if isinstance(dgd, list):
				loglikelihood = self.GetLogLikelihood(dgd)
				x = list()
				for i in dgd:
					x.append(i.GetNbEdges())
				matplotlib.pyplot.plot(x, loglikelihood, 'o-', **kwargs)
				matplotlib.pyplot.xlabel('Number of edges')
				matplotlib.pyplot.ylabel('Log-likelihood')

			else:
				return self.GetLogLikelihood(dgd)		

	def PenalizedLogLikelihood(self, dgd, criterion, **kwargs):
		"""
			Plot or get Penalized LogLikelihood of DiscreteGraphicalData using DiscreteGraphicalDistribution(s)

				:parameter:
					`dgd` (DiscreteGraphicalDistribution or list) - DiscreteGraphicalDistribuion(s) from which Log-Likelihood of DiscreteGraphicalData is computed
					`criterion 

				:usage:
					data.PenalizedLogLikelihood(dgd, Criterions.BIC)
		"""
		try:
			if not(DiscreteGraphicalDistribution.__instancecheck__(dgd) or isinstance(dgd, list)):
				raise Exception('Error: Invalid `dgd` parameter !')
			elif isinstance(dgd, list):
				for i in dgd:
					if not(DiscreteGraphicalDistribution.__instancecheck__(i)):
						raise Exception('Error: Invalid `dgd` parameter !')
			elif not(Criterion.__instancecheck__(criterion)):
				raise Exception('Error: Invalid `criterion` parameter !')
		except Exception as e:
			print e
		else:
			if isinstance(dgd, list):
				ploglikelihood = self.GetPenalizedLogLikelihood(dgd, criterion)
				x = list()
				for i in dgd:
					x.append(i.GetNbEdges())
				matplotlib.pyplot.plot(x, ploglikelihood, 'o-', **kwargs)
				matplotlib.pyplot.xlabel('Number of edges')
				matplotlib.pyplot.ylabel(criterion.name)
			else:
				return self.GetPenalizedLogLikelihood(dgd, criterion)		


def Graphs(dgds):
	"""
		Plot graphs for list of DiscreteGraphicalDistribution

			:parameter:
				`dgds` (list of DiscreteGraphicalDistribution) - list from which graphs needs to be plotted

			:usage:
				Graphs(distributions)
	"""
	try:
		if not(isinstance(dgds, list)):
			raise Exception('Error: Invalid `dgds` parameter !')
		else:
			for i in dgds:
				if not(DiscreteGraphicalDistribution.__instancecheck__(i)):
					raise Exception('Error: Invalid `dgds` parameter !')
	except Exception as e:
		print e
	else:
		cols = math.ceil(math.sqrt(len(dgds)))
		rows = cols
		for i in range(0, len(dgds)):
			matplotlib.pyplot.subplot(rows, cols, i+1)
			matplotlib.pyplot.title('Graph '+str(i))
			dgds[i].Graph()

