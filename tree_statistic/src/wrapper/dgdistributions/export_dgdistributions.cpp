#include <boost/shared_ptr.hpp>
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

using namespace stat_tool;
using namespace boost::python;

#include "tree_statistic/graphical_distribution.h"

#define WRAP DiscreteGraphicalDistributionWrap

class WRAP
{
	public:
		/*static boost::shared_ptr< DiscreteGraphicalDistribution > DGD_wrapper_init(str graphstring, const DiscreteMultivariateHistogram& hist)
		{
			typedef boost::property < boost::vertex_name_t, int> VertexProperty;
			typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, VertexProperty> GraphFromFile;
			GraphFromFile g;
			boost::dynamic_properties dp;
			boost::property_map<GraphFromFile, boost::vertex_name_t>::type name = get(boost::vertex_name, g);
			dp.property("node_id",name);
			boost::read_graphviz(extract<std::string>(graphstring), g, dp);
			UndirectedGraph ugraph(boost::num_vertices(g));
			UndirectedGraphTraits::edge_iterator ei, ei_end;
			for(boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
				boost::add_edge(boost::source(*ei,g), boost::target(*ei,g), ugraph);
			DiscreteGraphicalDistribution *dgd = new DiscreteGraphicalDistribution(ugraph, hist);
			return boost::shared_ptr< DiscreteGraphicalDistribution >(dgd);
		}*/
		static boost::python::list DGD_wrapper_compare(const DiscreteGraphicalDistribution& dgd, str graphstring)
		{
			typedef boost::property < boost::vertex_name_t, std::string> VertexProperty;
			typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, VertexProperty> GraphFromString;
			GraphFromString g;
			boost::dynamic_properties dp;
			boost::property_map<GraphFromString, boost::vertex_name_t>::type name = get(boost::vertex_name, g);
			dp.property("node_id",name); 
			boost::read_graphviz(extract<std::string>(graphstring), g, dp);
			UndirectedGraph rgraph(boost::num_vertices(g));
			UndirectedGraphTraits::edge_iterator ei, ei_end;
			for(boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
				boost::add_edge(boost::source(*ei,g), boost::target(*ei,g), rgraph);
			std::pair<double, double> comparison = dgd.compare(rgraph);
			boost::python::list l_comparison;
			l_comparison.append(comparison.first);
			l_comparison.append(comparison.second);
			return l_comparison;
		}
		static float DGD_wrapper_get_mass(DiscreteGraphicalDistribution& dgd, boost::python::list l_event)
		{
			std::vector<int> event(len(l_event));
			for(unsigned int i = 0; i < len(l_event); ++i)
				event[i] = extract<int>(l_event[i]);
			return dgd.get_mass(event,false);
		}
		static int DGD_wrapper_get_nb_variables(const DiscreteGraphicalDistribution& dgd)
		{
			return dgd.get_nb_variables();
		}
		static int DGD_wrapper_get_nb_edges(const DiscreteGraphicalDistribution& dgd)
		{
			return dgd.get_nb_edges();
		}
		static void DGD_wrapper_get_graph(const DiscreteGraphicalDistribution& dgd, str filename)
		{
			UndirectedGraph ugraph = dgd.get_ugraph();
			std::ofstream file;
			file.open(extract<char*>(filename));
			boost::write_graphviz(file, ugraph);
			file.close();
		}
		static int DGD_wrapper_get_nb_parameters(const DiscreteGraphicalDistribution& dgd)
		{
			return dgd.get_nb_parameters();
		}
		static str DGD_wrapper_get_formula(const DiscreteGraphicalDistribution& dgd)
		{
			return str(dgd.get_formula());
		}
};

void class_discrete_graphical_distribution()
{
	class_< DiscreteGraphicalDistribution >
		("_DiscreteGraphicalDistribution", "Discrete Graphical Distribution", init< const DiscreteGraphicalDistribution& >())
		//.def("__init__", make_constructor(WRAP::DGD_wrapper_init))
		.def("_Compare", WRAP::DGD_wrapper_compare, "_Compare(self, str) -> list"
			"Return comparison between str graph and DiscreteGraphicalDistribution")
		.def("_GetMass", WRAP::DGD_wrapper_get_mass, "_GetMass(self, list) -> float"
			"Return mass of event")
		.def("_GetNbVariables", WRAP::DGD_wrapper_get_nb_variables, "_GetNbVariables(self) -> int"
			"Return the number of vertices of DiscreteGRaphicalDistribution")
		.def("_GetNbEdges", WRAP::DGD_wrapper_get_nb_edges, "_GetNbEdges(self) -> int"
			"Return the number of edges of DiscreteGRaphicalDistribution")
		.def("_GetGraph", WRAP::DGD_wrapper_get_graph, "_GetGraph(self, str)"
			"Write DiscreteGraphicalDistribution graph in str file")
		.def("_GetNbParameters", WRAP::DGD_wrapper_get_nb_parameters, "_GetNbParameters(self) -> int"
			"Return number of parameters for DiscreteGraphicalDistribution")
		.def("_GetFormula", WRAP::DGD_wrapper_get_formula, "_GetFormula(self) ->str"
			"Return decomposition formula of DiscreteGraphicaldisrtibution")
	;
};

#undef WRAP

#define WRAP DiscreteGraphicalParametricWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteGraphicalParametric > DGP_wrapper_init(str graphstring, boost::python::list distributions)
		{
			typedef boost::property < boost::vertex_name_t, int> VertexProperty;
			typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperty> GraphFromFile;
			GraphFromFile g;
			boost::dynamic_properties dp;
			boost::property_map<GraphFromFile, boost::vertex_name_t>::type name = get(boost::vertex_name, g);
			dp.property("node_id",name);
			boost::read_graphviz(extract<std::string>(graphstring), g, dp);
			BidirectedGraph bgraph(boost::num_vertices(g));
			BidirectedGraphTraits::edge_iterator ei, ei_end;
			for(boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
				boost::add_edge(boost::source(*ei,g), boost::target(*ei,g), bgraph);
			DistributionCollection distribution_collection;
			Neighborhood covariables, variables;
			for(unsigned int i = 0; i < len(distributions); ++i)
			{
				covariables.clear();
				variables.clear();
				for(unsigned int j = 0; j < len(distributions[i][0]); ++j)
					covariables.insert(extract<unsigned int>(distributions[i][0][j]));
				for(unsigned int j = 0; j < len(distributions[i][1]); ++j)
					covariables.insert(extract<unsigned int>(distributions[i][1][j]));
				distribution_collection.decomposition.push_back(std::pair<Neighborhood, Neighborhood>(covariables, variables));
				if(covariables.size() == 0)
				{
					if(variables.size() == 1)
						(distribution_collection.DUD).push_back(extract< DiscreteParametric* >(distributions[i][2]));
					else
						(distribution_collection.DMD).push_back(extract< Stat_trees::DiscreteMultivariateParametric* >(distributions[i][2]));
				} else {
					if(variables.size() == 1)
						(distribution_collection.DUCD).push_back(extract< DiscreteUnivariateConditionalParametric* >(distributions[i][2]));
					else
						(distribution_collection.DMCD).push_back(extract< DiscreteMultivariateConditionalParametric* >(distributions[i][2]));
				}
			}
			DiscreteGraphicalParametric *dgp = new DiscreteGraphicalParametric(boost::num_vertices(bgraph), bgraph, distribution_collection);
			return boost::shared_ptr< DiscreteGraphicalParametric >(dgp);
		}
};

void class_discrete_graphical_parametric()
{
	class_< DiscreteGraphicalParametric, bases< DiscreteGraphicalDistribution > >
			("_DiscreteGraphicalParametric", "Discrete Graphical Parametric Distributions", init< const DiscreteGraphicalParametric& >())
		.def("__init__", make_constructor(WRAP::DGP_wrapper_init))
	;
};

#undef WRAP

template<int num> struct DGEAlgorithms { int v; enum { value = num };
DGEAlgorithms(int _v) : v(_v) { } operator int() const { return v; } };

#define WRAP DGEEnumAWrap

void enum_dge_algorithms_ids()
{
	enum_<DGEAlgorithms<5> >("DGEAlgorithms")
		.value("Relevance", DiscreteGraphicalData::RELEVANCE)
		.value("CLT",DiscreteGraphicalData::CLT)
		.value("CLR", DiscreteGraphicalData::CLR)
		.value("MRNET", DiscreteGraphicalData::MRNET)
		.value("ML", DiscreteGraphicalData::ML)
		.export_values()
	;
};

#undef WRAP

template<int num> struct DGECriterions { int v; enum { value = num };
DGECriterions(int _v) : v(_v) { } operator int() const { return v; } };

#define WRAP DGEEnumCWrap

void enum_dge_criterions_ids()
{
	enum_<DGECriterions<2> >("DGECriterions")
		.value("AIC", DiscreteGraphicalData::AIC)
		.value("BIC",DiscreteGraphicalData::BIC)
		.export_values()
	;
};
#undef WRAP

#define WRAP DiscreteGraphicalDataWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteGraphicalData > DGDa_wrapper_init(int nb_variables, boost::python::list l_event)
		{
			std::vector< std::vector<int> > event(len(l_event), std::vector<int>(nb_variables,0));
			for(unsigned int i = 0; i < len(l_event); ++i)
			{
				for(unsigned int j = 0; j < len(l_event[i]); ++j)
					event[i][j] = extract<int>(l_event[i][j]);
			}
			DiscreteGraphicalData* dgd = new DiscreteGraphicalData(nb_variables, event);
			return boost::shared_ptr< DiscreteGraphicalData >(dgd);
		}
		static  DiscreteGraphicalDistribution* DGDa_wrapper_distribution_estimation0(const DiscreteGraphicalData& dge, int algo, int max_step)
		{
			DiscreteGraphicalDistribution *dgd = dge.distribution_estimation(algo, max_step);
			return dgd;
		}
		static DiscreteGraphicalDistribution* DGDa_wrapper_distribution_estimation1(const DiscreteGraphicalData& dge, int algo, double threshold)
		{
			DiscreteGraphicalDistribution *dgd = dge.distribution_estimation(algo, threshold);
			return dgd;
		}
		static  DiscreteGraphicalParametric* DGDa_wrapper_parametric_estimation0(const DiscreteGraphicalData& dge, int algo, int max_step)
		{
			DiscreteGraphicalParametric *dgd = dge.type_parametric_estimation(algo, max_step);
			return dgd;
		}
		static DiscreteGraphicalParametric* DGDa_wrapper_parametric_estimation1(const DiscreteGraphicalData& dge, int algo, double threshold)
		{
			DiscreteGraphicalParametric *dgd = dge.type_parametric_estimation(algo, threshold);
			return dgd;
		}
		static boost::python::list DGDa_wrapper_stepwise_distribution_estimation(const DiscreteGraphicalData& dge, int algo)
		{
			std::vector< DiscreteGraphicalDistribution* > dgds = dge.stepwise_distribution_estimation(algo);
			boost::python::list l_dgds;
			for(std::vector< DiscreteGraphicalDistribution* >::iterator it = dgds.begin(); it != dgds.end(); ++it)
				l_dgds.append(*it);
			return l_dgds;
		}
		/*static list DGE_wrapper_ROC(const DiscreteGraphicalEstimation& dge, str graphstring, int algo)
		{
			typedef boost::property < boost::vertex_name_t, int> VertexProperty;
			typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, VertexProperty> GraphFromFile;
			GraphFromFile g;
			boost::dynamic_properties dp;
			boost::property_map<GraphFromFile, boost::vertex_name_t>::type name = get(boost::vertex_name, g);
			dp.property("node_id",name);
			boost::read_graphviz(extract<std::string>(graphstring), g, dp);
			UndirectedGraph rgraph(boost::num_vertices(g));
			UndirectedGraphTraits::edge_iterator ei, ei_end;
			for(boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
				boost::add_edge(boost::source(*ei,g), boost::target(*ei,g), rgraph);
			std::vector< DiscreteGraphicalDistribution* > dgds = dge.StepwiseDistributionEstimation(algo);
			std::pair<double, double> point;
			list x;
			list y;
			list l_points;
			for(std::vector< DiscreteGraphicalDistribution* >::iterator it = dgds.begin(); it != dgds.end(); ++it)
			{
				point = (*it)->Compare(rgraph);
				x.append(point.first);
				y.append(point.second);
			}
			l_points.append(x);
			l_points.append(y);
			return l_points;
		}*/

		static float DGDa_wrapper_loglikelihood(const DiscreteGraphicalData& dge, DiscreteGraphicalDistribution& dgd)
		{
			return dge.likelihood_computation(dgd);
		}

		static float DGDa_wrapper_penalized_loglikelihood(const DiscreteGraphicalData& dge, DiscreteGraphicalDistribution& dgd, int criterion)
		{
			return dge.penalized_likelihood_computation(dgd, criterion);
		}

		static float DGDa_wrapper_kullback_distance(const DiscreteGraphicalData& dge, DiscreteGraphicalDistribution& dgd)
		{
			return dge.kullback_distance_computation(dgd);
		}
		static str  DGDa_wrapper_display(const DiscreteGraphicalData& dge)
		{
			std::stringstream s;
      str res;
      dge.display(s);
      res = str(s.str());
      return res;
		}
		/*static list DGE_wrapper_BP(const DiscreteGraphicalEstimation& dge, int algo)
		{
			DiscreteMultivariateHistogram *hist = dge.GetHistogram();
			std::vector< DiscreteGraphicalDistribution* > dgds = dge.StepwiseDistributionEstimation(algo);
			list x;
			list y;
			list l_points;
			for(std::vector< DiscreteGraphicalDistribution* >::iterator it = dgds.begin(); it != dgds.end(); ++it)
			{
				if((*it)->TriangulatedGraphChecking())
				{
					y.append(2*(*it)->ComputeLogLikelihood(*hist)-log(hist->GetNbEvents())*(*it)->GetNbParameters());
					x.append(distance(dgds.begin(), it));
				}
			}
			l_points.append(x);
			l_points.append(y);
			return l_points;
		}*/
};

void class_discrete_graphical_data()
{
	class_< DiscreteGraphicalData >
		("_DiscreteGraphicalData", "Discrete Graphical Data", init< const DiscreteGraphicalData& >())
		.def("__init__", make_constructor(WRAP::DGDa_wrapper_init))
		.def("_DistributionEstimation", WRAP::DGDa_wrapper_distribution_estimation0, return_value_policy< manage_new_object >(), "_DistributionEstimation(self, int) -> _DiscreteGraphicalDistribution"
			"_DiscreteGraphicalDistribution estimation")
		.def("_DistributionEstimationThreshold", WRAP::DGDa_wrapper_distribution_estimation1, return_value_policy< manage_new_object >(), "_DistributionEstimationThreshold(self, int, double) -> _DiscreteGraphicalDistribution"
			"_DiscreteGraphicalDistribution estimation by thresholds algorithms")
		.def("_StepwiseDistributionEstimation", WRAP::DGDa_wrapper_stepwise_distribution_estimation, "_StepwiseDistributionEstimation(self, int) -> list"
			"Return list of estimated DiscreteGraphicalDistribution with int algorithm")
		.def("_ParametricEstimation", WRAP::DGDa_wrapper_parametric_estimation0, return_value_policy< manage_new_object >(), "_ParametricEstimation(self, int) -> _DiscreteGraphicalParametric"
			"_DiscreteGraphicalParametric estimation")
		.def("_ParametricEstimationThreshold", WRAP::DGDa_wrapper_parametric_estimation1, return_value_policy< manage_new_object >(), "_ParametricEstimationThreshold(self, int, double) -> _DiscreteGraphicalParametric"
			"_DiscreteGraphicalParametric estimation by thresholds algorithms")
		.def("_LogLikelihood", WRAP::DGDa_wrapper_loglikelihood, "_LogLikelihood(self, _DiscreteGraphicalDistribution) ->float"
			"Return LogLikelihood of data considering _DiscreteGraphicalDistribution")
		.def("_PenalizedLogLikelihood", WRAP::DGDa_wrapper_penalized_loglikelihood, "_PenalizedLogLikelihood(self, _DiscreteGraphicalDistribution, int) ->float"
			"Return Penalized LogLikelihood of data considering _DiscreteGraphicalDistribution with int criterion")
		.def("_KullbackDistance", WRAP::DGDa_wrapper_kullback_distance, "_KullbackDistance(self, _DiscreteGraphicalDistribution) ->float"
			"Return Kullback distance between DiscreteGraphicalData and _DiscreteGraphicalDistribution")
		.def("_Display", WRAP::DGDa_wrapper_display, "_Display(self) -> str"
			"Return histogram as a string")
		//.def("_ROC", WRAP::DGE_wrapper_ROC, "_ROC(self, _DiscreteGraphicalEstimation, int) -> list"
		//	"Return plotable ROC curves for DiscreteGraphicalDistribution estimation with int algorithm")
	;
};

#undef WRAP

