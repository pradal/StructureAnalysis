#ifndef DISCRETE_GRAPHICAL_DISTRIBUTION_H
#define DISCRETE_GRAPHICAL_DISTRIBUTION_H

#include "stat_tool/stat_tools.h"
#include "univariate_conditional_distribution.h"
#include "multivariate_distribution.h"
#include "multivariate_conditional_distribution.h"

#ifndef D_INF
#define D_INF std::numeric_limits<double>::infinity()
#endif

#ifndef MAX
#define MAX(x,y)  ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#endif

#ifndef I_MAX
#define I_MAX std::numeric_limits<unsigned int>::max()
#endif

#ifndef I_MIN
#define I_MIN std::numeric_limits<unsigned int>::min()
#endif

#include <boost/math/special_functions/gamma.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/undirected_graph.hpp>

typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS> UndirectedGraph;
typedef boost::graph_traits< UndirectedGraph > UndirectedGraphTraits;

struct estimated_edge {
	unsigned int step;
	double score;
};
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, estimated_edge > EstimationGraph;
typedef boost::graph_traits< EstimationGraph > EstimationGraphTraits;

typedef std::pair<unsigned int, unsigned int> Edge;
typedef std::vector< Edge > Edges;

typedef std::set<unsigned int> Neighborhood;
typedef std::set< Neighborhood > Cliques;

struct clique_vertex
{
	Neighborhood clique;
};
struct separator_edge
{
	Neighborhood separator;
};
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, clique_vertex, separator_edge> CliqueGraph;
typedef boost::graph_traits< CliqueGraph > CliqueGraphTraits;

typedef std::multimap<double, std::pair<Edge, Neighborhood> > EdgesScoresML;
typedef std::pair<double, std::pair<Edge, Neighborhood> > EdgeScoreML;
typedef std::multimap<double, std::pair< Edge, double> > EdgesScoresMRNET;
typedef std::pair<double, std::pair<Edge, double> > EdgeScoreMRNET;
typedef std::multimap<double, Edge > EdgesScores;
typedef std::pair<double, Edge > EdgeScore;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> BidirectedGraph;
typedef boost::graph_traits< BidirectedGraph > BidirectedGraphTraits;

template<typename T>
struct ReestimationCollection
{
	std::vector< std::pair<Neighborhood, Neighborhood> > decomposition;
	std::vector< Reestimation<T>* > DUR;
	std::vector< DiscreteUnivariateConditionalReestimation<T>* > DUCR;
	std::vector< Stat_trees::DiscreteMultivariateReestimation<T>* > DMR;
	std::vector< DiscreteMultivariateConditionalReestimation<T>* > DMCR;
};

struct DistributionCollection
{
	std::vector< std::pair< Neighborhood, Neighborhood > > decomposition;
	std::vector< Distribution* > DUD;
	std::vector< DiscreteUnivariateConditionalDistribution* > DUCD;
	std::vector< Stat_trees::DiscreteMultivariateDistribution* > DMD;
	std::vector< DiscreteMultivariateConditionalDistribution* > DMCD;
};

typedef std::pair<unsigned int, int> Vertice;
typedef std::vector< Vertice > Vertices;

class DiscreteGraphicalDistribution;
class DiscreteGraphicalParametric;

template<typename T>
class DiscreteGraphicalReestimation
{
	public:
		enum
		{
			RELEVANCE,
			CLT,
			CLR,
			MRNET,
			ML
		};
		enum
		{
			AIC,
			BIC
		};
		DiscreteGraphicalReestimation();
		~DiscreteGraphicalReestimation();
		DiscreteGraphicalReestimation(unsigned int nb_variables);
		DiscreteGraphicalReestimation(const DiscreteGraphicalReestimation<T>& dgr);
		void update(const std::vector<int>& event, T frequency);
		void nb_elements_computation();
		std::ostream& display(std::ostream& os) const;
		EstimationGraph graph_estimation(int algo) const;
		EstimationGraph graph_estimation(int algo, int max_step) const;
		EstimationGraph graph_estimation(int algo, double min_threshold) const;
		DiscreteGraphicalDistribution* distribution_estimation(int algo, int max_step=-1) const;
		DiscreteGraphicalDistribution* distribution_estimation(int algo, double min_threshold=-1*D_INF) const;
		std::vector< DiscreteGraphicalDistribution* > stepwise_distribution_estimation(int algo) const;
		//DiscreteGraphicalDistribution* distribution_estimation(const UndirectedGraph& ugraph) const;
		DiscreteGraphicalParametric* type_parametric_estimation(int algo, int max_step=-1) const;
		DiscreteGraphicalParametric* type_parametric_estimation(int algo, double min_threshold=-1*D_INF) const;
		//std::vector< DiscreteGraphicalParametric* > stepwise_type_parametric_estimation(int algo) const;
		//DiscreteGraphicalParametric* type_parametric_estimation(const BidirectedGraph& bgraph) const;
		double likelihood_computation(DiscreteGraphicalDistribution& dgd) const;
		double penalized_likelihood_computation(DiscreteGraphicalDistribution& dgd, int criterion) const;
		double kullback_distance_computation(DiscreteGraphicalDistribution& dgd) const;

	protected:
		unsigned int nb_variables;
		T nb_elements;
		std::map<std::vector<int>, T> histogram;

	private:
		EstimationGraph algorithm_CLT(int max_step) const;
		EstimationGraph algorithms_threshold(int algo, double min_threshold) const;
		EstimationGraph algorithm_ML(int max_step) const;
		EstimationGraph algorithm_MRNET(double min_threshold) const;

		EdgesScoresML initialize_scores_ML() const;
		EdgesScoresMRNET initialize_scores_MRNET() const;
		EdgesScores initialize_scores(int algo) const;

		bool continue_CLT(const EdgesScores& edges_scores, const EstimationGraph& egraph, int max_step, int step) const;
		bool continue_ML(const EdgesScoresML& edges_scores, const EstimationGraph& egraph, int max_step, int step) const;
		bool continue_MRNET(const EdgesScoresMRNET& edges_scores, const EstimationGraph& egraph, double min_threshold) const;
		bool continue_threshold_algorithms(const EdgesScores& edges_scores, const EstimationGraph& egraph, double min_threshold) const;

		void execute_algorithm_CLT(EdgesScores& edges_scores, EstimationGraph& egraph, int& step) const;
		void execute_algorithm_ML(EdgesScoresML& edges_scores, EstimationGraph& egraph, int& step) const;
		void execute_algorithm_MRNET(EdgesScoresMRNET& edges_scores, EstimationGraph& egraph, int& step) const;
		void execute_algorithms(EdgesScores& edges_scores, EstimationGraph& egraph, int& step) const;

		double compute_conditional_mutual_information(const Edge& edge, const Neighborhood& neighborhood) const;
		Neighborhood get_neighborhood(const Edge& edge, const EstimationGraph& egraph) const;
		Edges get_edges() const;

		UndirectedGraph undirected_graph_computation(const EstimationGraph& egraph) const;
		std::vector<UndirectedGraph> undirected_graphs_computation(const EstimationGraph& egraph) const;

		BidirectedGraph bidirected_graph_computation(const std::vector< std::pair<Neighborhood, Neighborhood> >& decomposition) const;

		std::vector< Neighborhood > get_ps(const UndirectedGraph& ugraph) const;
		void bron_kerbosch_all_cliques(const UndirectedGraph& ugraph, const Neighborhood& r, Neighborhood p, Neighborhood x, Cliques& cliques, const int& superior_to = 0) const;
		std::vector< std::pair<Neighborhood, Neighborhood> > get_component_decomposition(const Neighborhood& start, const Neighborhood& vertices, const UndirectedGraph& ugraph) const;

		DistributionCollection distribution_collection_computation(const UndirectedGraph& ugraph, bool parametric) const;
		void component_decompositions_computation(unsigned int component_nb_variables, const Neighborhood& setted_vertices, const Cliques& cliques, const std::vector< std::pair<Neighborhood, Neighborhood> >& decomposition, const UndirectedGraph& ugraph, std::vector< std::vector< std::pair<Neighborhood, Neighborhood> > >& decompositions) const;
		ReestimationCollection<T> decomposition_reestimation_collection_computation(const std::vector< std::pair<Neighborhood, Neighborhood> >& decomposition) const;
		std::pair<double, DistributionCollection > distribution_estimation(const ReestimationCollection<T>& reestimation_collection, bool parametric) const;
		//std::pair<double, DistributionCollection > type_parametric_estimation(const ReestimationCollection<T>& reestimation_collection) const;
		DistributionCollection distribution_collections_agglomerate(const std::vector< DistributionCollection >& distribution_collections) const;

		Neighborhood get_neighborhood(const Neighborhood& vertices, const UndirectedGraph& ugraph) const;
		Neighborhood get_neighborhood(unsigned int vertex, const UndirectedGraph& ugraph) const;
		Neighborhood get_intersection(const Neighborhood& vertices0, const Neighborhood& vertices1) const;
		Neighborhood get_complement(const Neighborhood& intersection, const Neighborhood& vertices) const;
		Neighborhood get_setminus(const Neighborhood& set, const Neighborhood& minus) const;
		Neighborhood get_setminus(const Neighborhood& set, unsigned int minus) const;
		Neighborhood get_setplus(const Neighborhood& set, unsigned int plus) const;

		Reestimation<T>* get_reestimation(unsigned int vertex) const;
		Stat_trees::DiscreteMultivariateReestimation<T>* get_reestimation(const Neighborhood& vertices) const;
		DiscreteUnivariateConditionalReestimation<T>* get_reestimation(const Neighborhood& covertices, unsigned int vertex) const;
		DiscreteMultivariateConditionalReestimation<T>* get_reestimation(const Neighborhood& covertices, const Neighborhood& vertices) const;

		//DistributionCollection distribution_estimation(const std::vector< ReestimationCollection<T> >& reestimation_collections) const;
		bool same_neighborhood(const Neighborhood& n0, const Neighborhood& n1) const;
};

class DiscreteGraphicalData : public DiscreteGraphicalReestimation<int>
{
	public:
		DiscreteGraphicalData();
		~DiscreteGraphicalData();
		DiscreteGraphicalData(unsigned int inb_variables, const std::vector< std::vector<int> >& vectors);
		DiscreteGraphicalData(const DiscreteGraphicalData& dgd);
		DiscreteGraphicalData(unsigned int inb_variables);
};

class DiscreteGraphicalDistribution
{
	public:
		DiscreteGraphicalDistribution();
		~DiscreteGraphicalDistribution();
		DiscreteGraphicalDistribution(const DiscreteGraphicalDistribution& dgd);
		DiscreteGraphicalDistribution(unsigned int inb_variables, const UndirectedGraph& iugraph, const DistributionCollection& idistribution_collection);
		double get_mass(const std::vector<int>& event, bool log_computation);
		std::pair<double, double> compare(const UndirectedGraph& graph) const;
		unsigned int get_nb_variables() const;
		unsigned int get_nb_parameters() const;
		unsigned int get_nb_edges() const;
		UndirectedGraph get_ugraph() const;
		std::string get_formula() const;
		std::vector<int> simulation() const;
		DiscreteGraphicalData* simulation(unsigned int number) const;

	protected:
		unsigned int nb_variables;
		unsigned int nb_parameters;
		std::map<std::vector<int>, double> mass;
		DistributionCollection distribution_collection;
		UndirectedGraph ugraph;
		double mass_computation(const std::vector<int>& event) const;
};


class DiscreteGraphicalParametric : public DiscreteGraphicalDistribution
{
	public:
		DiscreteGraphicalParametric();
		~DiscreteGraphicalParametric();
		DiscreteGraphicalParametric(const DiscreteGraphicalParametric& dgp);
		DiscreteGraphicalParametric(unsigned int inb_variables, const BidirectedGraph& ibgraph, const DistributionCollection& idistribution_collection);
		BidirectedGraph get_bgraph() const;

	protected:
		BidirectedGraph bgraph;

	private:
		UndirectedGraph moral_graph_computation(const BidirectedGraph& ibgraph) const;	
};

#include "graphical_distribution.hpp"
#endif
