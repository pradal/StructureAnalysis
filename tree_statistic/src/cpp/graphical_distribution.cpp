#ifndef DISCRETE_GRAPHICAL_DISTRIBUTION_CPP
#define DISCRETE_GRAPHICAL_DISTRIBUTION_CPP

#include "graphical_distribution.h"

DiscreteGraphicalData::DiscreteGraphicalData() : DiscreteGraphicalReestimation<int>()
{
}

DiscreteGraphicalData::~DiscreteGraphicalData()
{
}

DiscreteGraphicalData::DiscreteGraphicalData(unsigned int inb_variables, const std::vector< std::vector<int> >& vectors)
{
	nb_variables = inb_variables;
	std::map< std::vector<int>, int>::iterator ith;
	nb_elements = 0;
	for(std::vector< std::vector<int> >::const_iterator it = vectors.begin(); it != vectors.end(); ++it)
	{
		nb_elements++;
		ith = histogram.find(*it);
		if(ith == histogram.end())
			histogram.insert(std::pair<std::vector<int>, int>(*it, 1));
		else
			(*ith).second += 1;
	}
}

DiscreteGraphicalData::DiscreteGraphicalData(const DiscreteGraphicalData& dgd) : DiscreteGraphicalReestimation<int>(dgd)
{
}

DiscreteGraphicalData::DiscreteGraphicalData(unsigned int inb_variables) : DiscreteGraphicalReestimation<int>(inb_variables)
{
}

DiscreteGraphicalDistribution::DiscreteGraphicalDistribution()
{
	nb_variables = 0;
	nb_parameters = 0;
}

DiscreteGraphicalDistribution::~DiscreteGraphicalDistribution()
{
}

DiscreteGraphicalDistribution::DiscreteGraphicalDistribution(const DiscreteGraphicalDistribution& dgd)
{
	nb_variables = dgd.nb_variables;
	nb_parameters = dgd.nb_parameters;
	mass = dgd.mass;
	ugraph = dgd.ugraph;
	distribution_collection.decomposition = dgd.distribution_collection.decomposition;
	distribution_collection.DUD = dgd.distribution_collection.DUD;
	distribution_collection.DUCD = dgd.distribution_collection.DUCD;
	distribution_collection.DMD = dgd.distribution_collection.DMD;
	distribution_collection.DMCD = dgd.distribution_collection.DMCD;
}

DiscreteGraphicalDistribution::DiscreteGraphicalDistribution(unsigned int inb_variables, const UndirectedGraph& iugraph, const DistributionCollection& idistribution_collection)
{
	nb_variables = inb_variables;
	ugraph = iugraph;
	distribution_collection.decomposition = idistribution_collection.decomposition;
	nb_parameters = 0;
	distribution_collection.DUD = idistribution_collection.DUD;
	for(std::vector< Distribution* >::iterator it = (distribution_collection.DUD).begin(); it != (distribution_collection.DUD).end(); ++it)
		nb_parameters += (*it)->nb_parameter;
	distribution_collection.DMD = idistribution_collection.DMD;
	for(std::vector< Stat_trees::DiscreteMultivariateDistribution* >::iterator it = (distribution_collection.DMD).begin(); it != (distribution_collection.DMD).end(); ++it)
		nb_parameters += (*it)->nb_parameter;
	distribution_collection.DUCD = idistribution_collection.DUCD;
	for(std::vector< DiscreteUnivariateConditionalDistribution* >::iterator it = (distribution_collection.DUCD).begin(); it != (distribution_collection.DUCD).end(); ++it)
		nb_parameters += (*it)->get_nb_parameters();
	distribution_collection.DMCD = idistribution_collection.DMCD;
	for(std::vector< DiscreteMultivariateConditionalDistribution* >::iterator it = (distribution_collection.DMCD).begin(); it != (distribution_collection.DMCD).end(); ++it)
		nb_parameters += (*it)->get_nb_parameters();
}

unsigned int DiscreteGraphicalDistribution::get_nb_edges() const
{
	return boost::num_edges(ugraph);
}

std::string DiscreteGraphicalDistribution::get_formula() const
{
	std::string decomposition = "";
	Neighborhood::const_iterator ni_end;
	for(std::vector< std::pair<Neighborhood, Neighborhood> >::const_iterator it = (distribution_collection.decomposition).begin(); it != (distribution_collection.decomposition).end(); ++it)
	{
		decomposition += " P[";
		ni_end = ((*it).second).end();
		ni_end--;
		for(Neighborhood::const_iterator ni = ((*it).second).begin(); ni != ni_end; ++ni)
			decomposition += "v" + toString(*ni) +", ";
		if((*it).first.size() == 0)
		{
			decomposition += "v" + toString(*ni_end) + "]";
		} else {
			decomposition += "v" + toString(*ni_end) + " | ";
			ni_end = ((*it).first).end();
			ni_end--;
			for(Neighborhood::const_iterator ni = ((*it).first).begin(); ni != ni_end; ++ni)
				decomposition += "v" + toString(*ni)+ ", ";
			decomposition += "v" + toString(*ni_end) + "]";
		}
	}
	return decomposition;
}

std::vector<int> DiscreteGraphicalDistribution::simulation() const
{
	unsigned int dud = 0, ducd = 0, dmd = 0, dmcd = 0;
	std::vector<int> event(nb_variables);
	for(std::vector< std::pair<Neighborhood, Neighborhood > >::const_iterator it = (distribution_collection.decomposition).begin(); it != (distribution_collection.decomposition).end(); ++it)
	{
		if(((*it).first).size() == 0)
		{
			if(((*it).second).size() == 1)
			{
				event[(*((*it).second).begin())] = distribution_collection.DUD[dud]->simulation();
				dud++;
			} else {
				std::vector<unsigned int> mevent = distribution_collection.DMD[dmd]->simulation();
				Neighborhood::const_iterator itn = ((*it).second).begin();
				for(std::vector<unsigned int>::iterator ite = mevent.begin(); ite != mevent.end(); ++ite)
				{
					event[*itn] = *ite;
					++itn;
				}
				dmd++;
			}
		} else {
			std::vector<int> cevent(((*it).first).size());
			Neighborhood::const_iterator itn = ((*it).first).begin();
			for(std::vector<int>::iterator ite = cevent.begin(); ite != cevent.end(); ++ite)
			{
				*ite = event[*itn];
				++itn;
			}
			if(((*it).second).size() == 1)
			{
				event[(*((*it).second).begin())] = distribution_collection.DUCD[ducd]->simulation(cevent);
				ducd++;
			} else {
				std::vector<int> mevent = distribution_collection.DMCD[dmcd]->simulation(cevent);
				itn = ((*it).second).begin();
				for(std::vector<int>::iterator ite = mevent.begin(); ite != mevent.end(); ++ite)
				{
					event[*itn] = *ite;
					++itn;
				}
				dmcd++;
			}
		}
	}
	return event;
}

DiscreteGraphicalData* DiscreteGraphicalDistribution::simulation(unsigned int number) const
{
	DiscreteGraphicalData *dgd = new DiscreteGraphicalData(nb_variables);
	for(unsigned int i = 0; i < number; ++i)
		dgd->update(simulation(), 1);
	dgd->nb_elements_computation();
	return dgd;
}

double DiscreteGraphicalDistribution::get_mass(const std::vector<int>& event, bool log_computation)
{
	std::map<std::vector<int>, double>::iterator it = mass.find(event);
	double p;
	if(it == mass.end())
	{
		p = mass_computation(event);
		mass.insert(std::pair<std::vector<int>, double>(event, p));
	} else {
		p = (*it).second;
	}
	if(log_computation)
		return log(p);
	else
		return p;
}

double DiscreteGraphicalDistribution::mass_computation(const std::vector<int>& event) const
{
	unsigned int dud = 0, ducd = 0, dmd = 0, dmcd = 0;
	double log_p = 0;
	for(std::vector< std::pair<Neighborhood, Neighborhood > >::const_iterator it = (distribution_collection.decomposition).begin(); it != (distribution_collection.decomposition).end(); ++it)
	{
		if(((*it).first).size() == 0)
		{
			if(((*it).second).size() == 1)
			{
				log_p += log(distribution_collection.DUD[dud]->mass[event[(*((*it).second).begin())]]);
				dud++;
			} else {
				std::vector<unsigned int> mevent(((*it).second).size());
				Neighborhood::const_iterator itn = ((*it).second).begin();
				for(std::vector<unsigned int>::iterator ite = mevent.begin(); ite != mevent.end(); ++ite)
				{
					*ite = event[*itn];
					++itn;
				}
				log_p += distribution_collection.DMD[dmd]->get_mass(mevent,true);
				dmd++;
			}
		} else {
			std::vector<int> cevent(((*it).first).size());
			Neighborhood::const_iterator itn = ((*it).first).begin();
			for(std::vector<int>::iterator ite = cevent.begin(); ite != cevent.end(); ++ite)
			{
				*ite = event[*itn];
				++itn;
			}
			if(((*it).second).size() == 1)
			{
				log_p += distribution_collection.DUCD[ducd]->get_mass(cevent, event[(*((*it).second).begin())], true);
				ducd++;
			} else {
				std::vector<int> mevent(((*it).second).size());
				itn = ((*it).second).begin();
				for(std::vector<int>::iterator ite = mevent.begin(); ite != mevent.end(); ++ite)
				{
					*ite = event[*itn];
					++itn;
				}
				log_p += distribution_collection.DMCD[dmcd]->get_mass(cevent, mevent,true);
				dmcd++;
			}
		}
	}
	return exp(log_p);
}

unsigned int DiscreteGraphicalDistribution::get_nb_variables() const
{
	return nb_variables;
}

unsigned int DiscreteGraphicalDistribution::get_nb_parameters() const
{
	return nb_parameters;
}

UndirectedGraph DiscreteGraphicalDistribution::get_ugraph() const
{
	return ugraph;
}

std::pair<double, double> DiscreteGraphicalDistribution::compare(const UndirectedGraph& graph) const 
{
	unsigned int tp = 0, fp = 0, fn = 0;
	UndirectedGraphTraits::edge_iterator ei, ei_end;
	for(boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
	{
		if(boost::edge(boost::source(*ei, graph), boost::target(*ei, graph), ugraph).second)
			tp++;
		else
			fn++;
	}
	for(boost::tie(ei, ei_end) = boost::edges(ugraph); ei != ei_end; ++ei)
	{
		if(!(boost::edge(boost::source(*ei, ugraph), boost::target(*ei, ugraph), graph).second))
			fp++;
	}
	if(tp+fp == 0 && tp+fn == 0)
		return std::pair<double, double>(0, 1);
	if(tp+fp == 0 && tp+fn != 0)
		return std::pair<double, double>(((double)tp)/((double)(tp+fn)),1);
	if(tp+fp != 0 && tp+fn == 0)
		return std::pair<double, double>(0, ((double)tp)/((double)(tp+fp)));
	if(tp+fp != 0 && tp+fn != 0)
		return std::pair<double, double>(((double)tp)/((double)(tp+fn)),((double)tp)/((double)(tp+fp)));
}

DiscreteGraphicalParametric::DiscreteGraphicalParametric() : DiscreteGraphicalDistribution()
{
}

DiscreteGraphicalParametric::~DiscreteGraphicalParametric()
{
}

DiscreteGraphicalParametric::DiscreteGraphicalParametric(const DiscreteGraphicalParametric& dgp) : DiscreteGraphicalDistribution(dgp)
{
	bgraph = dgp.bgraph;
}

DiscreteGraphicalParametric::DiscreteGraphicalParametric(unsigned int inb_variables, const BidirectedGraph& ibgraph, const DistributionCollection& idistribution_collection) : DiscreteGraphicalDistribution(inb_variables, moral_graph_computation(ibgraph), idistribution_collection)
{
	bgraph = ibgraph;
}

BidirectedGraph DiscreteGraphicalParametric::get_bgraph() const
{
	return bgraph;
}

UndirectedGraph DiscreteGraphicalParametric::moral_graph_computation(const BidirectedGraph& ibgraph) const
{
	UndirectedGraph iugraph(boost::num_vertices(ibgraph));
	BidirectedGraphTraits::vertex_iterator iv, iv_end;
	BidirectedGraphTraits::in_edge_iterator iei0, iei0_begin, iei0_end, iei1, iei1_begin, iei1_end;
	for(boost::tie(iv, iv_end) = boost::vertices(ibgraph); iv != iv_end; ++iv)
	{
		if(boost::in_degree(*iv, ibgraph) > 1)
		{
			boost::tie(iei0_begin, iei1_end) = boost::in_edges(*iv, ibgraph);
			iei0_end = iei1_end;
			--iei0_end;
			for(iei0 = iei0_begin; iei0 != iei0_end; ++iei0)
			{
				iei1_begin = iei0;
				++iei1_begin;
				for(iei1 = iei1_begin; iei1 != iei1_end; ++iei1)
				{
					if(!boost::edge(boost::source(*iei0, ibgraph), boost::source(*iei1, ibgraph), iugraph).second)
						boost::add_edge(boost::source(*iei0, ibgraph), boost::source(*iei1, ibgraph), iugraph);
				}
			}
		}
	}
	BidirectedGraphTraits::edge_iterator ei, ei_end;
	for(boost::tie(ei, ei_end) = boost::edges(ibgraph); ei != ei_end; ++ei)
	{
		if(!boost::edge(boost::source(*ei, ibgraph), boost::target(*ei, ibgraph), iugraph).second)
			boost::add_edge(boost::source(*ei, ibgraph), boost::target(*ei, ibgraph), iugraph);
	}
	return iugraph;
}
#endif
