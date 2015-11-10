#ifndef DISCRETE_GRAPHICAL_DISTRIBUTION_CPP
#define DISCRETE_GRAPHICAL_DISTRIBUTION_CPP

#include "graphical_distribution.h"
#include "graphical_distribution.hpp"

DiscreteGraphicalDistribution::DiscreteGraphicalDistribution()
{
}

DiscreteGraphicalDistribution::~DiscreteGraphicalDistribution()
{
}

DiscreteGraphicalDistribution::DiscreteGraphicalDistribution(const DiscreteGraphicalDistribution& dgd)
{
	means = dgd.means;
	cvm = dgd.cvm;
	histogram = dgd.histogram;
	ugraph = dgd.ugraph;
	nparameters = dgd.nparameters;
}

DiscreteGraphicalDistribution::DiscreteGraphicalDistribution(const UndirectedGraph& iugraph, const DiscreteMultivariateHistogram& hist)
{
	ugraph = iugraph;
	nparameters = 0;
	if(TriangulatedGraphChecking())
		TriangulateHistogramDecomposition(hist);
	else
		NonTriangulateHistogramDecomposition(hist);
	means = hist.GetMeans();
	CVMComputation(hist);
}

/*void DiscreteGraphicalDistribution::GraphTriangulation()
{
	Neighborhood p, separator;
	std::set< Neighborhood > separators;
	for(unsigned int i = 0; i < boost::num_vertices(ugraph); ++i)
		p.insert(i);
	BronKerbosch(Neighborhood(), p, Neighborhood(), cliques);
	Cliques::iterator itb;
	for(Cliques::iterator it = cliques.begin(); it != cliques.end(); ++it)
	{
		itb = it;
		++itb;
		while(itb != cliques.end())
		{
		}
	}
}*/

unsigned int DiscreteGraphicalDistribution::GetNbParameters() const
{
	return nparameters;
}

bool DiscreteGraphicalDistribution::TriangulatedGraphChecking() const
{
	UndirectedGraphTraits::vertex_iterator vi, vi_end;
	UndirectedGraphTraits::adjacency_iterator ai, ai_end;
	Neighborhood labelled_vertices, non_labelled_vertices, i_labelled_vertices, adj_labelled_vertices;
	Neighborhood::iterator itl0, itl1;
	labelled_vertices.insert(0);
	unsigned int max, i_max;
	bool success = true, end;
	for(unsigned int i = 1; i < boost::num_vertices(ugraph); ++i)
	{
		non_labelled_vertices.insert(i);
	}
	for(unsigned int i = 1; i < boost::num_vertices(ugraph); ++i)
	{
		max = 0;
		for(Neighborhood::iterator it = non_labelled_vertices.begin(); it != non_labelled_vertices.end(); ++it)
		{
			adj_labelled_vertices.clear();
			for(boost::tie(ai, ai_end) = boost::adjacent_vertices(*it, ugraph); ai != ai_end; ++ai)
			{
				if(labelled_vertices.find(*ai) != labelled_vertices.end())
					adj_labelled_vertices.insert(*ai);
			}
			if(adj_labelled_vertices.size() > max)
			{
				max = adj_labelled_vertices.size();
				i_max = *it;
				i_labelled_vertices = adj_labelled_vertices;
			}
		}
		end = (i_labelled_vertices.size() < 2);
		labelled_vertices.insert(i_max);
		non_labelled_vertices.erase(i_max);
		if(!end)
		{
			itl0 = i_labelled_vertices.begin();
			itl1 = itl0;
			++itl1;
			while(!end && success)
			{
				while(itl1 != i_labelled_vertices.end() && success)
				{
					if(!boost::edge(*itl0, *itl1, ugraph).second)
						success = false;
					++itl1;
				}
				++itl0;
				itl1 = itl0;
				++itl1;
				if(itl1 == i_labelled_vertices.end())
					end = true;
			}
			if(!success)
				return false;
		}
	}
	return success;
}

void DiscreteGraphicalDistribution::TriangulateHistogramDecomposition(const DiscreteMultivariateHistogram& hist)
{
	Cliques cliques, ccliques;
	Neighborhood separator, p;
	for(unsigned int i = 0; i < boost::num_vertices(ugraph); ++i)
		p.insert(i);
	BronKerbosch(Neighborhood(), p, Neighborhood(), cliques);
	CliqueGraph cgraph;
	clique_vertex cv;
	separator_edge se;
	DiscreteMultivariateHistogram *dmh = NULL;
	for(Cliques::iterator it = cliques.begin(); it != cliques.end(); ++it)
	{
		cv.clique = (*it).second;
		boost::add_vertex(cv, cgraph);
		dmh = hist.GetMarginalHistogram(cv.clique);
		nparameters += dmh->GetNbBreaks()-1;
		(histogram.first).push_back(new DiscreteMultivariateHistogram(*dmh));
		delete dmh;
	}
	std::multimap<unsigned int, boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge> > card_sep_score;
	CliqueGraphTraits::vertex_iterator v0, v1, v0_begin, v1_begin, v0_end, v1_end;
	boost::tie(v0_begin, v1_end) = boost::vertices(cgraph);
	v0_end = v1_end;
	--v0_end;
	for(v0 = v0_begin; v0 != v0_end; ++v0)
	{
		v1_begin = v0;
		++v1_begin;
		for(v1 = v1_begin; v1 != v1_end; ++v1)
		{
			(se.separator).clear();
			for(Neighborhood::iterator it = (cgraph[*v0].clique).begin(); it != (cgraph[*v0].clique).end(); ++it)
			{
				if((cgraph[*v1].clique).find(*it) != (cgraph[*v1].clique).end())
					(se.separator).insert(*it);
			}
			card_sep_score.insert(std::pair<unsigned int, boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge> >((se.separator).size(), boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge>(v0, v1, se)));
		}
	}
	std::multimap<unsigned int, boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge> >::iterator itcss;
	CliqueGraphTraits::vertices_size_type d[boost::num_vertices(cgraph)];
	while(boost::num_edges(cgraph) != boost::num_vertices(cgraph)-1)
	{
		std::fill_n(d, boost::num_vertices(cgraph), 0);
		itcss = card_sep_score.end();
		--itcss;
		boost::breadth_first_search(cgraph, *(((*itcss).second).get<0>()), boost::visitor(boost::make_bfs_visitor(boost::record_distances(d, boost::on_tree_edge()))));
		while(d[*(((*itcss).second).get<1>())] != 0)
		{
			std::fill_n(d, boost::num_vertices(cgraph), 0);
			card_sep_score.erase(itcss);
			itcss = card_sep_score.end();
			--itcss;
			boost::breadth_first_search(cgraph, *(((*itcss).second).get<0>()), boost::visitor(boost::make_bfs_visitor(boost::record_distances(d, boost::on_tree_edge()))));
		}
		boost::add_edge(*(((*itcss).second).get<0>()), *(((*itcss).second).get<1>()), ((*itcss).second).get<2>(), cgraph);
		if(((((*itcss).second).get<2>()).separator).size() > 0)
		{
			dmh = hist.GetMarginalHistogram((((*itcss).second).get<2>()).separator);
			nparameters += dmh->GetNbBreaks()-1;
			(histogram.second).push_back(new DiscreteMultivariateHistogram(*dmh));
			delete dmh;
		}
		card_sep_score.erase(itcss);
	}
}

void DiscreteGraphicalDistribution::NonTriangulateHistogramDecomposition(const DiscreteMultivariateHistogram& hist)
{
	Cliques cliques, ccliques;
	Neighborhood separator, p;
	for(unsigned int i = 0; i < boost::num_vertices(ugraph); ++i)
		p.insert(i);
	BronKerbosch(Neighborhood(), p, Neighborhood(), cliques);
	CliqueGraph cgraph;
	clique_vertex cv;
	separator_edge se;
	DiscreteMultivariateHistogram *dmh = NULL;
	for(Cliques::iterator it = cliques.begin(); it != cliques.end(); ++it)
	{
		cv.clique = (*it).second;
		boost::add_vertex(cv, cgraph);
		dmh = hist.GetMarginalHistogram(cv.clique);
		nparameters += dmh->GetNbBreaks()-1;
		(histogram.first).push_back(new DiscreteMultivariateHistogram(*dmh));
		delete dmh;
	}
	std::multimap<unsigned int, boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge> > card_sep_score;
	CliqueGraphTraits::vertex_iterator v0, v1, v0_begin, v1_begin, v0_end, v1_end;
	boost::tie(v0_begin, v1_end) = boost::vertices(cgraph);
	v0_end = v1_end;
	--v0_end;
	for(v0 = v0_begin; v0 != v0_end; ++v0)
	{
		v1_begin = v0;
		++v1_begin;
		for(v1 = v1_begin; v1 != v1_end; ++v1)
		{
			(se.separator).clear();
			for(Neighborhood::iterator it = (cgraph[*v0].clique).begin(); it != (cgraph[*v0].clique).end(); ++it)
			{
				if((cgraph[*v1].clique).find(*it) != (cgraph[*v1].clique).end())
					(se.separator).insert(*it);
			}
			card_sep_score.insert(std::pair<unsigned int, boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge> >((se.separator).size(), boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge>(v0, v1, se)));
		}
	}
	std::multimap<unsigned int, boost::tuple<CliqueGraphTraits::vertex_iterator, CliqueGraphTraits::vertex_iterator, separator_edge> >::iterator itcss;
	CliqueGraphTraits::vertices_size_type d[boost::num_vertices(cgraph)];
	while(boost::num_edges(cgraph) != boost::num_vertices(cgraph)-1)
	{
		std::fill_n(d, boost::num_vertices(cgraph), 0);
		itcss = card_sep_score.end();
		--itcss;
		boost::breadth_first_search(cgraph, *(((*itcss).second).get<0>()), boost::visitor(boost::make_bfs_visitor(boost::record_distances(d, boost::on_tree_edge()))));
		while(d[*(((*itcss).second).get<1>())] != 0)
		{
			std::fill_n(d, boost::num_vertices(cgraph), 0);
			card_sep_score.erase(itcss);
			itcss = card_sep_score.end();
			--itcss;
			boost::breadth_first_search(cgraph, *(((*itcss).second).get<0>()), boost::visitor(boost::make_bfs_visitor(boost::record_distances(d, boost::on_tree_edge()))));
		}
		boost::add_edge(*(((*itcss).second).get<0>()), *(((*itcss).second).get<1>()), ((*itcss).second).get<2>(), cgraph);
		if(((((*itcss).second).get<2>()).separator).size() > 0)
		{
			dmh = hist.GetMarginalHistogram((((*itcss).second).get<2>()).separator);
			nparameters += dmh->GetNbBreaks()-1;
			(histogram.second).push_back(new DiscreteMultivariateHistogram(*dmh));
			delete dmh;
		}
		card_sep_score.erase(itcss);
	}
}

std::pair<double, double> DiscreteGraphicalDistribution::Compare(const UndirectedGraph& graph) const
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

Moments DiscreteGraphicalDistribution::GetMeans() const
{
	return means;
}

PairedMoments DiscreteGraphicalDistribution::GetCVM() const
{
	return cvm;
}

double DiscreteGraphicalDistribution::ComputeLogLikelihood(const DiscreteMultivariateHistogram& hist) const
{
	double loglikelihood = 0;
	std::vector< std::vector<int> > vectors = hist.GetEvents();
	for(std::vector< std::vector<int> >::const_iterator it = vectors.begin(); it != vectors.end(); ++it)
		loglikelihood += log(GetMass(*it))*hist.GetNbEvent(*it);
	return loglikelihood;
}

UndirectedGraph DiscreteGraphicalDistribution::GetGraph() const
{
	return ugraph;
}

std::ostream& DiscreteGraphicalDistribution::Summary(std::ostream& os) const
{
	return os;
}
/*
DiscreteGraphicalEstimation* DiscreteGraphicalDistribution::Simulate(unsigned int number) const
{
	DiscreteMultivariateHistogram* hist = new DiscreteMultivariateHistogram(boost::num_vertices(ugraph));
	double gen;
	unsigned int which;
	DiscreteMultivariateMass::const_iterator it;
	double r;
	while(hist->GetNbEvents() != number)
	{
		it = mass.begin();
		r = generator();
		advance(it, (unsigned int)(r * mass.size()));
		r = generator();
		if(r < (*it).second)
			hist->AddEvent((*it).first);
	}
	DiscreteGraphicalEstimation* dge = new DiscreteGraphicalEstimation(*hist);
	delete hist;
	return dge;
}*/

double DiscreteGraphicalDistribution::GetMass(const std::vector<int>& event) const
{
	double logp = 0;
	for(std::vector< DiscreteMultivariateHistogram* >::const_iterator it = (histogram.first).begin(); it != (histogram.first).end(); ++it)
		logp += log((*it)->GetNbEvent(event)) - log((*it)->GetNbEvents());
	for(std::vector< DiscreteMultivariateHistogram* >::const_iterator it = (histogram.second).begin(); it != (histogram.second).end(); ++it)
		logp -= (log((*it)->GetNbEvent(event)) - log((*it)->GetNbEvents()));
	return exp(logp);
}

void DiscreteGraphicalDistribution::CVMComputation(const DiscreteMultivariateHistogram& hist)
{
	UndirectedGraphTraits::edge_iterator ei, ei_end;
	cvm.resize(boost::num_vertices(ugraph), Moments(boost::num_vertices(ugraph), 0));
	PairedMoments mcvm = hist.GetCVM();
	for(boost::tie(ei, ei_end) = boost::edges(ugraph); ei != ei_end; ++ei)
	{
		cvm[boost::source(*ei, ugraph)][boost::source(*ei, ugraph)] = mcvm[boost::source(*ei, ugraph)][boost::source(*ei, ugraph)];
		cvm[boost::source(*ei, ugraph)][boost::target(*ei, ugraph)] = mcvm[boost::source(*ei, ugraph)][boost::target(*ei, ugraph)];
		cvm[boost::target(*ei, ugraph)][boost::source(*ei, ugraph)] = mcvm[boost::target(*ei, ugraph)][boost::source(*ei, ugraph)];
		cvm[boost::target(*ei, ugraph)][boost::target(*ei, ugraph)] = mcvm[boost::target(*ei, ugraph)][boost::target(*ei, ugraph)];
	}
}

unsigned int DiscreteGraphicalDistribution::GetNbVertices() const
{
	return boost::num_vertices(ugraph);
}

BidirectedGraph DiscreteGraphicalDistribution::GraphOrientation() const
{
	Cliques cliques;
	Neighborhood p,x;
	for(unsigned int i = 0; i < boost::num_vertices(ugraph); ++i)
		p.insert(i);
	BronKerbosch(Neighborhood(), p, x, cliques);
	BidirectedGraph bgraph(boost::num_vertices(ugraph));
	for(Cliques::iterator itc = cliques.begin(); itc != cliques.end(); ++itc)
	{
		std::cout << "Clique de taille " << (*itc).first << " :";
		for(Neighborhood::iterator itn0 = ((*itc).second).begin(); itn0 != ((*itc).second).end(); ++itn0)
		{
			std::cout << *itn0 << " ";
			for(Neighborhood::iterator itn1 = ((*itc).second).begin(); itn1 != ((*itc).second).end(); ++itn1)
			{
				if(itn1 != itn0)
					boost::add_edge(*itn1, *itn0, bgraph);
			}
		}
		std::cout << std::endl;
	}
	return bgraph;
}

void DiscreteGraphicalDistribution::BronKerbosch(const Neighborhood& r, Neighborhood p, Neighborhood x, Cliques& cliques) const
{
	if(p.size() == 0 && x.size() == 0)
	{
		cliques.insert(std::pair<unsigned int, Neighborhood>(r.size(), r));
	} else {
		Neighborhood rold, rnew, pnew, xnew, vertices = p;
		UndirectedGraphTraits::adjacency_iterator ai, ai_end;
		for(Neighborhood::iterator it = vertices.begin(); it != vertices.end(); ++it)
		{
			rnew = r;
			rnew.insert(*it);
			xnew.clear();
			pnew.clear();
			p.erase(*it);
			for(boost::tie(ai, ai_end) = boost::adjacent_vertices(*it, ugraph); ai != ai_end; ++ai)
			{
				if(p.find(*ai) != p.end())
					pnew.insert(*ai);
				if(x.find(*ai) != x.end())
					xnew.insert(*ai);
			}
			BronKerbosch(rnew, pnew, xnew, cliques);
			x.insert(*it);
		}
	}
}

#endif
