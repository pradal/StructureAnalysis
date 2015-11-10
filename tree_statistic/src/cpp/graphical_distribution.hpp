#ifndef DISCRETE_GRAPHICAL_DISTRIBUTION_HPP
#define DISCRETE_GRAPHICAL_DISTRIBUTION_HPP

template<typename T> DiscreteGraphicalReestimation<T>::DiscreteGraphicalReestimation()
{
	nb_variables = 0;
	nb_elements = 0;
}

template<typename T> DiscreteGraphicalReestimation<T>::~DiscreteGraphicalReestimation()
{
}

template<typename T> DiscreteGraphicalReestimation<T>::DiscreteGraphicalReestimation(const DiscreteGraphicalReestimation<T>& dgr)
{
	histogram = dgr.histogram;
	nb_variables = dgr.nb_variables;
	nb_elements = dgr.nb_elements;
}

template<typename T> std::ostream& DiscreteGraphicalReestimation<T>::display(std::ostream& os) const
{
	std::vector< std::vector<std::string> > Output(histogram.size()+1, std::vector<std::string>(nb_variables+1,""));
	for(unsigned int i = 0; i < nb_variables; ++i)
		Output[0][i] = "Vertex "+ toString(i);
	Output[0][nb_variables] = "Occurancy";
	for(typename std::map<std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
	{
		for(unsigned int i = 0; i < nb_variables; ++i)
			Output[distance(histogram.begin(), it)+1][i] = toString(((*it).first)[i]);
		Output[distance(histogram.begin(), it)+1][nb_variables] = toString((*it).second);
	}
	for(std::vector< std::vector<std::string> >::iterator it = Output.begin(); it != Output.end(); ++it)
	{
		for(std::vector<std::string>::iterator itb = (*it).begin(); itb != (*it).end(); ++itb)
			os << " " << *itb << " ";
		os << std::endl;
	}
	return os;
}

template<typename T> DiscreteGraphicalReestimation<T>::DiscreteGraphicalReestimation(unsigned int inb_variables)
{
	nb_variables = inb_variables;
}

template<typename T> void DiscreteGraphicalReestimation<T>::update(const std::vector<int>& event, T frequency)
{
	assert(event.size() == nb_variables);
	typename std::map<std::vector<int>, T>::iterator it = histogram.find(event);
	if(it == histogram.end())
		histogram.insert(std::pair<std::vector<int>, T>(event, frequency));
	else
		(*it).second += frequency;
}

template<typename T> void DiscreteGraphicalReestimation<T>::nb_elements_computation()
{
	nb_elements = 0;
	for(typename std::map<std::vector<int>, T>::iterator it = histogram.begin(); it != histogram.end(); ++it)
		nb_elements += (*it).second;
}

template<typename T> EstimationGraph DiscreteGraphicalReestimation<T>::graph_estimation(int algo, int max_step) const
{
	switch(algo)
	{
		case CLT:
			return algorithm_CLT(max_step);
			break;
		case ML :
			return algorithm_ML(max_step);
			break;
	}
}

template<typename T> EstimationGraph DiscreteGraphicalReestimation<T>::graph_estimation(int algo, double min_threshold) const
{
	switch(algo)
	{
		case MRNET :
			return algorithm_MRNET(min_threshold);
			break;
		case RELEVANCE :
		case CLR :
			return algorithms_threshold(algo, min_threshold);
			break;
	}
}

template<typename T> EstimationGraph DiscreteGraphicalReestimation<T>::graph_estimation(int algo) const
{
	switch(algo)
	{
		case CLT :
		case ML :
			return graph_estimation(algo, -1);
			break;
		case MRNET :
		case RELEVANCE :
		case CLR :
			return graph_estimation(algo, -1*D_INF);
	}
}


template<typename T> DiscreteGraphicalDistribution* DiscreteGraphicalReestimation<T>::distribution_estimation(int algo, int max_step) const
{
	UndirectedGraph ugraph = undirected_graph_computation(graph_estimation(algo, max_step));
	DiscreteGraphicalDistribution* dgd = new DiscreteGraphicalDistribution(nb_variables, ugraph, distribution_collection_computation(ugraph, false));
	return dgd;
}

template<typename T> DiscreteGraphicalDistribution* DiscreteGraphicalReestimation<T>::distribution_estimation(int algo, double min_threshold) const
{
	UndirectedGraph ugraph = undirected_graph_computation(graph_estimation(algo, min_threshold));
	DiscreteGraphicalDistribution* dgd = new DiscreteGraphicalDistribution(nb_variables, ugraph, distribution_collection_computation(ugraph, false));
	return dgd;
}

template<typename T> std::vector< DiscreteGraphicalDistribution* > DiscreteGraphicalReestimation<T>::stepwise_distribution_estimation(int algo) const
{
	std::vector< UndirectedGraph > ugraphs = undirected_graphs_computation(graph_estimation(algo));
	std::vector< DiscreteGraphicalDistribution* > dgds(ugraphs.size());
	for(std::vector< UndirectedGraph >::iterator it = ugraphs.begin(); it != ugraphs.end(); ++it)
		dgds[distance(ugraphs.begin(), it)] = new DiscreteGraphicalDistribution(nb_variables, *it, distribution_collection_computation(*it, false));
	return dgds;
}

template<typename T> DiscreteGraphicalParametric* DiscreteGraphicalReestimation<T>::type_parametric_estimation(int algo, int max_step) const
{
	UndirectedGraph ugraph = undirected_graph_computation(graph_estimation(algo, max_step));
	DistributionCollection distribution_collection = distribution_collection_computation(ugraph, true);
	BidirectedGraph bgraph = bidirected_graph_computation(distribution_collection.decomposition);
	DiscreteGraphicalParametric* dgd = new DiscreteGraphicalParametric(nb_variables, bgraph, distribution_collection);
	return dgd;
}

template<typename T> DiscreteGraphicalParametric* DiscreteGraphicalReestimation<T>::type_parametric_estimation(int algo, double min_threshold) const
{
	UndirectedGraph ugraph = undirected_graph_computation(graph_estimation(algo, min_threshold));
	DistributionCollection distribution_collection = distribution_collection_computation(ugraph, true);
	BidirectedGraph bgraph = bidirected_graph_computation(distribution_collection.decomposition);
	DiscreteGraphicalParametric* dgd = new DiscreteGraphicalParametric(nb_variables, bgraph, distribution_collection);
	return dgd;
}

/*
template<typename T> std::vector< DiscreteGraphicalDistribution* > DiscreteGraphicalReestimation<T>::stepwise_distribution_estimation(int algo) const
{
	std::vector< UndirectedGraph > ugraphs = undirected_graphs_computation(graph_estimation(algo));
	std::vector< DiscreteGraphicalDistribution* > dgds(ugraphs.size());
	for(std::vector< UndirectedGraph >::iterator it = ugraphs.begin(); it != ugraphs.end(); ++it)
		dgds[distance(ugraphs.begin(), it)] = new DiscreteGraphicalDistribution(nb_variables, *it, distribution_estimation(get_reestimation_collections(*it)));
	return dgds;
}
*/
template<typename T> EstimationGraph DiscreteGraphicalReestimation<T>::algorithm_ML(int max_step) const
{
	EdgesScoresML edges_scores = initialize_scores_ML();
	EstimationGraph egraph(nb_variables);
	int step = 1;
	while(continue_ML(edges_scores, egraph, max_step, step))
		execute_algorithm_ML(edges_scores, egraph, step);
	return egraph;
}

template<typename T> EstimationGraph DiscreteGraphicalReestimation<T>::algorithm_CLT(int max_step) const
{
	EdgesScores edges_scores = initialize_scores(CLT);
	EstimationGraph egraph(nb_variables);
	int step = 1;
	while(continue_CLT(edges_scores, egraph, max_step, step))
		execute_algorithm_CLT(edges_scores, egraph, step);
	return egraph;
}

template<typename T> EstimationGraph DiscreteGraphicalReestimation<T>::algorithm_MRNET(double min_threshold) const
{
	EdgesScoresMRNET edges_scores = initialize_scores_MRNET();
	EstimationGraph egraph(nb_variables);
	int step = 1;
	while(continue_MRNET(edges_scores, egraph, min_threshold))
		execute_algorithm_MRNET(edges_scores, egraph, step);
	return egraph;
}

template<typename T> EstimationGraph DiscreteGraphicalReestimation<T>::algorithms_threshold(int algo, double min_threshold) const
{
	EdgesScores edges_scores = initialize_scores(algo);
	EstimationGraph egraph(nb_variables);
	int step = 1;
	while(continue_threshold_algorithms(edges_scores, egraph, min_threshold))
		execute_algorithms(edges_scores, egraph, step);
	return egraph;
}

template<typename T> EdgesScoresML DiscreteGraphicalReestimation<T>::initialize_scores_ML() const
{
	EdgesScoresML edges_scores;
	Edges edges = get_edges();
	for(Edges::iterator it = edges.begin(); it != edges.end(); ++it)
		edges_scores.insert(EdgeScoreML(compute_conditional_mutual_information(*it, Neighborhood()), std::pair<Edge, Neighborhood>(*it, Neighborhood())));
	return edges_scores;
}

template<typename T> EdgesScoresMRNET DiscreteGraphicalReestimation<T>::initialize_scores_MRNET() const
{
	EdgesScoresMRNET edges_scores;
	Edges edges = get_edges();
	double score;
	for(Edges::iterator it = edges.begin(); it != edges.end(); ++it)
	{
		score = compute_conditional_mutual_information(*it, Neighborhood());
		edges_scores.insert(EdgeScoreMRNET(score, std::pair<Edge, double>(*it, score)));
	}
	return edges_scores;
}

template<typename T> EdgesScores DiscreteGraphicalReestimation<T>::initialize_scores(int algo) const
{
	EdgesScores edges_scores;
	Edges edges = get_edges();
	switch(algo)
	{
		case CLT :
		case RELEVANCE :
			for(Edges::iterator it = edges.begin(); it != edges.end(); ++it)
				edges_scores.insert(EdgeScore(compute_conditional_mutual_information(*it, Neighborhood()), *it));
			break;
		case CLR :
			std::map<unsigned int, std::vector<double> > vertices_mi;
			std::map<unsigned int, std::vector<double> >::iterator itb;
			double score;
			std::map<Edge, double> tedges_scores;
			for(Edges::iterator it = edges.begin(); it != edges.end(); ++it)
			{
				score = compute_conditional_mutual_information(*it, Neighborhood());
				tedges_scores.insert(std::pair<Edge, double>(*it, score));
				itb = vertices_mi.find((*it).first);
				if(itb == vertices_mi.end())
				{
					vertices_mi.insert(std::pair<unsigned int, std::vector<double> >((*it).first, std::vector<double>(1, score)));
				} else {
					((*itb).second).push_back(score);
				}
				itb = vertices_mi.find((*it).second);
				if(itb == vertices_mi.end())
				{
					vertices_mi.insert(std::pair<unsigned int, std::vector<double> >((*it).second, std::vector<double>(1, score)));
				} else {
					((*itb).second).push_back(score);
				}
			}
			double m0, m1, v0, v1;
			for(Edges::iterator it = edges.begin(); it != edges.end(); ++it)
			{
				score = (*(tedges_scores.find((*it)))).second;
				itb = vertices_mi.find((*it).first);
				m0 = 0;
				for(std::vector<double>::iterator itt = ((*itb).second).begin(); itt != ((*itb).second).end(); ++itt)
				 m0 += *itt;
				m0 /= ((*itb).second).size();
				v0 = 0;
				for(std::vector<double>::iterator itt = ((*itb).second).begin(); itt != ((*itb).second).end(); ++itt)
				 v0 += pow(*itt - m0,2);
				v0 /= ((*itb).second).size()-1;
				itb = vertices_mi.find((*it).second);
				m1 = 0;
				for(std::vector<double>::iterator itt = ((*itb).second).begin(); itt != ((*itb).second).end(); ++itt)
				 m1 += *itt;
				m1 /= ((*itb).second).size();
				v1 = 0;
				for(std::vector<double>::iterator itt = ((*itb).second).begin(); itt != ((*itb).second).end(); ++itt)
				 v1 += pow(*itt - m1,2);
				v1 /= ((*itb).second).size()-1;
				edges_scores.insert(EdgeScore(sqrt(pow(MAX(score-m0,0),2)/v0+pow(MAX(score-m1,0),2)/v1),*it));
			}
			break;
	}
	return edges_scores;
}

template<typename T> bool DiscreteGraphicalReestimation<T>::continue_CLT(const EdgesScores& edges_scores, const EstimationGraph& egraph, int max_step, int step) const
{
	if(boost::num_edges(egraph) == nb_variables-1)
	{
		return false;
	} else {
		if(max_step >= 0 && step > max_step)
		{
			return false;
		} else {
			return true;
			/*if(max_step >= 0 && nb_elements*((*(edges_scores.rbegin())).first) <= 10e-4)
			{
				return false;
			} else {
				return true;
			}*/
		}
	}
}

template<typename T> bool DiscreteGraphicalReestimation<T>::continue_ML(const EdgesScoresML& edges_scores, const EstimationGraph& egraph, int max_step, int step) const
{
	if(boost::num_edges(egraph) == nb_variables*(nb_variables-1)/2.)
	{
		return false;
	} else {
		if(max_step >= 0 && step > max_step)
		{
			return false;
		} else {
			return true;
			/*if(max_step >= 0 && ((*(edges_scores.rbegin())).first)*nb_elements <= 10e-4)
			{
				return false;
			} else {
				return true;
			}*/
		}
	}
}

template<typename T> bool DiscreteGraphicalReestimation<T>::continue_MRNET(const EdgesScoresMRNET& edges_scores, const EstimationGraph& egraph, double min_threshold) const
{
	if(boost::num_edges(egraph) == nb_variables*(nb_variables-1)/2.)
	{
		return false;
	} else {
		if((*(edges_scores.rbegin())).first < min_threshold)
		{
			return false;
		}	else {
			return true;
		}
	}
}

template<typename T> bool DiscreteGraphicalReestimation<T>::continue_threshold_algorithms(const EdgesScores& edges_scores, const EstimationGraph& egraph, double min_threshold) const
{
	if(boost::num_edges(egraph) == nb_variables*(nb_variables-1)/2.)
	{
		return false;
	} else {
		if((*(edges_scores.rbegin())).first < min_threshold)
		{
			return false;
		}	else {
			return true;
		}
	}
}

template<typename T> void DiscreteGraphicalReestimation<T>::execute_algorithm_CLT(EdgesScores& edges_scores, EstimationGraph& egraph, int& step) const
{
	EstimationGraphTraits::vertices_size_type d[boost::num_vertices(egraph)];
	std::fill_n(d, boost::num_vertices(egraph), 0);
	EdgesScores::iterator it = edges_scores.end();
	--it;
	boost::breadth_first_search(egraph, ((*it).second).first, boost::visitor(boost::make_bfs_visitor(boost::record_distances(d, boost::on_tree_edge()))));
	while(d[((*it).second).second] != 0)
	{
		std::fill_n(d, boost::num_vertices(egraph), 0);
		edges_scores.erase(it);
		it = edges_scores.end();
		--it;
		boost::breadth_first_search(egraph, ((*it).second).first, boost::visitor(boost::make_bfs_visitor(boost::record_distances(d, boost::on_tree_edge()))));
	}
	estimated_edge ep;
	ep.step = step;
	ep.score = (*it).first;
	boost::add_edge(((*it).second).first, ((*it).second).second, ep, egraph);
	edges_scores.erase(it);
	step++;
}

template<typename T> void DiscreteGraphicalReestimation<T>::execute_algorithm_ML(EdgesScoresML& edges_scores, EstimationGraph& egraph, int& step) const
{
	EdgesScoresML t_edges_scores = edges_scores;
	edges_scores.clear();
	//EdgesScoresML::iterator it;
	EdgeScoreML es;
	Neighborhood n;
	bool updated;
	estimated_edge ep;
	for(EdgesScoresML::iterator it = t_edges_scores.begin(); it != t_edges_scores.end(); ++it)
	{
		es = *it;
		n = get_neighborhood((es.second).first, egraph);
		if(same_neighborhood(n, ((es).second).second))
		{
			edges_scores.insert(es);
		} else {
			(es.second).second = n;
			es.first = compute_conditional_mutual_information((es.second).first, (es.second).second);
			edges_scores.insert(es);
		}
	}
	EdgesScoresML::iterator it = edges_scores.end();
	it--;
	ep.score = (*it).first;
	ep.step = step;
	boost::add_edge((((*it).second).first).first, (((*it).second).first).second, ep, egraph);
	edges_scores.erase(it);
	/*if(step > 2 && step < nb_variables*(nb_variables-1)/2.)
	{
		updated = false;
		while(!updated)
		{
			it = edges_scores.end();
			--it;
			es = *it;
			edges_scores.erase(it);
			n = get_neighborhood((es.second).first, egraph);
			if(same_neighborhood(n, ((es).second).second))
			{
				updated = true;
				ep.score = es.first;
				ep.step = step;
				boost::add_edge(((es.second).first).first, ((es.second).first).second, ep, egraph);
				std::cout << ep.score << std::endl;
			} else {
				(es.second).second = n;
				es.first = compute_conditional_mutual_information((es.second).first, (es.second).second);
				if(es.first >= (*(edges_scores.rbegin())).first)
				{
					updated = true;
					ep.score = es.first;
					ep.step = step;
					boost::add_edge(((es.second).first).first, ((es.second).first).second, ep, egraph);
					std::cout << ep.score << std::endl;
				} else {
					edges_scores.insert(es);
				}
			}
		}
	} else {
		it = edges_scores.end();
		--it;
		ep.score = es.first;
		ep.step = step;
		boost::add_edge((((*it).second).first).first, (((*it).second).first).second, ep, egraph);
		edges_scores.erase(it);
	}*/
	step++;
}

template<typename T> bool DiscreteGraphicalReestimation<T>::same_neighborhood(const Neighborhood& n0, const Neighborhood& n1) const
{
	if(n0.size() != n1.size())
	{
		return false;
	} else {
		if(n0.size() == 0)
		{
			return true;
		} else {
			bool end = (n0.size() > 0);
			bool one_different = false;
			if(!end)
			{
				Neighborhood::iterator it0 = n0.begin(), it1 = n1.begin();
				while(!end && !one_different)
				{
					if(*it0 != *it1)
						one_different = true;
					++it0;
					++it1;
					if(it0 == n0.end() || it1 == n1.end())
						end = true;
				}
			}
			return !(one_different);
		}
	}
}

template<typename T> void DiscreteGraphicalReestimation<T>::execute_algorithm_MRNET(EdgesScoresMRNET& edges_scores, EstimationGraph& egraph, int& step) const
{
	EdgesScoresMRNET::iterator it;
	EdgeScoreMRNET es;
	bool updated = false;
	EstimationGraphTraits::adjacency_iterator ai, ai_end;
	double fscore, sscore;
	estimated_edge ep;
	ep.step = step;
	if(step < nb_variables*(nb_variables-1)/2.)
	{
		while(!updated)
		{
			it = edges_scores.end();
			--it;
			es = *it;
			edges_scores.erase(it);
			ep.score =  (es.second).second;
			fscore = ep.score;
			for(boost::tie(ai, ai_end) = boost::adjacent_vertices(((es.second).first).first, egraph); ai != ai_end; ++ai)
				fscore -= (egraph[boost::edge(*ai, ((es.second).first).first, egraph).first]).score;
			sscore = ep.score;
			for(boost::tie(ai, ai_end) = boost::adjacent_vertices(((es.second).first).second, egraph); ai != ai_end; ++ai)
				sscore -= (egraph[boost::edge(*ai, ((es.second).first).second, egraph).first]).score;
			es.first = MAX(fscore, sscore);
			if(es.first >= (*(edges_scores.rbegin())).first)
			{
				updated = true;
				boost::add_edge((((es).second).first).first, (((es).second).first).second, ep, egraph);
			} else {
				edges_scores.insert(es);
			}
		}
	} else {
		it = edges_scores.end();
		--it;
		es = *it;
		edges_scores.erase(it);
		boost::add_edge((((es).second).first).first, (((es).second).first).second, ep, egraph);
	}
	step++;
}

template<typename T> void DiscreteGraphicalReestimation<T>::execute_algorithms(EdgesScores& edges_scores, EstimationGraph& egraph, int& step) const
{
	EdgesScores::iterator it;
	std::pair< EdgesScores::iterator, EdgesScores::iterator > pit;
	estimated_edge ep;
	pit = edges_scores.equal_range((*(edges_scores.rbegin())).first);
	for(it = pit.first; it != pit.second; ++it)
	{
		ep.step = step;
		ep.score = (*it).first;
		boost::add_edge(((*it).second).first, ((*it).second).second, ep, egraph);
	}
	edges_scores.erase(pit.first, pit.second);
	step++;
}

template<typename T> Neighborhood DiscreteGraphicalReestimation<T>::get_neighborhood(const Edge& edge, const EstimationGraph& egraph) const
{
	Neighborhood n, nt;
	EstimationGraphTraits::adjacency_iterator ai, ai_end;
	for(boost::tie(ai, ai_end) = boost::adjacent_vertices(edge.first, egraph); ai != ai_end; ++ai)
		nt.insert(*ai);
	for(boost::tie(ai, ai_end) = boost::adjacent_vertices(edge.second, egraph); ai != ai_end; ++ai)
	{
		if(nt.find(*ai) != nt.end())
			n.insert(*ai);
	}
	return n;
}

template<typename T> double DiscreteGraphicalReestimation<T>::compute_conditional_mutual_information(const Edge& edge, const Neighborhood& neighborhood) const
{
	DiscreteUnivariateConditionalReestimation<T> *hf = new DiscreteUnivariateConditionalReestimation<T>(histogram, neighborhood, edge.first);
	double cmi = hf->entropy_computation();
	DiscreteUnivariateConditionalReestimation<T> *hs = new DiscreteUnivariateConditionalReestimation<T>(histogram, neighborhood, edge.second);
	cmi += hs->entropy_computation();
	Neighborhood variables;
	variables.insert(edge.first);
	variables.insert(edge.second);
	DiscreteMultivariateConditionalReestimation<T> *h = new DiscreteMultivariateConditionalReestimation<T>(histogram, neighborhood, variables);
	cmi -= h->entropy_computation();
	delete hf;
	delete hs;
	delete h;
	return cmi;
}

template<typename T> Edges DiscreteGraphicalReestimation<T>::get_edges() const
{
	Edges edges;
	for(unsigned int i = 0; i < nb_variables-1; ++i)
	{
		for(unsigned int j = i+1; j < nb_variables; ++j)
			edges.push_back(Edge(i,j));
	}
	return edges;
}

template<typename T> UndirectedGraph DiscreteGraphicalReestimation<T>::undirected_graph_computation(const EstimationGraph& egraph) const
{
	UndirectedGraph ugraph(nb_variables);
	EstimationGraphTraits::edge_iterator ei, ei_end;
	for(boost::tie(ei, ei_end) = boost::edges(egraph); ei != ei_end; ++ei)
		boost::add_edge(boost::source(*ei, egraph), boost::target(*ei, egraph), ugraph);
	return ugraph;
}

template<typename T> void DiscreteGraphicalReestimation<T>::bron_kerbosch_all_cliques(const UndirectedGraph& ugraph, const Neighborhood& r, Neighborhood p, Neighborhood x, Cliques& cliques, const int& superior_to) const
{
	if(p.size() == 0 && x.size() == 0)
	{
		if(r.size() > superior_to)
			cliques.insert(r);
	} else {
		//unsigned int u = *(get_intersection(p, x).begin());
		Neighborhood rold, rnew, pnew, xnew,	vertices = p;// get_setminus(p, get_neighborhood(u, ugraph));
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
			bron_kerbosch_all_cliques(ugraph, rnew, pnew, xnew, cliques, superior_to);
			x.insert(*it);
		}
	}
}

template<typename T> std::vector< Neighborhood > DiscreteGraphicalReestimation<T>::get_ps(const UndirectedGraph& ugraph) const
{
	std::vector<int> component(nb_variables);
	int num = boost::connected_components(ugraph, &component[0]);
	std::vector< Neighborhood > ps = std::vector< Neighborhood >(num, Neighborhood());
	for(unsigned int i = 0; i < nb_variables; i++)
		ps[component[i]].insert(i);
	return ps;
}

template<typename T> std::vector< std::pair<Neighborhood, Neighborhood> > DiscreteGraphicalReestimation<T>::get_component_decomposition(const Neighborhood& start, const Neighborhood& vertices, const UndirectedGraph& ugraph) const
{
	std::vector< std::pair<Neighborhood, Neighborhood> > decomposition = std::vector< std::pair<Neighborhood, Neighborhood> >(1, std::pair<Neighborhood, Neighborhood>(Neighborhood(), start));
	Neighborhood front_vertices = start, t_front_vertices, setted_vertices = start, temp, i;
	Cliques cliques;
	while(setted_vertices.size() != vertices.size())
	{
		t_front_vertices.clear();
		for(Neighborhood::iterator it = front_vertices.begin(); it != front_vertices.end(); ++it)
		{
			temp = get_setminus(get_neighborhood(*it, ugraph), setted_vertices);
			t_front_vertices.insert(temp.begin(), temp.end());
			for(Neighborhood::iterator itb = temp.begin(); itb != temp.end(); ++itb)
			{
				i.clear();
				i.insert(*itb);
				decomposition.push_back(std::pair<Neighborhood, Neighborhood>(get_intersection(get_neighborhood(*itb, ugraph), setted_vertices), i));
				setted_vertices.insert(*itb);
			}
		
		}
		front_vertices = t_front_vertices;
	}
	return decomposition;
}

/*template<typename T> ReestimationCollection<T> DiscreteGraphicalReestimation<T>::get_decomposition_reestimation_collection(const std::vector<std::pair<Neighborhood, Neighborhood> >& decomposition) const
{
	ReestimationCollection<T> reestimation_collection;
	for(std::vector< std::pair< Neighborhood, Neighborhood > >::const_iterator it = decomposition.begin(); it != decomposition.end(); ++it)
	{
		(reestimation_collection.decomposition).push_back(*it);
		if(((*it).first).size() == 0)
		{
			if(((*it).second).size() == 1)
				(reestimation_collection.DUR).push_back(get_reestimation(*(((*it).second).begin())));
			else
				(reestimation_collection.DMR).push_back(get_reestimation((*it).second));
		} else {
			if(((*it).second).size() == 1)
				(reestimation_collection.DUCR).push_back(get_reestimation((*it).first, *(((*it).second).begin())));
			else
				(reestimation_collection.DMCR).push_back(get_reestimation((*it).first, (*it).second));
		}
	}
	return reestimation_collection;
}*/

template<typename T> Neighborhood DiscreteGraphicalReestimation<T>::get_neighborhood(const Neighborhood& vertices, const UndirectedGraph& ugraph) const
{
	Neighborhood neighborhood;
	UndirectedGraphTraits::adjacency_iterator ai, ai_end;
	for(Neighborhood::const_iterator it = vertices.begin(); it != vertices.end(); ++it)
	{
		for(boost::tie(ai, ai_end) = boost::adjacent_vertices(*it, ugraph); ai != ai_end; ++ai)
		{
			if(vertices.find(*ai) == vertices.end())
				neighborhood.insert(*ai);
		}
	}
	return neighborhood;
}

template<typename T> Neighborhood DiscreteGraphicalReestimation<T>::get_neighborhood(unsigned int vertex, const UndirectedGraph& ugraph) const
{
	Neighborhood neighborhood;
	UndirectedGraphTraits::adjacency_iterator ai, ai_end;
	for(boost::tie(ai, ai_end) = boost::adjacent_vertices(vertex, ugraph); ai != ai_end; ++ai)
		neighborhood.insert(*ai);
	return neighborhood;
}

template<typename T> Neighborhood DiscreteGraphicalReestimation<T>::get_intersection(const Neighborhood& vertices0, const Neighborhood& vertices1) const
{
	Neighborhood intersection;
	for(Neighborhood::const_iterator it = vertices0.begin(); it != vertices0.end(); ++it)
	{
		if(vertices1.find(*it) != vertices1.end())
			intersection.insert(*it);
	}
	return intersection;
}

template<typename T> Neighborhood DiscreteGraphicalReestimation<T>::get_complement(const Neighborhood& intersection, const Neighborhood& vertices) const
{
	Neighborhood complement;
	for(Neighborhood::const_iterator it = vertices.begin(); it != vertices.end(); ++it)
	{
		if(intersection.find(*it) == intersection.end())
			complement.insert(*it);
	}
	return complement;
}

template<typename T>  Neighborhood DiscreteGraphicalReestimation<T>::get_setminus(const Neighborhood& set, const Neighborhood& minus) const
{
	Neighborhood setminus = set;
	for(Neighborhood::const_iterator it = minus.begin(); it != minus.end(); ++it)
		setminus.erase(*it);
	return setminus;
}

template<typename T>  Neighborhood DiscreteGraphicalReestimation<T>::get_setminus(const Neighborhood& set, unsigned int minus) const
{
	Neighborhood setminus = set;
	setminus.erase(minus);
	return setminus;
}

template<typename T>  Neighborhood DiscreteGraphicalReestimation<T>::get_setplus(const Neighborhood& set, unsigned int plus) const
{
	Neighborhood setplus = set;
	setplus.insert(plus);
	return setplus;
}

template<typename T> Reestimation<T>* DiscreteGraphicalReestimation<T>::get_reestimation(unsigned int vertex) const
{
	std::map<int, T> marginal_histogram;
	typename std::map<int, T>::iterator itm;
	for(typename std::map< std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
	{
		itm = marginal_histogram.find(((*it).first)[vertex]);
		if(itm == marginal_histogram.end())
			marginal_histogram.insert(std::pair<int, T>(((*it).first)[vertex], (*it).second));
		else
			(*itm).second += (*it).second;
	}
	Reestimation<T> *reestimation = new Reestimation<T>((*(marginal_histogram.rbegin())).first+1);
	for(typename std::map<int, T>::iterator itm = marginal_histogram.begin(); itm != marginal_histogram.end(); ++itm)
		reestimation->frequency[(*itm).first] = (*itm).second;
	reestimation->nb_element_computation();
  reestimation->offset_computation();
  reestimation->max_computation();
  reestimation->mean_computation();
  reestimation->variance_computation();
	return reestimation;
}

template<typename T> Stat_trees::DiscreteMultivariateReestimation<T>* DiscreteGraphicalReestimation<T>::get_reestimation(const Neighborhood& vertices) const
{
	std::map<std::vector<unsigned int>, T> marginal_histogram;
	int min = I_MAX;
	int max = I_MIN;
	std::vector<unsigned int> new_vector(vertices.size());
	typename std::map< std::vector<unsigned int>, T>::iterator itm;
	for(typename std::map< std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
	{
		for(Neighborhood::const_iterator itb = vertices.begin(); itb != vertices.end(); ++itb)
			new_vector[distance(vertices.begin(), itb)] = ((*it).first)[*itb];
		min = MIN(min, *(min_element(new_vector.begin(), new_vector.end())));
		max = MAX(max, *(max_element(new_vector.begin(), new_vector.end())));
		itm = marginal_histogram.find(new_vector);
		if(itm == marginal_histogram.end())
			marginal_histogram.insert(std::pair< std::vector<unsigned int>, T>(new_vector, (*it).second));
		else
			(*itm).second += (*it).second;
	}
	Stat_trees::DiscreteMultivariateReestimation<T> *reestimation = new Stat_trees::DiscreteMultivariateReestimation<T>(vertices.size(), min, max-min+1);
	for(typename std::map< std::vector<unsigned int>, T>::iterator it = marginal_histogram.begin(); it != marginal_histogram.end(); ++it)
		reestimation->update((*it).first, (*it).second);
	reestimation->nb_value_computation();
  reestimation->offset_computation();
  reestimation->max_computation();
  reestimation->mean_computation();
  reestimation->variance_computation();
	reestimation->sum_computation();
	return reestimation;
}

template<typename T> DiscreteUnivariateConditionalReestimation<T>* DiscreteGraphicalReestimation<T>::get_reestimation(const Neighborhood& covertices, unsigned int vertex) const
{
	return new DiscreteUnivariateConditionalReestimation<T>(histogram, covertices, vertex);
}

template<typename T> DiscreteMultivariateConditionalReestimation<T>* DiscreteGraphicalReestimation<T>::get_reestimation(const Neighborhood& covertices, const Neighborhood& vertices) const
{
	return new DiscreteMultivariateConditionalReestimation<T>(histogram, covertices, vertices);
}
/*
template<typename T> DistributionCollection DiscreteGraphicalReestimation<T>::distribution_estimation(const std::vector< ReestimationCollection<T> >& reestimation_collections) const
{
	DistributionCollection distribution_collection;
	(distribution_collection.decomposition).resize(0);
	(distribution_collection.DUD).resize(0);
	(distribution_collection.DUCD).resize(0);
	(distribution_collection.DMD).resize(0);
	(distribution_collection.DMCD).resize(0);
	for(typename std::vector< ReestimationCollection<T> >::const_iterator it = reestimation_collections.begin(); it != reestimation_collections.end(); ++it)
	{
		for(std::vector< std::pair<Neighborhood, Neighborhood> >::const_iterator itd = ((*it).decomposition).begin(); itd != ((*it).decomposition).end(); ++itd)
			(distribution_collection.decomposition).push_back((*itd));
		for(typename std::vector< Reestimation<T>* >::const_iterator itr = ((*it).DUR).begin(); itr != ((*it).DUR).end(); ++itr)
		{
			(distribution_collection.DUD).push_back(new Distribution((*itr)->nb_value));
			(*itr)->distribution_estimation((distribution_collection.DUD).back());
		}
		for(typename std::vector< Stat_trees::DiscreteMultivariateReestimation<T>* >::const_iterator itr = ((*it).DMR).begin(); itr != ((*it).DMR).end(); ++itr)
			(distribution_collection.DMD).push_back(new Stat_trees::DiscreteMultivariateDistribution(*((*itr)->distribution_estimation())));
		for(typename std::vector< DiscreteUnivariateConditionalReestimation<T>* >::const_iterator itr = ((*it).DUCR).begin(); itr != ((*it).DUCR).end(); ++itr)
			(distribution_collection.DUCD).push_back(new DiscreteUnivariateConditionalDistribution(*((*itr)->distribution_estimation())));
		for(typename std::vector< DiscreteMultivariateConditionalReestimation<T>* >::const_iterator itr = ((*it).DMCR).begin(); itr != ((*it).DMCR).end(); ++itr)
			(distribution_collection.DMCD).push_back(new DiscreteMultivariateConditionalDistribution(*((*itr)->distribution_estimation())));
	}
	return distribution_collection;
}
*/
template<typename T> double DiscreteGraphicalReestimation<T>::likelihood_computation(DiscreteGraphicalDistribution& dgd) const
{
	double lh = 0;
	for(typename std::map<std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
		lh += (*it).second*dgd.get_mass((*it).first, true);
	return lh;
}

template<typename T> double DiscreteGraphicalReestimation<T>::penalized_likelihood_computation(DiscreteGraphicalDistribution& dgd, int criterion) const
{
	double p = 2*likelihood_computation(dgd);
	unsigned int dim = dgd.get_nb_parameters();
	switch(criterion)
	{
		case AIC :
			p -= 2*dim;
			break;
		case BIC :
			p -= dim*log(nb_elements);
			break;
	}
	return p;
}

template<typename T> double DiscreteGraphicalReestimation<T>::kullback_distance_computation(DiscreteGraphicalDistribution& dgd) const
{
	double d = 0;
	for(typename std::map<std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
		d += (*it).second/((double)nb_elements)*( log((*it).second) - log(nb_elements) -dgd.get_mass((*it).first,true));
	return d;
}

template<typename T> std::vector< UndirectedGraph > DiscreteGraphicalReestimation<T>::undirected_graphs_computation(const EstimationGraph& egraph) const
{
	std::vector< UndirectedGraph > ugraphs;
	UndirectedGraph ugraph(nb_variables);
	EstimationGraphTraits::edge_iterator ei, ei_end;
	unsigned int step = 0;
	for(boost::tie(ei, ei_end) = boost::edges(egraph); ei != ei_end; ++ei)
	{
		if(egraph[*ei].step > step)
		{
			ugraphs.push_back(ugraph);
			step++;
		}
		boost::add_edge(boost::source(*ei, egraph), boost::target(*ei, egraph), ugraph);
	}
	ugraphs.push_back(ugraph);
	return ugraphs;
}

/*template<typename T> double DiscreteGraphicalReestimation<T>::get_mutual_information(const std::set<unsigned int>& vertices) const
{
}*/

/*template<typename T> ReestimationCollection DiscreteGraphicalReestimation<T>::decomposition_computation(const UndirectedGraph& ugraph) const
{
	BidirectedGraph bgraph = directed_graph_computation(ugraph,0);
	std::vector<unsigned int> vertices_order = topological_sort(bgraph, 0);
	Neighborhood neighborhood;
	for(std::vector<unsigned int>::iterator it = vertices_order.begin(); it != vertices_order.end(); ++it)
	{
		neighborhood = get_neighborhood(*it, bgraph);
		if(neighborhood.size() == 0)
		{
		} else {
		}
	}
	if(triangulated_graph(ugraph))
	{
		CliqueGraph ctree = clique_tree_computation(ugraph);
		return clique_tree_mass_computation(ctree);
	} else {
		return non_triangulated_mass_table_computation(ugraph);
	}
	return mass;
}*/
/*
template<typename T> std::map< std::vector<int>, T> get_marginal_histogram(const Neighborhood& vertices) const
{
	std::map< std::vector<int>, T> marginal_histogram;
	typename std::map<std::vector<int>, T>::iterator it_new;
	std::vector<int> new_vector(vertices.size());
	Neighborhood::const_iterator itv;
	for(typename std::map<std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
	{
		for(std::vector<int>::const_iterator itb = ((*it).first).begin(); itb != ((*it).first).end(); ++itb)
		{
			itv = vertices.find(distance(((*it).first).begin(), itb))
			if(itv != vertices.end())
				new_vector[distance(vertices.begin(), itv)] = *itb;
		}
		it_new = marginal_histogram.find(new_vector);
		if(it_new == marginal_histogram.end())
			marginal_histogram.insert(std::pair< std::vector<int>, T>(new_vector, (*it).second));
		else
			(*it_new).second += (*it).second;
	}
	return marginal_histogram;
}*/

/*template<typename T> BidirectedGraph DiscreteGraphicalReestimation<T>::directed_graph_computation(const UndirectedGraph& ugraph, int from_vertex)
{
	BidirectedGraph bgraph(nb_variables);
	Neighborhood unsetted_vertices, setted_vertices;
	Neighborhood::iterator itu, it;
	setted_vertices.insert(from_vertex);
	for(unsigned int i = 0; i < nb_variables; ++i)
	{
		if(i != from_vertex)
			unsetted_vertices.insert(i);
	}
	while(unsetted_vertices.size() != 0)
	{
		itu = unsetted_vertices.begin();
		for(it = setted_vertices.begin(); it != setted_vertices.end(); ++it)
		{
			if(boost::edge(*itu, *it, ugraph).second)
				boost::add_edge(*it, *itu, bgraph);
		}
		setted_vertices.insert(*itu);
		unsetted_vertices.erase(itu);
	}
	return bgraph;
}*/

/*template<typename T> Neighborhood DiscreteGraphicalReestimation<T>::get_neighborhood(unsigned int vertex, const BidirectedGraph& bgraph) const
{
	BidirectedGraphTraits::adjacency_iterator ai, ai_end;
	Neighborhooh neighborhood;
	for(boost::tie(ai, ai_end) = boost::adjacent_vertices(vertex, bgraph); ai != ei_end; ++ai)
		neighborhood.insert(*ai);
	return neighborhood;
}*/

/*template<> std::map< std::vector<int>, double > DiscreteGraphicalReestimation<unsigned int>::non_triangulated_mass_table_computation(const UndirectedGraph& ugraph)
{
	std::map<std::vector<int>, double> mass;
	for(std::map< std::vector<int>, unsigned int>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
	{
		mass.insert(std::pair<std::vector<int>, double>((*it).first, (double)((*it).second)));
	}
	return mass;
}

template<> std::map< std::vector<int>, double > DiscreteGraphicalReestimation<double>::non_triangulated_mass_table_computation(const UndirectedGraph& ugraph)
{
	return histogram;
}

template<typename T> bool DiscreteGraphicalReestimation<T>::triangulated_graph(const UndirectedGraph& ugraph) const
{
	UndirectedGraphTraits::vertex_iterator vi, vi_end;
	UndirectedGraphTraits::adjacency_iterator ai, ai_end;
	Neighborhood labelled_vertices, non_labelled_vertices, i_labelled_vertices, adj_labelled_vertices;
	Neighborhood::iterator itl0, itl1;
	labelled_vertices.insert(0);
	unsigned int max, i_max;
	bool success = true, end;
	for(unsigned int i = 1; i < nb_variables; ++i)
	{
		non_labelled_vertices.insert(i);
	}
	for(unsigned int i = 1; i < nb_variables; ++i)
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



template<typename T> CliqueGraph DiscreteGraphicalReestimation<T>::clique_tree_computation(const UndirectedGraph& ugraph)
{
	Cliques cliques, ccliques;
	Neighborhood separator, p;
	for(unsigned int i = 0; i < boost::num_vertices(ugraph); ++i)
		p.insert(i);
	bron_kerbosch_all_cliques(Neighborhood(), p, Neighborhood(), cliques);
	CliqueGraph cgraph;
	clique_vertex cv;
	separator_edge se;
	DiscreteMultivariateHistogram *dmh = NULL;
	for(Cliques::iterator it = cliques.begin(); it != cliques.end(); ++it)
	{
		cv.clique = (*it).second;
		boost::add_vertex(cv, cgraph);
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
		card_sep_score.erase(itcss);
	}
	return cgraph;
}

template<typename T> std::map< std::vector<int>, double> clique_tree_decomposition(const CliqueGraph& ctree, int from_vertex) const
{
}*/

/*
bool DiscreteGraphicalEstimation::TriangulatedGraphChecking(EstimationGraph& egraph, const Edge& edge) const
{
	boost::add_edge(edge.first, edge.second, egraph);
	EstimationGraphTraits::vertex_iterator vi, vi_end;
	EstimationGraphTraits::adjacency_iterator ai, ai_end;
	Neighborhood labelled_vertices, non_labelled_vertices, i_labelled_vertices, adj_labelled_vertices;
	Neighborhood::iterator itl0, itl1;
	labelled_vertices.insert(0);
	unsigned int max, i_max;
	bool success = true, end;
	for(unsigned int i = 1; i < boost::num_vertices(egraph); ++i)
	{
		non_labelled_vertices.insert(i);
	}
	for(unsigned int i = 1; i < boost::num_vertices(egraph); ++i)
	{
		max = 0;
		for(Neighborhood::iterator it = non_labelled_vertices.begin(); it != non_labelled_vertices.end(); ++it)
		{
			adj_labelled_vertices.clear();
			for(boost::tie(ai, ai_end) = boost::adjacent_vertices(*it, egraph); ai != ai_end; ++ai)
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
					if(!boost::edge(*itl0, *itl1, egraph).second)
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
			{
				boost::remove_edge(edge.first, edge.second, egraph);
				return false;
			}
		}
	}
	boost::remove_edge(edge.first, edge.second, egraph);
	return success;
}
*/

template<typename T> DistributionCollection DiscreteGraphicalReestimation<T>::distribution_collection_computation(const UndirectedGraph& ugraph, bool parametric) const
{
	std::vector< Neighborhood > ps = get_ps(ugraph);
	std::vector< DistributionCollection > distribution_collections(ps.size());
	Cliques cliques;
	std::vector< std::vector< std::pair< Neighborhood, Neighborhood> > > decompositions;
	std::pair<double, DistributionCollection> distribution_collection, b_distribution_collection;
	for(std::vector<Neighborhood>::iterator it = ps.begin(); it != ps.end(); ++it)
	{
		cliques.clear();
		decompositions.clear();
		b_distribution_collection.first = -1*D_INF;
		bron_kerbosch_all_cliques(ugraph, Neighborhood(), *it, Neighborhood(), cliques);
		component_decompositions_computation((*it).size(), Neighborhood(), cliques, std::vector< std::pair<Neighborhood, Neighborhood> >(), ugraph, decompositions);
		for(std::vector< std::vector< std::pair<Neighborhood, Neighborhood> > >::iterator itb = decompositions.begin(); itb != decompositions.end(); ++itb)
		{
			Neighborhood::const_iterator ni_end;
			for(std::vector< std::pair<Neighborhood, Neighborhood> >::const_iterator itt = (*itb).begin(); itt != (*itb).end(); ++itt)
			{
				std::cout << " P[";
				ni_end = ((*itt).second).end();
				if(ni_end != ((*itt).second).begin())
				{
					ni_end--;
					for(Neighborhood::const_iterator ni = ((*itt).second).begin(); ni != ni_end; ++ni)
						std::cout << "v" << toString(*ni) << ", ";
					if((*itt).first.size() == 0)
					{
						std::cout << "v" << toString(*ni_end) << "]";
					} else {
						std::cout << "v" << toString(*ni_end) << " | ";
						ni_end = ((*itt).first).end();
						ni_end--;
						for(Neighborhood::const_iterator ni = ((*itt).first).begin(); ni != ni_end; ++ni)
							std::cout << "v" << toString(*ni) << ", ";
						std::cout << "v" << toString(*ni_end) << "]";
					}
				} else {
					std::cout << "]";
				}
			}
			std::cout << std::endl;
			distribution_collection = distribution_estimation(decomposition_reestimation_collection_computation(*itb), parametric);
			if(distribution_collection.first >= b_distribution_collection.first)
			{
				b_distribution_collection.first = distribution_collection.first;
				b_distribution_collection.second = distribution_collection.second;
			}
		}
		distribution_collections[distance(ps.begin(), it)] = b_distribution_collection.second;
	}
	return distribution_collections_agglomerate(distribution_collections);
}

template<typename T> void DiscreteGraphicalReestimation<T>::component_decompositions_computation(unsigned int component_nb_variables, const Neighborhood& setted_vertices, const Cliques& cliques, const std::vector< std::pair<Neighborhood, Neighborhood> >& decomposition, const UndirectedGraph& ugraph, std::vector< std::vector< std::pair<Neighborhood, Neighborhood> > >& decompositions) const
{
	if(setted_vertices.size() == component_nb_variables)
	{
		decompositions.push_back(decomposition);
	} else {
		Neighborhood intersection, variables, covariables, it_setted_vertices;
		Cliques it_cliques;
		Cliques::iterator itc;
		std::vector< std::pair<Neighborhood, Neighborhood> > it_decomposition;
		for(Cliques::const_iterator it = cliques.begin(); it != cliques.end(); ++it)
		{
			intersection = get_intersection(setted_vertices, *it);
			if(intersection.size() > 0 || setted_vertices.size() == 0)
			{
				variables = get_complement(intersection, *it);
				if(variables.size() > 0)
				{
					covariables = get_intersection(get_neighborhood(variables, ugraph), setted_vertices);
					it_setted_vertices = setted_vertices;
					it_cliques = cliques;
					it_decomposition = decomposition;
					it_decomposition.push_back(std::pair<Neighborhood, Neighborhood>(covariables, variables));
					it_setted_vertices.insert(variables.begin(), variables.end());
					itc = it_cliques.begin();
					advance(itc, distance(cliques.begin(), it));
					it_cliques.erase(itc);
					component_decompositions_computation(component_nb_variables, it_setted_vertices, it_cliques, it_decomposition, ugraph, decompositions);
				}
			}
		}
	}
}

template<typename T> ReestimationCollection<T> DiscreteGraphicalReestimation<T>::decomposition_reestimation_collection_computation(const std::vector< std::pair<Neighborhood, Neighborhood> >& decomposition) const
{
	ReestimationCollection<T> reestimation_collection;
	for(std::vector< std::pair< Neighborhood, Neighborhood > >::const_iterator it = decomposition.begin(); it != decomposition.end(); ++it)
	{
		(reestimation_collection.decomposition).push_back(*it);
		if(((*it).first).size() == 0)
		{
			if(((*it).second).size() == 1)
				(reestimation_collection.DUR).push_back(get_reestimation(*(((*it).second).begin())));
			else
				(reestimation_collection.DMR).push_back(get_reestimation((*it).second));
		} else {
			if(((*it).second).size() == 1)
				(reestimation_collection.DUCR).push_back(get_reestimation((*it).first, *(((*it).second).begin())));
			else
				(reestimation_collection.DMCR).push_back(get_reestimation((*it).first, (*it).second));
		}
	}
	return reestimation_collection;
}

template<typename T> std::pair<double, DistributionCollection > DiscreteGraphicalReestimation<T>::distribution_estimation(const ReestimationCollection<T>& reestimation_collection, bool parametric) const
{
	std::pair<double, DistributionCollection> distribution_collection;
	distribution_collection.first = 0;
	distribution_collection.second.decomposition = reestimation_collection.decomposition;
	(distribution_collection.second.DUD).resize(reestimation_collection.DUR.size());
	for(typename std::vector< Reestimation<T>* >::const_iterator it = reestimation_collection.DUR.begin(); it != reestimation_collection.DUR.end(); ++it)
	{	
		if(!parametric)
		{
			distribution_collection.second.DUD[distance(reestimation_collection.DUR.begin(), it)] = new Distribution((*it)->nb_value);
			(*it)->distribution_estimation(distribution_collection.second.DUD[distance(reestimation_collection.DUR.begin(), it)]);
		} else {
			distribution_collection.second.DUD[distance(reestimation_collection.DUR.begin(), it)] = new DiscreteParametric(*((*it)->type_parametric_estimation()));
		}
		std::cout << (*it)->likelihood_computation(*(distribution_collection.second.DUD[distance(reestimation_collection.DUR.begin(), it)])) << " ";
		distribution_collection.first += (*it)->likelihood_computation(*(distribution_collection.second.DUD[distance(reestimation_collection.DUR.begin(), it)]));
	}
	(distribution_collection.second.DMD).resize(reestimation_collection.DMR.size());
	StatError error;
	double l;
	for(typename std::vector< Stat_trees::DiscreteMultivariateReestimation<T>* >::const_iterator it = reestimation_collection.DMR.begin(); it != reestimation_collection.DMR.end(); ++it)
	{
		if(!parametric)
		{
			distribution_collection.second.DMD[distance(reestimation_collection.DMR.begin(), it)] = new Stat_trees::DiscreteMultivariateDistribution(*((*it)->distribution_estimation()));
		} else {
			distribution_collection.second.DMD[distance(reestimation_collection.DMR.begin(), it)] = new Stat_trees::DiscreteMultivariateParametric(*((*it)->type_parametric_estimation(error, l)));
		}
		std::cout << (*it)->likelihood_computation(*(distribution_collection.second.DMD[distance(reestimation_collection.DMR.begin(), it)])) << " ";
		distribution_collection.first += (*it)->likelihood_computation(*(distribution_collection.second.DMD[distance(reestimation_collection.DMR.begin(), it)]));
	}
	(distribution_collection.second.DUCD).resize(reestimation_collection.DUCR.size());
	for(typename std::vector< DiscreteUnivariateConditionalReestimation<T>* >::const_iterator it = reestimation_collection.DUCR.begin(); it != reestimation_collection.DUCR.end(); ++it)
	{
		if(!parametric)
		{
			distribution_collection.second.DUCD[distance(reestimation_collection.DUCR.begin(), it)] = new DiscreteUnivariateConditionalDistribution(*((*it)->distribution_estimation()));
		} else {
			distribution_collection.second.DUCD[distance(reestimation_collection.DUCR.begin(), it)] = new DiscreteUnivariateConditionalDistribution(*((*it)->type_parametric_estimation()));
		}
		distribution_collection.first += (*it)->likelihood_computation(*(distribution_collection.second.DUCD[distance(reestimation_collection.DUCR.begin(), it)]));
	}
	(distribution_collection.second.DMCD).resize(reestimation_collection.DMCR.size());
	for(typename std::vector< DiscreteMultivariateConditionalReestimation<T>* >::const_iterator it = reestimation_collection.DMCR.begin(); it != reestimation_collection.DMCR.end(); ++it)
	{
		if(!parametric)
		{
			distribution_collection.second.DMCD[distance(reestimation_collection.DMCR.begin(), it)] = new DiscreteMultivariateConditionalDistribution(*((*it)->distribution_estimation()));
		} else {
			distribution_collection.second.DMCD[distance(reestimation_collection.DMCR.begin(), it)] = new DiscreteMultivariateConditionalDistribution(*((*it)->type_parametric_estimation()));
		}
		distribution_collection.first += (*it)->likelihood_computation(*(distribution_collection.second.DMCD[distance(reestimation_collection.DMCR.begin(), it)]));
	}
	std::cout << distribution_collection.first << " end !" << std::endl;
	return distribution_collection;
}

template<typename T> DistributionCollection DiscreteGraphicalReestimation<T>::distribution_collections_agglomerate(const std::vector< DistributionCollection >& distribution_collections) const
{
	DistributionCollection distribution_collection;
	for(std::vector< DistributionCollection >::const_iterator it = distribution_collections.begin(); it != distribution_collections.end(); ++it)
	{
		for(std::vector< std::pair<Neighborhood, Neighborhood> >::const_iterator itb = (*it).decomposition.begin(); itb != (*it).decomposition.end(); ++itb)
			distribution_collection.decomposition.push_back(*itb);
		for(std::vector< Distribution* >::const_iterator itb = (*it).DUD.begin(); itb != (*it).DUD.end(); ++itb)
			distribution_collection.DUD.push_back(new Distribution(*(*itb)));
		for(std::vector< Stat_trees::DiscreteMultivariateDistribution* >::const_iterator itb = (*it).DMD.begin(); itb != (*it).DMD.end(); ++itb)
			distribution_collection.DMD.push_back(new Stat_trees::DiscreteMultivariateDistribution(*(*itb)));
		for(std::vector< DiscreteUnivariateConditionalDistribution* >::const_iterator itb = (*it).DUCD.begin(); itb != (*it).DUCD.end(); ++itb)
			distribution_collection.DUCD.push_back(new DiscreteUnivariateConditionalDistribution(*(*itb)));
		for(std::vector< DiscreteMultivariateConditionalDistribution* >::const_iterator itb = (*it).DMCD.begin(); itb != (*it).DMCD.end(); ++itb)
			distribution_collection.DMCD.push_back(new DiscreteMultivariateConditionalDistribution(*(*itb)));
	}
	return distribution_collection;
}

template<typename T> BidirectedGraph DiscreteGraphicalReestimation<T>::bidirected_graph_computation(const std::vector< std::pair<Neighborhood, Neighborhood> >& decomposition) const
{
	BidirectedGraph bgraph(nb_variables);
	//Neighborhood::const_iterator itn0_begin, itn0_end, itn1_begin, itn1_end;
	for(std::vector< std::pair<Neighborhood, Neighborhood> >::const_iterator it = decomposition.begin(); it != decomposition.end(); ++it)
	{
		if((*it).second.size() > 1)
		{
			for(Neighborhood::const_iterator itn0 = (*it).second.begin(); itn0 != (*it).second.end(); ++itn0)
			{
				for(Neighborhood::const_iterator itn1 = (*it).second.begin(); itn1 != (*it).second.end(); ++itn1)
				{
					if(itn1 != itn0)
						boost::add_edge(*itn0, *itn1, bgraph);
				}
			}
		}
		if((*it).first.size() > 0)
		{
			for(Neighborhood::const_iterator itn0 = (*it).first.begin(); itn0 != (*it).first.end(); ++itn0)
			{
				for(Neighborhood::const_iterator itn1 = (*it).second.begin(); itn1 != (*it).second.end(); ++itn1)
					boost::add_edge(*itn0, *itn1, bgraph);
			}
		}
	}
	return bgraph;
}

#endif
