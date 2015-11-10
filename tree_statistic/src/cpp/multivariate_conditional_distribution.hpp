#ifndef DISCRETE_MULTIVARIATE_CONDITIONAL_DISTRIBUTION_HPP
#define DISCRETE_MULTIVARIATE_CONDITIONAL_DISTRIBUTION_HPP

template<typename T> DiscreteMultivariateConditionalReestimation<T>::DiscreteMultivariateConditionalReestimation()
{
	nb_variables = 0;
	nb_covariables = 0;
	nb_elements = 0;
}


template<typename T> DiscreteMultivariateConditionalReestimation<T>::DiscreteMultivariateConditionalReestimation(const DiscreteMultivariateConditionalReestimation<T>& dmcr)
{
	conditional_histogram = dmcr.conditional_histogram;
	nb_variables = dmcr.nb_variables;
	nb_covariables = dmcr.nb_covariables;
	nb_elements = dmcr.nb_elements;
}

template<typename T> DiscreteMultivariateConditionalReestimation<T>::DiscreteMultivariateConditionalReestimation(unsigned int inb_covariables, unsigned int inb_variables)
{
	nb_variables = inb_variables;
	nb_covariables = inb_covariables;
	nb_elements = 0;
}

template<typename T> DiscreteMultivariateConditionalReestimation<T>::DiscreteMultivariateConditionalReestimation(const std::map<std::vector<int>, T>& histogram, const std::set<unsigned int>& covariables, std::set<unsigned int> variables)
{
	nb_covariables = covariables.size();
	nb_variables = variables.size();
	std::vector<int> covariates(covariables.size());
	std::vector<int> variates(variables.size());
	typename std::map< std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::iterator ith;
	typename std::map<std::vector<int>, T>::iterator itmh;
	std::map<std::vector<int>, T> tmap;
	for(typename std::map<std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
	{
		for(std::set<unsigned int>::const_iterator itc = covariables.begin(); itc != covariables.end(); ++itc)
			covariates[distance(covariables.begin(), itc)] = ((*it).first)[*itc];
		for(std::set<unsigned int>::const_iterator itc = variables.begin(); itc != variables.end(); ++itc)
			variates[distance(variables.begin(), itc)] = ((*it).first)[*itc];
		ith = conditional_histogram.find(covariates);
		if(ith == conditional_histogram.end())
		{
			tmap.clear();
			tmap.insert(std::pair<std::vector<int>, T>(variates, (*it).second));
			conditional_histogram.insert(std::pair<std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >(covariates, std::pair<T, std::map<std::vector<int>, T> >((*it).second, tmap)));
		} else {
			itmh = (((*ith).second).second).find(variates);
			if(itmh == (((*ith).second).second).end())
				(((*ith).second).second).insert(std::pair<std::vector<int>, T>(variates, (*it).second));
			else
				(*itmh).second += (*it).second;
			((*ith).second).first += (*it).second;
		}
	}
	nb_elements = 0;
	for(typename std::map<std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
		nb_elements += ((*it).second).first;
}

template<typename T> DiscreteMultivariateConditionalReestimation<T>::~DiscreteMultivariateConditionalReestimation()
{
}

template<typename T> void DiscreteMultivariateConditionalReestimation<T>::update(const std::vector<int>& covariates, const std::map<std::vector<int>, T>& marginal_histogram)
{
	assert(covariates.size() == nb_covariables);
	typename std::map<std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::iterator itcm;
	typename std::map<std::vector<int>, T>::iterator itm;
	itcm = conditional_histogram.find(covariates);
	if(itcm == conditional_histogram.end())
	{
		conditional_histogram.insert(std::pair<std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >(covariates, std::pair<T, std::map<std::vector<int>, T> >(0, marginal_histogram)));
		itcm = conditional_histogram.find(covariates);
		for(typename std::map<std::vector<int>, T>::iterator it = (((*itcm).second).second).begin(); it != (((*itcm).second).second).end(); ++it)
			((*itcm).second).first += (*it).second;
	} else {
		for(typename std::map<std::vector<int>, T>::const_iterator it = marginal_histogram.begin(); it != marginal_histogram.end(); ++it)
		{
			itm = (((*itcm).second).second).find((*it).first);
			if(itm == (((*itcm).second).second).end())
				(((*itcm).second).second).insert(*it);
			else
				(*itm).second += (*it).second;
			((*itcm).second).first += (*it).second;
		}
	}
}

template<typename T> void DiscreteMultivariateConditionalReestimation<T>::nb_value_computation()
{
	nb_elements = 0;
	for(typename std::map<std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
		nb_elements += ((*it).second).first;
}

template<typename T> double DiscreteMultivariateConditionalReestimation<T>::entropy_computation() const
{
	double e = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<std::vector<int>, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			e -= (*itb).second * (log((*itb).second) - log(((*it).second).first));
	}
	return e;
}

template<typename T> std::ostream& DiscreteMultivariateConditionalReestimation<T>::display(std::ostream& os) const
{
	std::vector< std::vector<std::string> > Output(1, std::vector<std::string>(3,""));
	Output[0][0] = "Covariates";
	Output[0][1] = "Responses";
	Output[0][2] = "Occurancy";
	unsigned int p = 1;
	bool not_setted;
	for(typename std::map< std::vector<int>, std::pair<T, std::map< std::vector<int>, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		not_setted = true;
		for(typename std::map<std::vector<int>, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			Output.push_back(std::vector<std::string>(3,""));
			Output[p][1] = toString(toString((*itb).first));
			Output[p][2] = toString((*itb).second);
			if(distance((((*it).second).second).begin(), itb) >= (int)(distance((((*it).second).second).begin(), (((*it).second).second).end())/2.) && not_setted)
			{
				Output[p][0] = toString((*it).first);
				not_setted = false;
			}
			p++;
		}
	}
	for(std::vector< std::vector<std::string> >::iterator it = Output.begin(); it != Output.end(); ++it)
	{
		for(std::vector<std::string>::iterator itb = (*it).begin(); itb != (*it).end(); ++itb)
			os << " " << *itb << " ";
		os << std::endl;
	}
	return os;
}

/*
template<> bool DiscreteMultivariateConditionalReestimation<unsigned int>::check_vectors(const Vectors& vectors, int response)
{
	/*
	 * TODO
	 *
	return true;
}*/

/*
template<> bool DiscreteMultivariateConditionalReestimation<double>::check_vectors(const Vectors& vectors, int response)
{
	/*
	 * TODO
	 *
	return true;
}*/

/*
template<typename T> bool DiscreteMultivariateConditionalReestimation<T>::check_correspondance_link_family(int family, int link) const
{
	switch(family)
	{
		case BINOMIAL :
			if(link == DiscreteMultivariateConditionalParametric::LOGIT || link == DiscreteMultivariateConditionalParametric::LOGLOG)
				return true;
			else
				return false;
			break;
		case POISSON :
			if(link == DiscreteMultivariateConditionalParametric::LOG || link == DiscreteMultivariateConditionalParametric::IDENTITY)
				return true;
			else
				return false;
			break;
	}
}*/

template<typename T> DiscreteMultivariateConditionalDistribution* DiscreteMultivariateConditionalReestimation<T>::distribution_estimation() const 
{
	DiscreteMultivariateConditionalMass cmass;
	std::map< std::vector<int>, double> mass;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		mass.clear();
		for(typename std::map<std::vector<int>, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			mass.insert(std::pair<std::vector<int>, double>((*itb).first, (double)((*itb).second)));
		cmass.insert(std::pair< std::vector<int>, std::map< std::vector<int>, double> >((*it).first, mass));
	}
	return new DiscreteMultivariateConditionalDistribution(nb_covariables, nb_variables, cmass);
}

template<typename T> DiscreteUnivariateConditionalReestimation<T>* DiscreteMultivariateConditionalReestimation<T>::get_sum_reestimation() const
{
	DiscreteUnivariateConditionalReestimation<T> *sum_reestimation = new DiscreteUnivariateConditionalReestimation<T>(nb_covariables);
	std::map<int, T> sum_histogram;
	typename std::map<int, T>::iterator itsh;
	unsigned int sum;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		sum_histogram.clear();
		for(typename std::map<std::vector<int>, T>::const_iterator itm = (((*it).second).second).begin(); itm != (((*it).second).second).end(); ++itm)
		{
			sum = 0;
			for(std::vector<int>::const_iterator itv = ((*itm).first).begin(); itv != ((*itm).first).end(); ++itv)
				sum += *itv;
			itsh = sum_histogram.find(sum);
			if(itsh == sum_histogram.end())
				sum_histogram.insert(std::pair<int, T>(sum, (*itm).second));
			else
				(*itsh).second += (*itm).second;
		}
		sum_reestimation->update((*it).first, sum_histogram);
	}
	sum_reestimation->nb_value_computation();
	return sum_reestimation;
}


template<typename T> DiscreteMultivariateConditionalParametric* DiscreteMultivariateConditionalReestimation<T>::type_parametric_estimation(unsigned int maxits) const
{
	DiscreteUnivariateConditionalReestimation<T> *sum_reestimation = get_sum_reestimation();
	DiscreteMultivariateConditionalParametric *best = NULL;//, *current = NULL;
	//best = new DiscreteMultivariateConditionalParametric(*multinomial_sum_compound_estimation(*(sum_reestimation->type_parametric_estimation(maxits))));
	/*current = new DiscreteMultivariateConditionalParametric(*(parametric_estimation(BINOMIAL, DiscreteMultivariateConditionalParametric::LOGLOG, maxits)));
	if(likelihood_computation(*best) < likelihood_computation(*current))
	{
		delete best;
		best = new DiscreteMultivariateConditionalParametric(*current);
	}*/
	//delete current;
	return multinomial_sum_compound_estimation(*(sum_reestimation->type_parametric_estimation(maxits)));
}


template<typename T> DiscreteMultivariateConditionalParametric* DiscreteMultivariateConditionalReestimation<T>::parametric_estimation(int family, int link, int ident, unsigned int maxits) const
{
	DiscreteUnivariateConditionalReestimation<T> *sum_reestimation = get_sum_reestimation();
	return multinomial_sum_compound_estimation(*(sum_reestimation->parametric_estimation(family, link, maxits)));
}

template<typename T> DiscreteMultivariateConditionalParametric* DiscreteMultivariateConditionalReestimation<T>::multinomial_sum_compound_estimation(const DiscreteUnivariateConditionalParametric& sum_compound) const
{
	std::cout << nb_variables << std::endl;
	std::vector<T> sums = std::vector<T>(nb_variables,0);
	for(typename std::map< std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<std::vector<int>, T>::const_iterator itm = (((*it).second).second).begin(); itm != (((*it).second).second).end(); ++itm)
		{
			for(std::vector<int>::const_iterator itv = ((*itm).first).begin(); itv != ((*itm).first).end(); ++itv)
				sums[distance(((*itm).first).begin(), itv)] += *itv * (*itm).second;
		}
	}
	T sum = 0;
	for(typename std::vector<T>::iterator it = sums.begin(); it != sums.end(); ++it)
		sum += *it;
	Parameters parameters(nb_variables);
	for(typename std::vector<T>::iterator it = sums.begin(); it != sums.end(); ++it)
		parameters[distance(sums.begin(), it)] = *it/((double)sum);
	return new DiscreteMultivariateConditionalParametric(nb_covariables, nb_variables, DiscreteMultivariateConditionalParametric::MULTINOMIALSUMCOMPOUND, sum_compound, parameters);
}

template<typename T> double DiscreteMultivariateConditionalReestimation<T>::likelihood_computation(DiscreteMultivariateConditionalDistribution& ducd, bool log_computation) const
{
	double lh;
	if(log_computation)
		lh = 0;
	else
		lh = 1;
	for(typename std::map< std::vector<int>, std::pair<T, std::map< std::vector<int>, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map< std::vector<int>, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			if(log_computation)
				lh += (*itb).second * ducd.get_mass((*it).first, (*itb).first, true);
			else
				lh *= pow(ducd.get_mass((*it).first, (*itb).first, false), (*itb).second);
		}
	}
	return lh;
}

template<typename T> double DiscreteMultivariateConditionalReestimation<T>::kullback_distance_computation(DiscreteMultivariateConditionalDistribution& dmcd) const
{
	double d = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map< std::vector<int>, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map< std::vector<int>, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			d += (*itb).second/((double)nb_elements) * (log((*itb).second) - log(((*it).second).first) - dmcd.get_mass((*it).first, (*itb).first, true));
	}
	return d;
}

/*
template<typename T> std::multimap<double, std::pair<double, double> > DiscreteMultivariateConditionalReestimation<T>::get_pearson_residuals(DiscreteMultivariateConditionalDistribution& ducd) const
{
	std::multimap<double, std::pair<double, double> > pearson_residuals;
	double m, s;
	for(typename std::map<std::vector<int>, std::pair<unsigned int, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		m = ducd.get_mean((*it).first);
		s = pow(ducd.get_variance((*it).first), 0.5);
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			pearson_residuals.insert(std::pair<double, std::pair<double, double> >(m, std::pair<double, double>(((*itb).first-m)/s, (*itb).second/((double)(((*it).second).first)))));
	}
	return pearson_residuals;
}*/
#endif
