#ifndef DISCRETE_UNIVARIATE_CONDITIONAL_DISTRIBUTION_CPP
#define DISCRETE_UNIVARIATE_CONDITIONAL_DISTRIBUTION_CPP

#include "univariate_conditional_distribution.h"

using namespace stat_tool;

DiscreteUnivariateConditionalDistribution::DiscreteUnivariateConditionalDistribution()
{
	nb_covariables = 0;
	nb_parameters = 0;
}

DiscreteUnivariateConditionalDistribution::DiscreteUnivariateConditionalDistribution(const DiscreteUnivariateConditionalDistribution& ducd)
{
	nb_parameters = ducd.nb_parameters;
	nb_covariables = ducd.nb_covariables;
	mass = ducd.mass;
	means = ducd.means;
	variances = ducd.variances;
}

DiscreteUnivariateConditionalDistribution::DiscreteUnivariateConditionalDistribution(unsigned int inb_covariables, const DiscreteUnivariateConditionalMass& imass)
{
	nb_covariables = inb_covariables;
	mass = imass;
	double mass_sum;
	nb_parameters = 0;
	for(DiscreteUnivariateConditionalMass::iterator it = mass.begin(); it != mass.end(); ++it)
	{
		mass_sum = 0;
		nb_parameters += ((*it).second).size()-1;
		for(std::map<int, double>::iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
			mass_sum += (*itb).second;
		if(mass_sum != 1)
		{
			for(std::map<int, double>::iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
				(*itb).second /= mass_sum;
		}
	}
}

double DiscreteUnivariateConditionalDistribution::get_mass(const std::vector<int>& covariates, int response, bool log_computation)
{
	DiscreteUnivariateConditionalMass::iterator it = mass.find(covariates);
	std::map<int, double> tmap;
	double p;
	if(it == mass.end())
	{
		p = mass_computation(covariates, response);
		tmap.clear();
		tmap.insert(std::pair<int, double>(response, p));
		mass.insert(std::pair< std::vector<int>, std::map<int, double> >(covariates, tmap));
	} else {
		std::map<int, double>::iterator itb = ((*it).second).find(response);
		if(itb == ((*it).second).end())
		{
			p = mass_computation(covariates, response);
			((*it).second).insert(std::pair<int, double>(response, p));
		} else {
			p = (*itb).second;
		}
	}
	if(log_computation)
		return log(p);
	else
		return p;
}

double DiscreteUnivariateConditionalDistribution::get_mean(const std::vector<int>& covariates)
{
	std::map<std::vector<int>, double>::iterator it = means.find(covariates);
	double m;
	if(it == means.end())
	{
		m = mean_computation(covariates);
		means.insert(std::pair< std::vector<int>, double>(covariates, m));
		return m;
	} else {
		return (*it).second;
	}
}

double DiscreteUnivariateConditionalDistribution::get_variance(const std::vector<int>& covariates)
{
	std::map<std::vector<int>, double>::iterator it = variances.find(covariates);
	double v;
	if(it == variances.end())
	{
		v = variance_computation(covariates);
		variances.insert(std::pair< std::vector<int>, double>(covariates, v));
		return v;
	} else {
		return (*it).second;
	}
}

unsigned int DiscreteUnivariateConditionalDistribution::get_nb_parameters() const
{
	return nb_parameters;
}

unsigned int DiscreteUnivariateConditionalDistribution::get_nb_covariables() const
{
	return nb_covariables;
}

int DiscreteUnivariateConditionalDistribution::simulation(const std::vector<int>& covariates) const
{
	DiscreteUnivariateConditionalMass::const_iterator it = mass.find(covariates);
	//std::map<int, double>::const_iterator itb;
	if(it == mass.end())
	{
		return 0;
	} else {
		bool success = false;
		std::map<int, double>::const_iterator itb;
		do
		{
			itb = ((*it).second).begin();
			advance(itb, (int)(generator()*(((*it).second).size()+1)));
		} while ((*itb).second < generator());
		return (*itb).first;
	}
}

DiscreteUnivariateConditionalData* DiscreteUnivariateConditionalDistribution::simulation(const std::vector< std::vector<int> >& covariates) const
{
	std::vector<int> response(covariates.size());
	for(std::vector< std::vector<int> >::const_iterator it = covariates.begin(); it != covariates.end(); ++it)
		response[distance(covariates.begin(), it)] = simulation(*it);
	return new DiscreteUnivariateConditionalData(nb_covariables, covariates, response);
}

double DiscreteUnivariateConditionalDistribution::mass_computation(const std::vector<int>& covariates, int response) const
{
	return 0;
}

double DiscreteUnivariateConditionalDistribution::mean_computation(const std::vector<int>& covariates) const
{
	DiscreteUnivariateConditionalMass::const_iterator it = mass.find(covariates);
	if(it == mass.end())
	{
		return D_NAN;
	} else {
		double m=0;
		for(std::map<int, double>::const_iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
		{
			m += (*itb).first * (*itb).second; 
		}
		return m;
	}
}

double DiscreteUnivariateConditionalDistribution::variance_computation(const std::vector<int>& covariates)
{
	double m = get_mean(covariates);
	if(!((boost::math::isnan)(m)))
	{
		DiscreteUnivariateConditionalMass::const_iterator it = mass.find(covariates);
		if(it == mass.end())
		{
			return D_NAN;
		} else {
			double v=0;
			for(std::map<int, double>::const_iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
			{
				v += pow((*itb).first-m,2) * (*itb).second; 
			}
			return v;
		}
	} else {
		return D_NAN;
	}
}

double DiscreteUnivariateConditionalDistribution::generator() const
{
	return (double)rand()/((double)RAND_MAX + 1);
}

DiscreteUnivariateConditionalParametric::DiscreteUnivariateConditionalParametric() : DiscreteUnivariateConditionalDistribution()
{
}

DiscreteUnivariateConditionalParametric::DiscreteUnivariateConditionalParametric(const DiscreteUnivariateConditionalParametric& ducp) : DiscreteUnivariateConditionalDistribution(ducp)
{
	family = ducp.family;
	link = ducp.link;
	parameters = ducp.parameters;
}

DiscreteUnivariateConditionalParametric::DiscreteUnivariateConditionalParametric(unsigned int inb_covariables, int ifamily, int ilink, const Parameters& iparameters)
{
	nb_covariables = inb_covariables;
	family = ifamily;
	link = ilink;
	parameters = iparameters;
	nb_parameters = parameters.size();
}

int DiscreteUnivariateConditionalParametric::simulation(const std::vector<int>& covariates) const
{
	double e = mean_computation(covariates);
	DiscreteParametric *dp = new DiscreteParametric(family, 0, 1, e, e);
	int res = dp->simulation();
	delete dp;
	return res;
}

double DiscreteUnivariateConditionalParametric::mean_computation(const std::vector<int>& covariates) const
{
	double e = parameters[0];
	for(std::vector<int>::const_iterator it = covariates.begin(); it != covariates.end(); ++it)
		e += parameters[distance(covariates.begin(), it)+1]*(*it);
	switch(link)
	{
		case IDENTITY :
			return MAX(0,e);
			break;
		case LOG :
			return exp(e);
			break;
		case LOGIT :
			return exp(e)/(1+exp(e));
			break;
		case LOGLOG :
			return 1-exp(-exp(e));
			break;
	}
}

double DiscreteUnivariateConditionalParametric::variance_computation(const std::vector<int>& covariates)
{
	double e = parameters[0];
	for(std::vector<int>::const_iterator it = covariates.begin(); it != covariates.end(); ++it)
		e += parameters[distance(covariates.begin(), it)+1]*(*it);
	switch(link)
	{
		case IDENTITY :
			return MAX(0,e);
			break;
		case LOG :
			return exp(e);
			break;
		case LOGIT :
			return (exp(e)/(1+exp(e)))*(1-(exp(e)/(1+exp(e))));
			break;
		case LOGLOG :
			return (1-(1-exp(-exp(e))))*pow(log((1-(1-exp(-exp(e))))),2)/((1-exp(-exp(e))));
			break;
	}
}

double DiscreteUnivariateConditionalParametric::mass_computation(const std::vector<int>& covariates, int response) const
{
	double e = mean_computation(covariates);
	DiscreteParametric *dp = new DiscreteParametric(family, 0, 1, e, e);
	double res = dp->mass[response];
	delete dp;
	return res;
}

Parameters DiscreteUnivariateConditionalParametric::get_parameters() const
{
	return parameters;
}

int DiscreteUnivariateConditionalParametric::get_family() const
{
	return family;
}

int DiscreteUnivariateConditionalParametric::get_link() const
{
	return link;
}

DiscreteUnivariateConditionalData::DiscreteUnivariateConditionalData() : DiscreteUnivariateConditionalReestimation<unsigned int>()
{
}

DiscreteUnivariateConditionalData::DiscreteUnivariateConditionalData(const DiscreteUnivariateConditionalData& ducd) : DiscreteUnivariateConditionalReestimation<unsigned int>(ducd)
{
}

DiscreteUnivariateConditionalData::DiscreteUnivariateConditionalData(unsigned int inb_covariables, const std::vector< std::vector<int> >& covariates, const std::vector<int>& response)
{
	nb_covariables = inb_covariables;
	nb_elements = response.size();
	std::map< std::vector<int>, std::pair<unsigned int, std::map<int, unsigned int> > >::iterator it0;
	std::map<int, unsigned int> tmap;
	std::map<int, unsigned int>::iterator it1;
	for(unsigned int i = 0; i < response.size(); ++i)
	{
		it0 = conditional_histogram.find(covariates[i]);
		if(it0 == conditional_histogram.end())
		{
			tmap.clear();
			tmap.insert(std::pair<int, unsigned int>(response[i], 1));
			conditional_histogram.insert(std::pair< std::vector<int>, std::pair<unsigned int, std::map<int, unsigned int> > >(covariates[i], std::pair<unsigned int, std::map<int, unsigned int> >(1, tmap)));
		} else {
			it1 = (((*it0).second).second).find(response[i]);
			if(it1 == (((*it0).second).second).end())
			{
				(((*it0).second).second).insert(std::pair<int, unsigned int>(response[i], 1));
				(((*it0).second).first) += 1;
			} else {
				(*it1).second += 1;
				(((*it0).second).first) += 1;
			}
		}
	}
	nb_different_elements = 0;
	for(std::map<std::vector<int>, std::pair<unsigned int, std::map<int, unsigned int> > >::iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		nb_different_elements += (((*it).second).second).size();
	}
}


DiscreteUnivariateConditionalData::~DiscreteUnivariateConditionalData()
{
}

/*
bool DiscreteUnivariateConditionalData::check_vectors(const Vectors& vectors, int response)
{
	return DiscreteUnivariateConditionalReestimation<unsigned int>::check_vectors(vectors, response);
}*/

bool DiscreteUnivariateConditionalData::check_correspondance_link_family(int family, int link) const
{
	switch(family)
	{
		case BINOMIAL :
			if(link == DiscreteUnivariateConditionalParametric::LOGIT || link == DiscreteUnivariateConditionalParametric::LOGLOG)
				return true;
			else
				return false;
			break;
		case POISSON :
			if(link == DiscreteUnivariateConditionalParametric::LOG || link == DiscreteUnivariateConditionalParametric::IDENTITY)
				return true;
			else
				return false;
			break;
	}	
}

#endif
