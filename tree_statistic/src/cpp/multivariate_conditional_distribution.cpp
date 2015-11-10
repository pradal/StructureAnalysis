#ifndef DISCRETE_MULTIVARIATE_CONDITIONAL_DISTRIBUTION_CPP
#define DISCRETE_MULTIVARIATE_CONDITIONAL_DISTRIBUTION_CPP

#include "multivariate_conditional_distribution.h"

using namespace stat_tool;

DiscreteMultivariateConditionalDistribution::DiscreteMultivariateConditionalDistribution()
{
	nb_variables = 0;
	nb_covariables = 0;
	nb_parameters = 0;
}

DiscreteMultivariateConditionalDistribution::DiscreteMultivariateConditionalDistribution(const DiscreteMultivariateConditionalDistribution& dmcd)
{
	nb_variables = dmcd.nb_variables;
	nb_parameters = dmcd.nb_parameters;
	nb_covariables = dmcd.nb_covariables;
	mass = dmcd.mass;
	means = dmcd.means;
	variances_covariances = dmcd.variances_covariances;
}

DiscreteMultivariateConditionalDistribution::DiscreteMultivariateConditionalDistribution(unsigned int inb_covariables, unsigned int inb_variables, const DiscreteMultivariateConditionalMass& imass)
{
	nb_variables = inb_variables;
	nb_covariables = inb_covariables;
	mass = imass;
	double mass_sum;
	nb_parameters = 0;
	for(DiscreteMultivariateConditionalMass::iterator it = mass.begin(); it != mass.end(); ++it)
	{
		mass_sum = 0;
		nb_parameters += ((*it).second).size()-1;
		for(std::map< std::vector<int>, double>::iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
			mass_sum += (*itb).second;
		if(mass_sum != 1)
		{
			for(std::map< std::vector<int>, double>::iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
				(*itb).second /= mass_sum;
		}
	}
}

DiscreteMultivariateConditionalDistribution::~DiscreteMultivariateConditionalDistribution()
{
}

double DiscreteMultivariateConditionalDistribution::get_mass(const std::vector<int>& covariates, const std::vector<int>& responses, bool log_computation)
{
	DiscreteMultivariateConditionalMass::iterator it = mass.find(covariates);
	std::map<std::vector<int>, double> tmap;
	double p;
	if(it == mass.end())
	{
		p = mass_computation(covariates, responses);
		tmap.clear();
		tmap.insert(std::pair< std::vector<int>, double>(responses, p));
		mass.insert(std::pair< std::vector<int>, std::map<std::vector<int>, double> >(covariates, tmap));
	} else {
		std::map< std::vector<int>, double>::iterator itb = ((*it).second).find(responses);
		if(itb == ((*it).second).end())
		{
			p = mass_computation(covariates, responses);
			((*it).second).insert(std::pair<std::vector<int>, double>(responses, p));
		} else {
			p = (*itb).second;
		}
	}
	if(log_computation)
		return log(p);
	else
		return p;
}

std::vector<double> DiscreteMultivariateConditionalDistribution::get_means(const std::vector<int>& covariates)
{
	std::map<std::vector<int>, std::vector<double> >::iterator it = means.find(covariates);
	std::vector<double> m;
	if(it == means.end())
	{
		m = means_computation(covariates);
		means.insert(std::pair< std::vector<int>, std::vector<double> >(covariates, m));
		return m;
	} else {
		return (*it).second;
	}
}

std::vector< std::vector<double> > DiscreteMultivariateConditionalDistribution::get_variances_covariances(const std::vector<int>& covariates)
{
	std::map<std::vector<int>, std::vector< std::vector<double> > >::iterator it = variances_covariances.find(covariates);
	std::vector< std::vector<double> > v;
	if(it == variances_covariances.end())
	{
		v = variances_covariances_computation(covariates);
		variances_covariances.insert(std::pair< std::vector<int>, std::vector< std::vector<double> > >(covariates, v));
		return v;
	} else {
		return (*it).second;
	}
}

unsigned int DiscreteMultivariateConditionalDistribution::get_nb_parameters() const
{
	return nb_parameters;
}

unsigned int DiscreteMultivariateConditionalDistribution::get_nb_covariables() const
{
	return nb_covariables;
}

unsigned int DiscreteMultivariateConditionalDistribution::get_nb_variables() const
{
	return nb_variables;
}

std::vector<int> DiscreteMultivariateConditionalDistribution::simulation(const std::vector<int>& covariates) const
{
	DiscreteMultivariateConditionalMass::const_iterator it = mass.find(covariates);
	//std::map<int, double>::const_iterator itb;
	if(it == mass.end())
	{
		return std::vector<int>(nb_variables,0);
	} else {
		bool success = false;
		std::map<std::vector<int>, double>::const_iterator itb;
		do
		{
			itb = ((*it).second).begin();
			advance(itb, (int)(generator()*(((*it).second).size()+1)));
		} while ((*itb).second < generator());
		return (*itb).first;
	}
}

DiscreteMultivariateConditionalData* DiscreteMultivariateConditionalDistribution::simulation(const std::vector< std::vector<int> >& covariates) const
{
	std::vector< std::vector<int> > responses(covariates.size());
	for(std::vector< std::vector<int> >::const_iterator it = covariates.begin(); it != covariates.end(); ++it)
		responses[distance(covariates.begin(), it)] = simulation(*it);
	return new DiscreteMultivariateConditionalData(nb_covariables, nb_variables, covariates, responses);
}

double DiscreteMultivariateConditionalDistribution::mass_computation(const std::vector<int>& covariates, const std::vector<int>& responses) const
{
	return 0;
}

std::vector<double> DiscreteMultivariateConditionalDistribution::means_computation(const std::vector<int>& covariates) const
{
	DiscreteMultivariateConditionalMass::const_iterator it = mass.find(covariates);
	if(it == mass.end())
	{
		return std::vector<double>(nb_variables, D_NAN);
	} else {
		std::vector<double> m = std::vector<double>(nb_variables,0);
		for(std::map<std::vector<int>, double>::const_iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
		{
			for(unsigned int i = 0; i < nb_variables; ++i)
				m[i] += ((*itb).first)[i] * (*itb).second; 
		}
		return m;
	}
}

std::vector< std::vector<double> > DiscreteMultivariateConditionalDistribution::variances_covariances_computation(const std::vector<int>& covariates)
{
	std::vector<double> m = get_means(covariates);
	if(!((boost::math::isnan)(m[0])))
	{
		DiscreteMultivariateConditionalMass::const_iterator it = mass.find(covariates);
		if(it == mass.end())
		{
			return std::vector< std::vector<double> >(nb_variables, std::vector<double>(nb_variables, D_NAN));
		} else {
			std::vector< std::vector<double> > v = std::vector< std::vector<double> >(nb_variables, std::vector<double>(nb_variables, 0)) ;
			for(std::map<std::vector<int>, double>::const_iterator itb = ((*it).second).begin(); itb != ((*it).second).end(); ++itb)
			{
				for(unsigned int i = 0; i < nb_variables; ++i)
					v[i][i] += pow(((*itb).first)[i]-m[i],2) * (*itb).second; 
				for(unsigned int i = 0; i < nb_variables-1; ++i)
				{
					for(unsigned int j = i+1; j < nb_variables-1; ++j)
					{
						v[i][j] += ((*itb).first)[i]*((*itb).first)[j]*(*itb).second;
					}
				}
			}
			for(unsigned int i = 0; i < nb_variables-1; ++i)
			{
				for(unsigned int j = i+1; j < nb_variables; ++j)
				{
					v[i][j] -= m[i]*m[j];
					v[j][i] = v[i][j];
				}
			}
			return v;
		}
	} else {
		return std::vector< std::vector<double> >(nb_variables, std::vector<double>(nb_variables, D_NAN));
	}
}

double DiscreteMultivariateConditionalDistribution::generator() const
{
	return (double)rand()/((double)RAND_MAX + 1);
}

DiscreteMultivariateConditionalParametric::DiscreteMultivariateConditionalParametric() : DiscreteMultivariateConditionalDistribution()
{
}

DiscreteMultivariateConditionalParametric::DiscreteMultivariateConditionalParametric(const DiscreteMultivariateConditionalParametric& dmcp) : DiscreteMultivariateConditionalDistribution(dmcp)
{
	ident = dmcp.ident;
	compound = new DiscreteUnivariateConditionalParametric(*dmcp.compound);
	parameters = dmcp.parameters;
}

DiscreteMultivariateConditionalParametric::DiscreteMultivariateConditionalParametric(unsigned int inb_covariables, unsigned int inb_variables, int iident, const DiscreteUnivariateConditionalParametric& ducp, const Parameters& iparameters)
{
	nb_covariables = inb_covariables;
	nb_variables = inb_variables;
	ident = iident;
	compound = new DiscreteUnivariateConditionalParametric(ducp);
	parameters = iparameters;
	double sum=0;
	for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it)
		sum += *it;
	if(sum != 1)
	{
		for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it)
			*it /= sum;
	}
	nb_parameters = compound->get_nb_parameters()+parameters.size()-1;
}

DiscreteMultivariateConditionalParametric::~DiscreteMultivariateConditionalParametric()
{
	delete compound;
}

std::vector<int> DiscreteMultivariateConditionalParametric::simulation(const std::vector<int>& covariates) const
{
	if(covariates.size() != nb_covariables)
	{
		return std::vector<int>(nb_variables, 0);
	} else {
		int sum = compound->simulation(covariates);
		std::vector<int> event(nb_variables, 0);
		double sp = 0;
		DiscreteParametric *mdist = NULL;
		for(unsigned int i = 0; i < nb_variables-1; ++i)
		{
			mdist = new DiscreteParametric(BINOMIAL, 0, sum, parameters[i]/(1.-sp), parameters[i]/(1.-sp));
			event[i] = mdist->simulation();
			delete mdist;
			sp += parameters[i];
			sum -= event[i];
		}
		event[nb_variables-1] = sum;
		return event;
	}
}

DiscreteUnivariateConditionalParametric* DiscreteMultivariateConditionalParametric::get_compound() const
{
	return new DiscreteUnivariateConditionalParametric(*compound);
}

Parameters DiscreteMultivariateConditionalParametric::get_parameters() const
{
	return parameters;
}

double DiscreteMultivariateConditionalParametric::mass_computation(const std::vector<int>& covariates, const std::vector<int>& variates) const
{
	switch(ident)
	{
		case MULTINOMIALSUMCOMPOUND :
			return multinomial_sum_compound_mass_computation(covariates, variates);
			break;
		default:
			return 0.;
			break;
	}
}

std::vector<double> DiscreteMultivariateConditionalParametric::means_computation(const std::vector<int>& covariates) const
{
	switch(ident)
	{
		case MULTINOMIALSUMCOMPOUND :
			return multinomial_sum_compound_means_computation(covariates);
			break;
		default :
			return std::vector<double>(nb_variables, D_NAN);
			break;
	}
}

std::vector< std::vector<double> > DiscreteMultivariateConditionalParametric::variances_covariances_computation(const std::vector<int>& covariates) const
{
	switch(ident)
	{
		case MULTINOMIALSUMCOMPOUND :
			return multinomial_sum_compound_variances_covariances_computation(covariates);
			break;
		default :
			return std::vector< std::vector<double> >(nb_variables, std::vector<double>(nb_variables, D_NAN));
			break;
	}
}

double DiscreteMultivariateConditionalParametric::multinomial_sum_compound_mass_computation(const std::vector<int>& covariates, const std::vector<int>& responses) const
{
	if(covariates.size() != nb_covariables || responses.size() != nb_variables)
	{
		return 0.;
	} else {
		int sum = 0;
		for(std::vector<int>::const_iterator it = responses.begin(); it != responses.end(); ++it)
			sum += *it;
		double p = compound->get_mass(covariates, sum, true);
		p += boost::math::lgamma(sum+1);
		for(unsigned int i = 0; i < nb_variables; ++i)
			p += responses[i] * log(parameters[i]) - boost::math::lgamma(responses[i]+1);
		return exp(p);
	}
}

std::vector<double> DiscreteMultivariateConditionalParametric::multinomial_sum_compound_means_computation(const std::vector<int>& covariates) const
{
	if(covariates.size() != nb_covariables)
	{
		return std::vector<double>(nb_variables, D_NAN);
	} else {
		double mean = compound->get_mean(covariates);
		std::vector<double> m(nb_variables);
		for(unsigned int i = 0; i < nb_variables; ++i)
			m[i] = mean*parameters[i];
		return m;
	}
}

std::vector< std::vector<double> > DiscreteMultivariateConditionalParametric::multinomial_sum_compound_variances_covariances_computation(const std::vector<int>& covariates) const
{
	if(covariates.size() != nb_covariables)
	{
		return std::vector< std::vector<double> >(nb_variables, std::vector<double>(nb_variables, D_NAN));
	} else {
		double mean = compound->get_mean(covariates);
		double variance = compound->get_variance(covariates);
		std::vector< std::vector<double> > v(nb_variables, std::vector<double>(nb_variables));
		for(unsigned int i = 0; i < nb_variables; ++i)
			v[i][i] = parameters[i] * (1 - parameters[i]) * mean + pow(parameters[i], 2) * variance;
		for(unsigned int i = 0; i < nb_variables-1; ++i)
		{
			for(unsigned int j = i+1; j < nb_variables; ++j)
			{
				v[i][j] = parameters[i]*parameters[j]*(variance-mean);
				v[j][i] = v[i][j];
			}
		}
		return v;
	}
}

DiscreteMultivariateConditionalData::DiscreteMultivariateConditionalData() : DiscreteMultivariateConditionalReestimation<unsigned int>()
{
}

DiscreteMultivariateConditionalData::DiscreteMultivariateConditionalData(const DiscreteMultivariateConditionalData& dmcd) : DiscreteMultivariateConditionalReestimation<unsigned int>(dmcd)
{
}

DiscreteMultivariateConditionalData::DiscreteMultivariateConditionalData(unsigned int inb_covariables, unsigned int inb_variables, const std::vector< std::vector<int> >& covariates, const std::vector< std::vector<int> >& responses)
{
	nb_covariables = inb_covariables;
	nb_variables = inb_variables;
	nb_elements = responses.size();
	std::map< std::vector<int>, std::pair<unsigned int, std::map<std::vector<int>, unsigned int> > >::iterator it0;
	std::map<std::vector<int>, unsigned int> tmap;
	std::map<std::vector<int>, unsigned int>::iterator it1;
	for(unsigned int i = 0; i < responses.size(); ++i)
	{
		it0 = conditional_histogram.find(covariates[i]);
		if(it0 == conditional_histogram.end())
		{
			tmap.clear();
			tmap.insert(std::pair<std::vector<int>, unsigned int>(responses[i], 1));
			conditional_histogram.insert(std::pair< std::vector<int>, std::pair<unsigned int, std::map<std::vector<int>, unsigned int> > >(covariates[i], std::pair<unsigned int, std::map<std::vector<int>, unsigned int> >(1, tmap)));
		} else {
			it1 = (((*it0).second).second).find(responses[i]);
			if(it1 == (((*it0).second).second).end())
			{
				(((*it0).second).second).insert(std::pair<std::vector<int>, unsigned int>(responses[i], 1));
				(((*it0).second).first) += 1;
			} else {
				(*it1).second += 1;
				(((*it0).second).first) += 1;
			}
		}
	}
}

DiscreteMultivariateConditionalData::~DiscreteMultivariateConditionalData()
{
}
#endif
