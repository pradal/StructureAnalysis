#ifndef DISCRETE_MULTIVARIATE_CONDITIONAL_DISTRIBUTION_H
#define DISCRETE_MULTIVARIATE_CONDITIONAL_DISTRIBUTION_H

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include <map>
#include <set>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/gamma.hpp>

#ifndef D_INF
#define D_INF std::numeric_limits<double>::infinity()
#endif

#ifndef D_NAN
#define D_NAN std::numeric_limits<double>::quiet_NaN()
#endif

#include "functions.hpp"
#include  "univariate_conditional_distribution.h"

typedef std::vector<double> Parameters;
typedef std::map< std::vector<int>, std::map< std::vector<int>, double> > DiscreteMultivariateConditionalMass;
class DiscreteMultivariateConditionalDistribution;
class DiscreteMultivariateConditionalParametric;

template<typename T>
class DiscreteMultivariateConditionalReestimation
{
	public:
		DiscreteMultivariateConditionalReestimation();
		DiscreteMultivariateConditionalReestimation(const DiscreteMultivariateConditionalReestimation<T>& duce);
		DiscreteMultivariateConditionalReestimation(unsigned int inb_covariables, unsigned int inb_variables);
		DiscreteMultivariateConditionalReestimation(const std::map<std::vector<int>, T>& histogram, const std::set<unsigned int>& covariables, std::set<unsigned int> variables);
		void update(const std::vector<int>& covariates, const std::map<std::vector<int>, T>& marginal_histogram);
		void nb_value_computation();
		~DiscreteMultivariateConditionalReestimation();
		DiscreteMultivariateConditionalDistribution* distribution_estimation() const;
		DiscreteUnivariateConditionalReestimation<T>* get_sum_reestimation() const;
		DiscreteMultivariateConditionalParametric* type_parametric_estimation(unsigned int maxits = 10e2) const;
		DiscreteMultivariateConditionalParametric* parametric_estimation(int ident, int family, int link, unsigned int maxits = 10e2) const;
		double entropy_computation() const;
		double likelihood_computation(DiscreteMultivariateConditionalDistribution& dmcd, bool log_computation=true) const;
		double kullback_distance_computation(DiscreteMultivariateConditionalDistribution& dmcd) const;
		std::ostream& display(std::ostream& os) const;
		
	protected:
		unsigned int nb_variables;
		unsigned int nb_covariables;
		T nb_elements;
		std::map< std::vector<int>, std::pair<T, std::map<std::vector<int>, T> > > conditional_histogram;

		DiscreteMultivariateConditionalParametric* multinomial_sum_compound_estimation(const DiscreteUnivariateConditionalParametric& sum_compound) const;
};

class DiscreteMultivariateConditionalData : public DiscreteMultivariateConditionalReestimation<unsigned int>
{
	public:
		DiscreteMultivariateConditionalData();
		DiscreteMultivariateConditionalData(const DiscreteMultivariateConditionalData& dmcd);
		DiscreteMultivariateConditionalData(unsigned int inb_covariables, unsigned int inb_variables, const std::vector< std::vector<int> >& covariates, const std::vector< std::vector<int> >& responses);
		~DiscreteMultivariateConditionalData();
};

class DiscreteMultivariateConditionalDistribution
{
	public:
		DiscreteMultivariateConditionalDistribution();
		DiscreteMultivariateConditionalDistribution(const DiscreteMultivariateConditionalDistribution& dmcd);
		DiscreteMultivariateConditionalDistribution(unsigned int inb_covariables, unsigned int inb_variables, const DiscreteMultivariateConditionalMass& imass);
		~DiscreteMultivariateConditionalDistribution();
		double get_mass(const std::vector<int>& covariates, const std::vector<int>& responses, bool log_computation=false);
		std::vector<double> get_means(const std::vector<int>& covariates);
		std::vector< std::vector<double> > get_variances_covariances(const std::vector<int>& covariates);
		unsigned int get_nb_parameters() const;
		unsigned int get_nb_covariables() const;
		unsigned int get_nb_variables() const;
		virtual std::vector<int> simulation(const std::vector<int>& covariates) const;
		DiscreteMultivariateConditionalData* simulation(const std::vector< std::vector<int> >& covariates) const;

	protected:
		unsigned int nb_variables;
		unsigned int nb_covariables;
		unsigned int nb_parameters;
		DiscreteMultivariateConditionalMass mass;
		std::map< std::vector<int>, std::vector<double> > means;
		std::map< std::vector<int>, std::vector< std::vector<double> > > variances_covariances;

		virtual double mass_computation(const std::vector<int>& covariates, const std::vector<int>& responses) const;
		virtual std::vector<double> means_computation(const std::vector<int>& covariates) const;
		virtual std::vector< std::vector<double> > variances_covariances_computation(const std::vector<int>& covariates);
		double generator() const;
};

class DiscreteMultivariateConditionalParametric : public DiscreteMultivariateConditionalDistribution
{
	public:
		enum {
			MULTINOMIALSUMCOMPOUND
		};
		DiscreteMultivariateConditionalParametric();
		DiscreteMultivariateConditionalParametric(const DiscreteMultivariateConditionalParametric& dmcp);
		DiscreteMultivariateConditionalParametric(unsigned int inb_covariables, unsigned int inb_variables, int iident, const DiscreteUnivariateConditionalParametric& ducp, const Parameters& iparameters);
		~DiscreteMultivariateConditionalParametric();
		virtual std::vector<int> simulation(const std::vector<int>& covariates) const;
		DiscreteUnivariateConditionalParametric* get_compound() const;
		Parameters get_parameters() const;

	protected:
		int ident;
		DiscreteUnivariateConditionalParametric* compound;
		Parameters parameters;
		virtual double mass_computation(const std::vector<int>& covariates, const std::vector<int>& response) const;
		virtual std::vector<double> means_computation(const std::vector<int>& covariates) const;
		virtual std::vector< std::vector<double> > variances_covariances_computation(const std::vector<int>& covariates) const;

		double multinomial_sum_compound_mass_computation(const std::vector<int>& covariates, const std::vector<int>& responses) const;
		std::vector<double> multinomial_sum_compound_means_computation(const std::vector<int>& covariates) const;
		std::vector< std::vector<double> > multinomial_sum_compound_variances_covariances_computation(const std::vector<int>& covariates) const;
};


#include "multivariate_conditional_distribution.hpp"
#endif
