#ifndef DISCRETE_UNIVARIATE_CONDITIONAL_DISTRIBUTION_H
#define DISCRETE_UNIVARIATE_CONDITIONAL_DISTRIBUTION_H

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include <map>
#include <set>
#include <eigen3/Eigen/Dense>
#include <boost/math/special_functions/fpclassify.hpp>

#ifndef D_INF
#define D_INF std::numeric_limits<double>::infinity()
#endif

#ifndef D_NAN
#define D_NAN std::numeric_limits<double>::quiet_NaN()
#endif

#include "functions.hpp"

typedef std::vector<double> Parameters;
typedef std::map< std::vector<int>, std::map<int, double> > DiscreteUnivariateConditionalMass;

class DiscreteUnivariateConditionalDistribution;
class DiscreteUnivariateConditionalParametric;

template<typename T>
class DiscreteUnivariateConditionalReestimation
{
	public:
		DiscreteUnivariateConditionalReestimation();
		DiscreteUnivariateConditionalReestimation(unsigned int inb_covariables);
		DiscreteUnivariateConditionalReestimation(const std::map<std::vector<int>, T>& histogram, const std::set<unsigned int>& covariables, unsigned int variable);
		DiscreteUnivariateConditionalReestimation(const DiscreteUnivariateConditionalReestimation<T>& duce);
		~DiscreteUnivariateConditionalReestimation();
		void update(const std::vector<int>& covariates, const std::map<int, T>& marginal_histogram);
		void nb_value_computation();
		DiscreteUnivariateConditionalDistribution* distribution_estimation() const;
		DiscreteUnivariateConditionalParametric* type_parametric_estimation(unsigned int maxits = 10e2) const;
		DiscreteUnivariateConditionalParametric* parametric_estimation(int family, int link, unsigned int maxits = 10e2) const;
		double entropy_computation() const;
		double likelihood_computation(DiscreteUnivariateConditionalDistribution& ducd, bool log_computation=true) const;
		double kullback_distance_computation(DiscreteUnivariateConditionalDistribution& ducd) const;
		std::ostream& display(std::ostream& os) const;
		std::multimap<double, std::pair<double, double> > get_pearson_residuals(DiscreteUnivariateConditionalDistribution& ducd) const;
		//DiscreteUnivariateConditionalParametric* parametric_variables_selection(const DiscreteUnivariateConditionalParametric& ducd) const;
		
	protected:
		unsigned int nb_covariables;
		unsigned int nb_elements;
		unsigned int nb_different_elements;

		std::map< std::vector<int>, std::pair<T, std::map<int, T> > > conditional_histogram;
		Eigen::VectorXd y_init() const;
		bool make_binomial(Eigen::VectorXd& y) const;
		Eigen::VectorXd identity_link_z_init(const Eigen::VectorXd& y) const;
		Eigen::VectorXd log_link_z_init(const Eigen::VectorXd& y) const;
		Eigen::VectorXd logit_link_z_init(const Eigen::VectorXd& y) const;
		Eigen::VectorXd loglog_link_z_init(const Eigen::VectorXd& y) const;
		void identity_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		void log_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		void logit_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		void loglog_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		Eigen::MatrixXd X_init() const;
		Eigen::MatrixXd W_init(const Eigen::VectorXd& y) const;
		void poisson_identity_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		void poisson_log_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		void binomial_logit_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		void binomial_loglog_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const;
		Eigen::VectorXd beta_compute(const Eigen::MatrixXd& X, const Eigen::MatrixXd& W, const Eigen::VectorXd& z) const;
		Parameters to_stdvector(const Eigen::VectorXd& beta) const;
};

class DiscreteUnivariateConditionalData : public DiscreteUnivariateConditionalReestimation<unsigned int>
{
	public:
		DiscreteUnivariateConditionalData();
		DiscreteUnivariateConditionalData(const DiscreteUnivariateConditionalData& ducd);
		DiscreteUnivariateConditionalData(unsigned int inb_covariables, const std::vector< std::vector<int> >& covariates, const std::vector<int>& response);
		~DiscreteUnivariateConditionalData();
		//bool check_vectors(const Vectors& vectors, int response);
		bool check_correspondance_link_family(int family, int link) const;
};

class DiscreteUnivariateConditionalDistribution
{
	public:
		DiscreteUnivariateConditionalDistribution();
		DiscreteUnivariateConditionalDistribution(const DiscreteUnivariateConditionalDistribution& ducd);
		DiscreteUnivariateConditionalDistribution(unsigned int inb_covariables, const DiscreteUnivariateConditionalMass& imass);
		double get_mass(const std::vector<int>& covariates, int response, bool log_computation=false);
		double get_mean(const std::vector<int>& covariates);
		double get_variance(const std::vector<int>& covariates);
		unsigned int get_nb_parameters() const;
		unsigned int get_nb_covariables() const;
		virtual int simulation(const std::vector<int>& covariates) const;
		DiscreteUnivariateConditionalData* simulation(const std::vector< std::vector<int> >& covariates) const;

	protected:
		unsigned int nb_covariables;
		unsigned int nb_parameters;
		DiscreteUnivariateConditionalMass mass;
		std::map< std::vector<int>, double > means;
		std::map< std::vector<int>, double > variances;

		virtual double mass_computation(const std::vector<int>& covariates, int response) const;
		virtual double mean_computation(const std::vector<int>& covariates) const;
		virtual double variance_computation(const std::vector<int>& covariates);
		double generator() const;
};

class DiscreteUnivariateConditionalParametric : public DiscreteUnivariateConditionalDistribution
{
	public:
		enum
		{
			IDENTITY,
			LOG,
			LOGIT,
			LOGLOG
		};
		DiscreteUnivariateConditionalParametric();
		DiscreteUnivariateConditionalParametric(const DiscreteUnivariateConditionalParametric& ducp);
		DiscreteUnivariateConditionalParametric(unsigned int inb_covariables, int ifamily, int ilink, const Parameters& iparameters);
		virtual int simulation(const std::vector<int>& covariates) const;
		Parameters get_parameters() const;
		int get_family() const;
		int get_link() const;

	protected:
		int family;
		int link;
		Parameters parameters;
		virtual double mass_computation(const std::vector<int>& covariates, int response) const;
		virtual double mean_computation(const std::vector<int>& covariates) const;
		virtual double variance_computation(const std::vector<int>& covariates);
};

#include "univariate_conditional_distribution.hpp"
#endif
