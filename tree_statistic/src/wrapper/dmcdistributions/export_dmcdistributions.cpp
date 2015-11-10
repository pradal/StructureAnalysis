#include <boost/shared_ptr.hpp>
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;

#include "tree_statistic/multivariate_conditional_distribution.h"

#define WRAP DiscreteMultivariateConditionalDistributionWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteMultivariateConditionalDistribution > DMCD_wrapper_init(int nb_covariables, int nb_variables, boost::python::list l_mass)
		{
			std::map< std::vector<int>, std::map< std::vector<int>, double> > mass;
			std::map< std::vector<int>, std::map< std::vector<int>, double> >::iterator it0;
			std::map<std::vector<int>, double> tmass;
			std::map<std::vector<int>, double>::iterator it1;
			std::vector<int> covariates(nb_covariables);
			std::vector<int> variates(nb_variables);
			for(unsigned int i = 0; i < len(l_mass); ++i)
			{
				for(unsigned int j = 0; j < len(l_mass[i][0]); ++j)
					covariates[j] = extract<int>(l_mass[i][0][j]);
				it0 = mass.find(covariates);
				if(it0 != mass.end())
				{
					tmass = (*it0).second;
					mass.erase(it0);
				}	else {
					tmass.clear();
				}
				for(unsigned int j = 0; j < len(l_mass[i][1]); ++j)
				{
					for(unsigned int k = 0; k < len(l_mass[i][1][j]); ++k)
						variates[k] = extract<int>(l_mass[i][1][j][k]);
					it1 = tmass.find(variates);
					if(it1 == tmass.end())
						tmass.insert(std::pair<std::vector<int>, double>(variates, extract<double>(l_mass[i][2][j])));
					else
						(*it1).second += extract<double>(l_mass[i][2][j]);
				}
				mass.insert(std::pair<std::vector<int>, std::map<std::vector<int>, double> >(covariates, tmass));
			}
			DiscreteMultivariateConditionalDistribution *ducd = new DiscreteMultivariateConditionalDistribution(nb_covariables, nb_variables, mass);
			return boost::shared_ptr< DiscreteMultivariateConditionalDistribution >(ducd);
		}
		static float DMCD_wrapper_get_mass(DiscreteMultivariateConditionalDistribution& dmcd, boost::python::list l_covariates, boost::python::list l_responses)
		{
			std::vector<int> covariates(len(l_covariates));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
				covariates[i] = extract<int>(l_covariates[i]);
			std::vector<int> responses(len(l_responses));
			for(unsigned int i = 0; i < len(l_responses); ++i)
				responses[i] = extract<int>(l_responses[i]);
			return dmcd.get_mass(covariates, responses, false);
		}
		static int DMCD_wrapper_get_nb_parameters(const DiscreteMultivariateConditionalDistribution& dmcd)
		{
			return dmcd.get_nb_parameters();
		}
		static int DMCD_wrapper_get_nb_covariables(const DiscreteMultivariateConditionalDistribution& dmcd)
		{
			return dmcd.get_nb_covariables();
		}
		static int DMCD_wrapper_get_nb_variables(const DiscreteMultivariateConditionalDistribution& dmcd)
		{
			return dmcd.get_nb_variables();
		}
		static boost::python::list DMCD_wrapper_get_means(DiscreteMultivariateConditionalDistribution& dmcd, boost::python::list l_covariates)
		{
			std::vector<int> covariates(len(l_covariates));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
				covariates[i] = extract<int>(l_covariates[i]);
			std::vector<double> means = dmcd.get_means(covariates);
			boost::python::list l_means;
			for(std::vector<double>::iterator it = means.begin(); it != means.end(); ++it)
				l_means.append(*it);
			return l_means;
		}
		static boost::python::list DMCD_wrapper_get_variances_covariances(DiscreteMultivariateConditionalDistribution& dmcd, boost::python::list l_covariates)
		{
			std::vector<int> covariates(len(l_covariates));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
				covariates[i] = extract<int>(l_covariates[i]);
			std::vector< std::vector<double> > variances_covariances = dmcd.get_variances_covariances(covariates);
			boost::python::list l_vc, l_l_vc;
			for(std::vector< std::vector<double> >::iterator it = variances_covariances.begin(); it != variances_covariances.end(); ++it)
			{
				l_vc = boost::python::list();
				for(std::vector<double>::iterator itb = (*it).begin(); itb != (*it).end(); ++itb)
					l_vc.append(*itb);
				l_l_vc.append(l_vc);
			}
			return l_l_vc;
		}
		static DiscreteMultivariateConditionalData* DMCD_wrapper_simulation(const DiscreteMultivariateConditionalDistribution& dmcd, boost::python::list l_covariates)
		{
			std::vector< std::vector<int> > covariates(len(l_covariates), std::vector<int>(len(l_covariates[0])));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
			{
				for(unsigned int j = 0; j < len(l_covariates[i]); ++j)
				{
					covariates[i][j] = extract<int>(l_covariates[i][j]);
				}
			}
			return new DiscreteMultivariateConditionalData(*(dmcd.simulation(covariates)));
		}
};

void class_discrete_multivariate_conditional_distribution()
{
	class_< DiscreteMultivariateConditionalDistribution >
		("_DiscreteMultivariateConditionalDistribution", "Discrete Multivariate Conditional Distributions", init< const DiscreteMultivariateConditionalDistribution& >())
		.def("__init__", make_constructor(WRAP::DMCD_wrapper_init))
		.def("_GetMass", WRAP::DMCD_wrapper_get_mass, "_GetMass(self, list, list) -> float"
			"Return mass of list with list as covariates")
		.def("_GetNbParameters", WRAP::DMCD_wrapper_get_nb_parameters, "_GetNbParameters(self) -> int"
			"Return number of parameters of _DiscreteMultivariateConditionalDistribution")
		.def("_GetNbCovariables", WRAP::DMCD_wrapper_get_nb_covariables, "_GetNbCovariables(self) -> int"
			"Return number of covariables of _DiscreteMultivariateConditionalDistribution")
		.def("_GetNbVariables", WRAP::DMCD_wrapper_get_nb_variables, "_GetNbVariables(self) -> int"
			"Return number of variables of _DiscreteMultivariateConditionalDistribution")
		.def("_GetMeans", WRAP::DMCD_wrapper_get_means, "_GetMeans(self, list) -> list"
			"Return means of with list as covariates")
		.def("_GetVariancesCovariances", WRAP::DMCD_wrapper_get_variances_covariances, "_GetVariancesCovariances(self, list) -> list"
			"Return variance-covariance matrix with list as covariates")
		.def("_Simulate", WRAP::DMCD_wrapper_simulation, return_value_policy< manage_new_object >(), "_Simulate(self, list) -> _DiscreteUnivariateConditionalData"
			"Return _DiscreteUnivariateConditionalDistribution simulation with list as covariates")
	;
};

#undef WRAP

template<int num> struct DMCDistributions { int v; enum { value = num };
DMCDistributions(int _v) : v(_v) { } operator int() const { return v; } };

#define WRAP DUCPEnumWrap

void enum_dmcp_distributions_ids()
{
	enum_<DMCDistributions<1> >("DMCDistributions")
		.value("MultinomialSumCompound", DiscreteMultivariateConditionalParametric::MULTINOMIALSUMCOMPOUND)
		.export_values()
	;
};

#undef WRAP

#define WRAP DiscreteMultivariateConditionalParametricWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteMultivariateConditionalParametric > DMCP_wrapper_init(int nb_covariables, int nb_variables, int ident, const DiscreteUnivariateConditionalParametric& ducp, boost::python::list l_parameters)
		{
			Parameters parameters(len(l_parameters));
			for(unsigned int i = 0; i < len(l_parameters); ++i)
				parameters[i] = extract<double>(l_parameters[i]);
			DiscreteMultivariateConditionalParametric *dmcp = new DiscreteMultivariateConditionalParametric(nb_covariables, nb_variables, ident, ducp, parameters);
			return boost::shared_ptr< DiscreteMultivariateConditionalParametric >(dmcp);
		}
		static boost::python::list DMCP_wrapper_get_parameters(const DiscreteMultivariateConditionalParametric& dmcp)
		{
			Parameters parameters = dmcp.get_parameters();
			boost::python::list l_parameters;
			for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it)
				l_parameters.append(*it);
			return l_parameters;
		}
		static DiscreteUnivariateConditionalParametric* DMCP_wrapper_get_compound(const DiscreteMultivariateConditionalParametric& dmcp)
		{
			return dmcp.get_compound();
		}
};

void class_discrete_multivariate_conditional_parametric()
{
	class_< DiscreteMultivariateConditionalParametric, bases< DiscreteMultivariateConditionalDistribution > >
		("_DiscreteMultivariateConditionalParametric", "Discrete Multivariate Conditional Parametric Distributions", init< const DiscreteMultivariateConditionalParametric& >())
		.def("__init__", make_constructor(WRAP::DMCP_wrapper_init))
		.def("_GetParameters", WRAP::DMCP_wrapper_get_parameters, "_GetParameters(self) -> list"
			"Return list of _DiscreteMultivariateConditionalParametric parameters")
		.def("_GetCompound", WRAP::DMCP_wrapper_get_compound, return_value_policy< manage_new_object >(), "_GetCompound(self) -> DiscreteUnivariateConditionalParametric"
			"Return Compounding distribution for _DiscreteMultivariateConditionalParametric")
	;
};

#undef WRAP

#define WRAP DiscreteMultivariateConditionalDataWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteMultivariateConditionalData > DMCDa_wrapper_init_lists(int nb_covariables, int nb_variables, boost::python::list l_covariates,  boost::python::list l_responses)
		{
			std::vector< std::vector<int> > covariates(len(l_covariates), std::vector<int>(nb_covariables));
			std::vector< std::vector<int> > responses(len(l_responses), std::vector<int>(nb_variables));
			std::vector<int> i_covariates(nb_covariables);
			std::vector<int> i_responses(nb_variables);
			for(unsigned int i = 0; i < len(l_responses); ++i)
			{
				for(unsigned int j = 0; j < len(l_covariates[i]); ++j)
					i_covariates[j] = extract<int>(l_covariates[i][j]);
				for(unsigned int j = 0; j < len(l_responses[i]); ++j)
					i_responses[j] = extract<int>(l_responses[i][j]);
				covariates[i] = i_covariates;
				responses[i] = i_responses;
			}
			DiscreteMultivariateConditionalData *dmcr = new DiscreteMultivariateConditionalData(nb_covariables, nb_variables, covariates, responses);
			return boost::shared_ptr< DiscreteMultivariateConditionalData >(dmcr);
		}
		static DiscreteMultivariateConditionalDistribution* DMCDa_wrapper_distribution_estimation(const DiscreteMultivariateConditionalData& dmcr)
		{
			return new DiscreteMultivariateConditionalDistribution(*(dmcr.distribution_estimation()));
		}
		static DiscreteMultivariateConditionalParametric* DMCDa_wrapper_type_parametric_estimation(const DiscreteMultivariateConditionalData& dmcr, int maxits)
		{
			return new DiscreteMultivariateConditionalParametric(*(dmcr.type_parametric_estimation(maxits)));
		}
		static DiscreteMultivariateConditionalParametric* DMCDa_wrapper_parametric_estimation(const DiscreteMultivariateConditionalData& dmcr, int ident, int family, int link, int maxits)
		{
			return new DiscreteMultivariateConditionalParametric(*(dmcr.parametric_estimation(ident, family, link, maxits)));
		}
		static float DMCDa_wrapper_likelihood_computation(const DiscreteMultivariateConditionalData& dmcr, DiscreteMultivariateConditionalDistribution& dmcd, bool log_computation)
		{
			return dmcr.likelihood_computation(dmcd, log_computation);
		}
		static float DMCDa_wrapper_kullback_distance_computation(const DiscreteMultivariateConditionalData& dmcr, DiscreteMultivariateConditionalDistribution& dmcd)
		{
			return dmcr.kullback_distance_computation(dmcd);
		}
		/*static bool DUCDa_wrapper_check_correspondance_link_family(const DiscreteUnivariateConditionalData& ducr, int family, int link)
		{
			return ducr.check_correspondance_link_family(family, link);
		}*/
		static str DUCDa_wrapper_display(const DiscreteMultivariateConditionalData& dmcr)
		{
      std::stringstream s;
      str res;
      dmcr.display(s);
      res = str(s.str());
      return res;
		}/*
		static boost::python::list DUCDa_wrapper_get_pearson_residuals(const DiscreteUnivariateConditionalData& ducr, DiscreteUnivariateConditionalDistribution& ducd)
		{
			std::multimap<double, std::pair<double, double> > c_pearson_residuals = ducr.get_pearson_residuals(ducd);
			boost::python::list x;
			boost::python::list y;
			boost::python::list cex;
			boost::python::list pearson_residuals;
			for(std::multimap<double, std::pair<double, double> >::iterator it = c_pearson_residuals.begin(); it != c_pearson_residuals.end(); ++it)
			{
				x.append((*it).first);
				y.append(((*it).second).first);
				cex.append(100*((*it).second).second);
			}
			pearson_residuals.append(x);
			pearson_residuals.append(y);
			pearson_residuals.append(cex);
			return pearson_residuals;

		}*/
};

void class_discrete_multivariate_conditional_data()
{
	class_< DiscreteMultivariateConditionalData >
		("_DiscreteMultivariateConditionalData", "Discrete Multivariate Conditional Data", init< const DiscreteMultivariateConditionalData& >())
		.def("__init__", make_constructor(WRAP::DMCDa_wrapper_init_lists))
		.def("_DistributionEstimation", WRAP::DMCDa_wrapper_distribution_estimation, return_value_policy< manage_new_object >(), "_DistributionEstimation(self) -> _DiscreteMultivariateConditionalDistribution"
			"Estimation of a _DiscreteMultivariateConditionalDistribution from a _DiscreteMultivariateConditionalData")
		.def("_TypeParametricEstimation", WRAP::DMCDa_wrapper_type_parametric_estimation, return_value_policy< manage_new_object >(), "_TypeParametricEstimation(self, int) -> _DiscreteMultivariateConditionalParametric"
			"Estimation of a _DiscreteMultivariateConditionalParametric from a _DiscreteUnivariateConditionalData")
		.def("_ParametricEstimation", WRAP::DMCDa_wrapper_parametric_estimation, return_value_policy< manage_new_object >(), "_ParametricEstimation(self, int, int, int) -> _DiscreteMultivariateConditionalParametric"
			"Estimation of a _DiscreteMultivariateConditionalParametric with int family and int link function from a _DiscreteUnivariateConditionalData")
		.def("_Likelihood", WRAP::DMCDa_wrapper_likelihood_computation, "_Likelihood(self, _DiscreteMultivariateConditionalDistribution, bool) -> double"
			"Loglikelihood of _DiscreteMultivariateConditionalData with _DiscreteMultivariateConditionalDistribution")
		.def("_KullbackDistance", WRAP::DMCDa_wrapper_kullback_distance_computation, "_KullbackDistance(self, _DiscreteMultivariateConditionalDistribution) -> float"
			"Kullback Distance between _DiscreteMultivariateConditionalData and _DiscreteMultivariateConditionalDistribution")
		.def("_Display", WRAP::DUCDa_wrapper_display, "_Display(self) -> str"
			"Display conditional histogram of _DiscreteUnivariateConditionalData")
	;
};

#undef WRAP

