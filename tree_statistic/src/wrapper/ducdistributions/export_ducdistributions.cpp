#include <boost/shared_ptr.hpp>
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;

#include "tree_statistic/univariate_conditional_distribution.h"

#define WRAP DiscreteUnivariateConditionalDistributionWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteUnivariateConditionalDistribution > DUCD_wrapper_init(int nb_covariables, boost::python::list l_mass)
		{
			std::map< std::vector<int>, std::map<int, double> > mass;
			std::map< std::vector<int>, std::map<int, double> >::iterator it0;
			std::map<int, double> tmass;
			std::map<int, double>::iterator it1;
			std::vector<int> covariates(nb_covariables);
			for(unsigned int i = 0; i < len(l_mass); ++i)
			{
				for(unsigned int j = 0; j < len(l_mass[i][0]); ++j)
					covariates[j] = extract<int>(l_mass[i][0][j]);
				it0 = mass.find(covariates);
				if(it0 != mass.end())
				{
					tmass = (*it0).second;
					mass.erase(it0);
				} else {
					tmass.clear();
				}
				for(unsigned int j = 0; j < len(l_mass[i][1]); ++j)
				{
					it1 = tmass.find(extract<int>(l_mass[i][1][j]));
					if(it1 == tmass.end())
						tmass.insert(std::pair<unsigned int, double>(extract<int>(l_mass[i][1][j]), extract<double>(l_mass[i][2][j])));
					else
						(*it1).second += extract<double>(l_mass[i][2][j]);
				}
				mass.insert(std::pair<std::vector<int>, std::map<int, double> >(covariates, tmass));
			}
			DiscreteUnivariateConditionalDistribution *ducd = new DiscreteUnivariateConditionalDistribution(nb_covariables, mass);
			return boost::shared_ptr< DiscreteUnivariateConditionalDistribution >(ducd);
		}
		static float DUCD_wrapper_get_mass(DiscreteUnivariateConditionalDistribution& ducd, boost::python::list l_covariates, int response)
		{
			std::vector<int> covariates(len(l_covariates));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
				covariates[i] = extract<int>(l_covariates[i]);
			return ducd.get_mass(covariates, response, false);
		}
		static int DUCD_wrapper_get_nb_parameters(const DiscreteUnivariateConditionalDistribution& ducd)
		{
			return ducd.get_nb_parameters();
		}
		static int DUCD_wrapper_get_nb_covariables(const DiscreteUnivariateConditionalDistribution& ducd)
		{
			return ducd.get_nb_covariables();
		}
		static float DUCD_wrapper_get_mean(DiscreteUnivariateConditionalDistribution& ducd, boost::python::list l_covariates)
		{
			std::vector<int> covariates(len(l_covariates));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
				covariates[i] = extract<int>(l_covariates[i]);
			return ducd.get_mean(covariates);
		}
		static float DUCD_wrapper_get_variance(DiscreteUnivariateConditionalDistribution& ducd, boost::python::list l_covariates)
		{
			std::vector<int> covariates(len(l_covariates));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
				covariates[i] = extract<int>(l_covariates[i]);
			return ducd.get_variance(covariates);
		}
		static DiscreteUnivariateConditionalData* DUCD_wrapper_simulation(const DiscreteUnivariateConditionalDistribution& ducd, boost::python::list l_covariates)
		{
			std::vector< std::vector<int> > covariates(len(l_covariates), std::vector<int>(len(l_covariates[0])));
			for(unsigned int i = 0; i < len(l_covariates); ++i)
			{
				for(unsigned int j = 0; j < len(l_covariates[i]); ++j)
				{
					covariates[i][j] = extract<int>(l_covariates[i][j]);
				}
			}
			return ducd.simulation(covariates);
		}
};

void class_discrete_univariate_conditional_distribution()
{
	class_< DiscreteUnivariateConditionalDistribution >
		("_DiscreteUnivariateConditionalDistribution", "Discrete Univariate Conditional Distributions", init< const DiscreteUnivariateConditionalDistribution& >())
		.def("__init__", make_constructor(WRAP::DUCD_wrapper_init))
		.def("_GetMass", WRAP::DUCD_wrapper_get_mass, "_GetMass(self, list, int) -> float"
			"Return mass of int with list as covariates")
		.def("_GetNbParameters", WRAP::DUCD_wrapper_get_nb_parameters, "_GetNbParameters(self) -> int"
			"Return number of parameters of _DiscreteUnivariateConditionalDistribution")
		.def("_GetNbCovariables", WRAP::DUCD_wrapper_get_nb_covariables, "_GetNbCovariates(self) -> int"
			"Return number of covariables of _DiscreteUnivariateConditionalDistribution")
		.def("_GetMean", WRAP::DUCD_wrapper_get_mean, "_GetMean(self, list) -> float"
			"Return mean of with list as covariates")
		.def("_GetVariance", WRAP::DUCD_wrapper_get_variance, "_GetVariance(self, list) -> float"
			"Return variance of with list as covariates")
		.def("_Simulate", WRAP::DUCD_wrapper_simulation, return_value_policy< manage_new_object >(), "_Simulate(self, list) -> _DiscreteUnivariateConditionalReestimation"
			"Return _DiscreteUnivariateConditionalDistribution simulation with list as covariates")
	;
};

#undef WRAP

template<int num> struct DUCPLinks { int v; enum { value = num };
DUCPLinks(int _v) : v(_v) { } operator int() const { return v; } };

#define WRAP DUCPEnumWrap

void enum_ducp_link_ids()
{
	enum_<DUCPLinks<4> >("DUCPLinks")
		.value("Identity", DiscreteUnivariateConditionalParametric::IDENTITY)
		.value("Log", DiscreteUnivariateConditionalParametric::LOG)
		.value("Logit", DiscreteUnivariateConditionalParametric::LOGIT)
		.value("LogLog", DiscreteUnivariateConditionalParametric::LOGLOG)
		.export_values()
	;
};
#undef WRAP

#define WRAP DiscreteUnivariateConditionalParametricWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteUnivariateConditionalParametric > DUCP_wrapper_init(int nb_covariables, int family, int link, boost::python::list l_parameters)
		{
			Parameters parameters(len(l_parameters));
			for(unsigned int i = 0; i < len(l_parameters); ++i)
				parameters[i] = extract<double>(l_parameters[i]);
			DiscreteUnivariateConditionalParametric *ducp = new DiscreteUnivariateConditionalParametric(nb_covariables, family, link, parameters);
			return boost::shared_ptr< DiscreteUnivariateConditionalParametric >(ducp);
		}
		static boost::python::list DUCP_wrapper_get_parameters(const DiscreteUnivariateConditionalParametric& ducp)
		{
			Parameters parameters = ducp.get_parameters();
			boost::python::list l_parameters;
			for(Parameters::iterator it = parameters.begin(); it != parameters.end(); ++it)
				l_parameters.append(*it);
			return l_parameters;
		}
		static int DUCP_wrapper_get_family(const DiscreteUnivariateConditionalParametric& ducp)
		{
			return ducp.get_family();
		}
		static int DUCP_wrapper_get_link(const DiscreteUnivariateConditionalParametric& ducp)
		{
			return ducp.get_link();
		}

};

void class_discrete_univariate_conditional_parametric()
{
	class_< DiscreteUnivariateConditionalParametric, bases< DiscreteUnivariateConditionalDistribution > >
		("_DiscreteUnivariateConditionalParametric", "Discrete Univariate Conditional Parametric Distributions", init< const DiscreteUnivariateConditionalParametric& >())
		.def("__init__", make_constructor(WRAP::DUCP_wrapper_init))
		.def("_GetParameters", WRAP::DUCP_wrapper_get_parameters, "_GetParameters(self) -> list"
			"Return list of _DiscreteUnivariateConditionalParametric parameters")
		.def("_GetFamily", WRAP::DUCP_wrapper_get_family, "_GetFamily(self) -> int"
			"Return corresponding value of family for _DiscreteUnivariateConditionalParametric")
		.def("_GetLink", WRAP::DUCP_wrapper_get_link, "_GetLink(self) -> int"
			"Return corresponding value of link for _DiscreteUnivariateConditionalParametric")
	;
};

#undef WRAP

#define WRAP DiscreteUnivariateConditionalDataWrap

class WRAP
{
	public:
		static boost::shared_ptr< DiscreteUnivariateConditionalData > DUCDa_wrapper_init_lists(int nb_covariables, boost::python::list l_covariates,  boost::python::list l_response)
		{
			std::vector< std::vector<int> > covariates(len(l_covariates), std::vector<int>(len(l_covariates[0])));
			std::vector<int> response(len(l_response));
			std::vector<int> i_covariates(len(l_covariates[0]));
			for(unsigned int i = 0; i < len(l_response); ++i)
			{
				response[i] = extract<int>(l_response[i]);
				for(unsigned int j = 0; j < len(l_covariates[i]); ++j)
					i_covariates[j] = extract<int>(l_covariates[i][j]);
				covariates[i] = i_covariates;
			}
			DiscreteUnivariateConditionalData *ducr = new DiscreteUnivariateConditionalData(nb_covariables, covariates, response);
			return boost::shared_ptr< DiscreteUnivariateConditionalData >(ducr);
		}
		static DiscreteUnivariateConditionalDistribution* DUCDa_wrapper_distribution_estimation(const DiscreteUnivariateConditionalData& ducr)
		{
			return new DiscreteUnivariateConditionalDistribution(*(ducr.distribution_estimation()));
		}
		static DiscreteUnivariateConditionalParametric* DUCDa_wrapper_type_parametric_estimation(const DiscreteUnivariateConditionalData& ducr, int maxits)
		{
			return new DiscreteUnivariateConditionalParametric(*(ducr.type_parametric_estimation(maxits)));
		}
		static DiscreteUnivariateConditionalParametric* DUCDa_wrapper_parametric_estimation(const DiscreteUnivariateConditionalData& ducr, int family, int link, int maxits)
		{
			return new DiscreteUnivariateConditionalParametric(*(ducr.parametric_estimation(family, link, maxits)));
		}
		static float DUCDa_wrapper_likelihood_computation(const DiscreteUnivariateConditionalData& ducr, DiscreteUnivariateConditionalDistribution& ducd, bool log_computation)
		{
			return ducr.likelihood_computation(ducd, log_computation);
		}
		static float DUCDa_wrapper_kullback_distance_computation(const DiscreteUnivariateConditionalData& ducr, DiscreteUnivariateConditionalDistribution& ducd)
		{
			return ducr.kullback_distance_computation(ducd);
		}
		static bool DUCDa_wrapper_check_correspondance_link_family(const DiscreteUnivariateConditionalData& ducr, int family, int link)
		{
			return ducr.check_correspondance_link_family(family, link);
		}
		static str DUCDa_wrapper_display(const DiscreteUnivariateConditionalData& ducr)
		{
      std::stringstream s;
      str res;
      ducr.display(s);
      res = str(s.str());
      return res;
		}
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
		}
};

void class_discrete_univariate_conditional_data()
{
	class_< DiscreteUnivariateConditionalData >
		("_DiscreteUnivariateConditionalData", "Discrete Univariate Conditional Data", init< const DiscreteUnivariateConditionalData& >())
		.def("__init__", make_constructor(WRAP::DUCDa_wrapper_init_lists))
		.def("_DistributionEstimation", WRAP::DUCDa_wrapper_distribution_estimation, return_value_policy< manage_new_object >(), "_DistributionEstimation(self) -> _DiscreteUnivariateConditionalDistribution"
			"Estimation of a _DiscreteUnivariateConditionalDistribution from a _DiscreteUnivariateConditionalData")
		.def("_TypeParametricEstimation", WRAP::DUCDa_wrapper_type_parametric_estimation, return_value_policy< manage_new_object >(), "_TypeParametricEstimation(self, int) -> _DiscreteUnivariateConditionalParametric"
			"Estimation of a _DiscreteUnivariateConditionalParametric from a _DiscreteUnivariateConditionalData")
		.def("_ParametricEstimation", WRAP::DUCDa_wrapper_parametric_estimation, return_value_policy< manage_new_object >(), "_ParametricEstimation(self, int, int, int) -> _DiscreteUnivariateConditionalParametric"
			"Estimation of a _DiscreteUnivariateConditionalParametric with int family and int link function from a _DiscreteUnivariateConditionalData")
		.def("_Likelihood", WRAP::DUCDa_wrapper_likelihood_computation, "_Likelihood(self, _DiscreteUnivariateConditionalDistribution, bool) -> float"
			"Loglikelihood of _DiscreteUnivariateConditionalData with _DiscreteUnivariateConditionalDistribution")
		.def("_KullbackDistance", WRAP::DUCDa_wrapper_kullback_distance_computation, "_KullbackDistance(self, _DiscreteUnivariateConditionalDistribution) -> float"
			"Compute Kullback Distance between _DiscreteUnivariateConditionalData and _DiscreteUnivariateConditionalDistribution")
		.def("_CheckCorrespondanceLinkFamily", WRAP::DUCDa_wrapper_check_correspondance_link_family, "_CheckCorrespondanceLinkFamily(self, family, link)"
			"Check if link is valid with choose family")
		.def("_Display", WRAP::DUCDa_wrapper_display, "_Display(self) -> str"
			"Display conditional histogram of _DiscreteUnivariateConditionalData")
		.def("_GetPearsonResiduals", WRAP::DUCDa_wrapper_get_pearson_residuals, "_GetPearsonResiduals(self, _DiscreteUnivariateConditionalDistribution) -> list"
			"Return Pearson residuals when fitting _DiscreteUnivariateConditionalData with _DiscreteUnivariateConditionalDistribution")
	;
};

#undef WRAP

