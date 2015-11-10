/*------------------------------------------------------------------------------
 *
 *        VPlants.Tree-Statistic : VPlants Tree-Statistic module
 *        HiddenMarkovTrees
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_hmt.cpp 9099 2010-06-08 09:03:00Z pradal $
 *
 *-----------------------------------------------------------------------------*/
// Includes ====================================================================
#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"

// required by "tree_statistic/multivariate_distribution.h"
#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_label.h"
#include "tree_statistic/tree_labels.h"
#include "statiskit/core/data/marginal/multivariate.h"

#include "tree_statistic/multivariate_distribution.h"

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::len
#include <boost/python/make_constructor.hpp>
// definition of boost::python::make_constructor

#include "../errors.h"

#include "stat_tool/wrapper_util.h"
#include "stat_tool/boost_python_aliases.h"

//GLM
#include <boost/math/special_functions/gamma.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <algorithm>

// #include <eigen2/Eigen/Dense>
#include <eigen3/Eigen/Dense>
#include "tree_statistic/core.h"
// Using =======================================================================
using namespace boost::python;
using namespace stat_tool;
using namespace Stat_trees;

// Declarations ================================================================
template<int num> struct UniqueInt { int v; enum { value=num };
UniqueInt(int _v) : v(_v) { } operator int() const { return v; } };

#define WRAP DmdEnumWrap
void enum_dmd_ids()
{
    enum_<UniqueInt<6> >("PMDistributionIds")
        .value("DISCRETE_MULTIVARIATE", Stat_trees::DISCRETE_MULTIVARIATE)
        .value("MIID", Stat_trees::MIID)
        .value("MPOISSON", Stat_trees::MPOISSON)
        .value("MNEGATIVE_MULTINOMIAL", Stat_trees::MNEGATIVE_MULTINOMIAL)
        .value("MMULTINOMIAL", Stat_trees::MMULTINOMIAL)
        .value("MCOMPOUND_MULTINOMIAL", Stat_trees::MCOMPOUND_MULTINOMIAL)
        .export_values()
    ;
};
#undef WRAP

#define WRAP StatTreeErrorWrap

void class_errors()
{
    // Error initialisation
    object stat_tree_errors = import("openalea.tree_statistic._errors");
    // Import StatError
    tree_statistic::StatTreeError = stat_tree_errors.attr("StatTreeError");
};
#undef WRAP

#define WRAP DiscreteMultivariateWrap
class WRAP
{

public :

   static str Dmd_wrapper_display(const DiscreteMultivariateDistribution& dist,
                                  bool exhaustive)
   {
      std::stringstream s;
      str res;

      dist.ascii_write(s, exhaustive);
      res = str(s.str());
      return res;
   }

   static boost::shared_ptr<DiscreteMultivariateDistribution>
   Dmd_wrapper_ascii_read(const char * path, double cumul_threshold)
   {
      StatError error;
      DiscreteMultivariateDistribution *dist = NULL;
      ostringstream error_message;

      dist = discrete_multivariate_ascii_read(error, path,
                                              cumul_threshold);
      if (dist == NULL)
      {
         tree_statistic::throw_stat_tree_error(error);
         // error_message << error;
         // tree_statistic::tree_statistic::throw_stat_tree_error(error_message);
      }

      return boost::shared_ptr<DiscreteMultivariateDistribution>(dist);
   }

   static str Dmd_wrapper_ascii_write0(const DiscreteMultivariateDistribution& dist)
   {
      std::stringstream s;
      str res;

      dist.header_ascii_print(s);
      res = str(s.str());
      return res;
   }


/*   static bool Dmd_wrapper_file_ascii_write2(const DiscreteMultivariateDistribution& dist,
                                               const char * path, bool exhaustive)
   {
      StatError error;
      bool res;
      ostringstream error_message;

      res = dist.ascii_write(error, path, exhaustive);
      if (!res)
      {
         error_message << error;
         tree_statistic::throw_stat_tree_error(error_message);
      }
      return res;
   }

   static bool Dmd_wrapper_file_ascii_write1(const DiscreteMultivariateDistribution& dist,
                                             const char * path)
   { return Dmd_wrapper_file_ascii_write2(dist, path, false); } */

   static int Dmd_wrapper_get_nb_variable(const DiscreteMultivariateDistribution& dist)
   { return dist.nb_variable; }



   static boost::python::list Dmd_wrapper_simulate0(const DiscreteMultivariateDistribution& dist)
   {
      boost::python::list res;
      DiscreteMultivariateDistribution::index_array sim;
      int i;

      sim = dist.simulation();
      for(i = 0; i < dist.nb_variable; i++)
         res.append(sim[i]);
      return res;
   }

     static int Dmd_wrapper_ordinal_association(const Vectors& vec)
     {
         unsigned int i, j, s,f;
         DiscreteMultivariateReestimation<int>* histo = NULL;
         histo = new DiscreteMultivariateReestimation<int>(vec);
         std::vector< std::vector<unsigned int> > vectors = histo->get_elements();
         s = vectors.size();
         for(i = 0; i < s; i++){
             f = histo->get_frequency(vectors[i]);
             for(j = 1; j < f; j++){
                 vectors.push_back(vectors[i]);
             }
         }
         ContingencyTable<unsigned int> *CT = new ContingencyTable<unsigned int>(vectors);
         CT->Parsimoniest();
         return 0;
     }

   static Vectors* Dmd_wrapper_simulate1(const DiscreteMultivariateDistribution& dist,
                                         unsigned int nb_value)

   {
      ostringstream error_message;

      if (nb_value < 1)
      {
         error_message << "Not a DiscreteMultivariateParametric";
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      return dist.simulation(nb_value);
   }

   // to be used with return_value_policy< manage_new_object >()
   // use DEF_RETURN_VALUE otherwise ?
   static DiscreteMultivariateDistribution*
   Dmd_wrapper_ascii_read1(char *filename)
   {
      StatError error;
      ostringstream error_message;



      DiscreteMultivariateDistribution *dist = NULL;

      dist = discrete_multivariate_ascii_read(error, filename);

      if (dist == NULL)
      {
         // tree_statistic::throw_stat_tree_error(error);
         error_message << error;
         tree_statistic::throw_stat_tree_error(error_message);
      }

      return dist;
   }

}; // WRAP


void class_dmdistributions()
{
  class_< DiscreteMultivariateDistribution >
    ("_DiscreteMultivariateDistribution", "Discrete Multivariate Distributions", init< const DiscreteMultivariateDistribution& >())
        .def("__init__", make_constructor(WRAP::Dmd_wrapper_ascii_read))
        .def("__init__", make_constructor(WRAP::Dmd_wrapper_ascii_read1))

        .def("Display", WRAP::Dmd_wrapper_display,
                        "Display(self, bool) -> str \n\n"
                        "Print _DiscreteMultivariateDistribution definition")
        // .def("FileAsciiWrite", WRAP::Dmd_wrapper_file_ascii_write1)
        // .def("FileAsciiWrite", WRAP::Dmd_wrapper_file_ascii_write2)

        .def("NbVariable", WRAP::Dmd_wrapper_get_nb_variable,
                          "NbVariable(self) -> int \n\n"
                          "return the number of variables\n")

        .def("Simulate", WRAP::Dmd_wrapper_simulate0,
                         "Simulate(self) -> list \n\n"
                         "simulate one vector")

        .def("Simulate", WRAP::Dmd_wrapper_simulate1,
                         return_value_policy< manage_new_object >(),
                         "Simulate(self, int) -> Vectors \n\n"
                         "simulate several vectors")

        /* DEF_RETURN_VALUE("Simulate", WRAP::Dmd_wrapper_simulate1,
                         args("Size"),
                         "Simulate(self, Size) -> list \n\n") */

        // .def(self_ns::str(self)) //__str__
        .def("__str__", WRAP::Dmd_wrapper_ascii_write0)
    ;

    def("DmdAsciiRead", WRAP::Dmd_wrapper_ascii_read1,
                        return_value_policy< manage_new_object >());
        def("OrdinalAssociationModel", WRAP::Dmd_wrapper_ordinal_association,
                                                    "OrdinalAssociationModel(self) -> int \n\n"
                                                    "Linear by linear model for association");

};

#undef WRAP

#define WRAP DiscreteMultivariateParametricWrap
class WRAP
{

public :

   static boost::shared_ptr<DiscreteMultivariateParametric>
   Dmp_wrapper_ascii_read(const char * path, double cumul_threshold = CUMUL_THRESHOLD)
   {
      StatError error;
      ostringstream error_message;
      DiscreteMultivariateDistribution *sdist = NULL;
      DiscreteMultivariateParametric *dist = NULL;

      sdist = discrete_multivariate_ascii_read(error, path,
                                               cumul_threshold);

      if (sdist == NULL)
      {
         tree_statistic::throw_stat_tree_error(error);
         /* error_message << error;
         tree_statistic::tree_statistic::throw_stat_tree_error(error_message); */
      }
      else
      {
         dist = dynamic_cast<DiscreteMultivariateParametric*>(sdist);
         if (dist == NULL)
         {
            error_message << "Not a DiscreteMultivariateParametric";
            tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            // tree_statistic::tree_statistic::throw_stat_tree_error(error_message);
         }
      }

      return boost::shared_ptr<DiscreteMultivariateParametric>(dist);
   }


   static boost::shared_ptr<DiscreteMultivariateParametric>
   Dmp_wrapper_ascii_read1(const char * path)
   { return Dmp_wrapper_ascii_read(path); }

   static boost::shared_ptr<DiscreteMultivariateParametric>
   Dmp_wrapper_build_param(int inb_variable, int iident, unsigned int iinf_bound,
                           boost::python::list isup_bound,
                           boost::python::list vparameter,
                           boost::python::list vprobability,
                           double cumul_threshold = CUMUL_THRESHOLD)
   {
      StatError error;
      DiscreteMultivariateParametric *dist = NULL;
      std::vector<unsigned int> lsup_bound;
      std::vector<double> lparameter, lprobability;
      ostringstream error_message;
      unsigned int plength;
      int i = 0;

      if (inb_variable < 2)
      {
         error_message << "Bad number of variables";
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      switch (iident)
      {
         case MNEGATIVE_MULTINOMIAL:
         {
            double rparameter = 1.0;
            if (boost::python::len(isup_bound) > 0)
            {
               error_message << "Bad number of superior bounds";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }

            if (boost::python::len(vparameter) != 1)
            {
               error_message << "Bad number of parameters";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            else
            {
               lparameter.resize(1);
               lparameter[0] = boost::python::extract<double>(vparameter[0]);
            }
            if (lparameter[0] <= 0.)
            {
               error_message << "Bad value of parameter";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }

            plength = boost::python::len(vprobability);

            if (plength != inb_variable)
            {
               error_message << "Bad number of probabilities";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            else
            {
               lprobability.resize(plength);

               for(i = 0; i < plength; i++)
               {
                  lprobability[i] = boost::python::extract<double>(vprobability[i]);
                  rparameter -= (1-lprobability[i]);
                  if ((lprobability[i] > 1.) || (lprobability[i] < 0.))
                  {
                     error_message << "Bad value for probability " << i+1;
                     tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
                  }
               }
            }
            if ((rparameter >= 1.) || (rparameter <= 0.))
            {
               error_message << "Bad values for probabilities";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            break;
         }
         case MPOISSON :
         {
            if (boost::python::len(isup_bound) > 0)
            {
                error_message << "Bad number of superior bounds";
                tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            if (boost::python::len(vparameter) != inb_variable+1)
            {
                error_message << "Bad number of parameters";
                tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            else
            {
               lparameter.resize(inb_variable+1);
               for(i = 0; i < inb_variable+1; i++)
               {
                  lparameter[i] = boost::python::extract<double>(vparameter[i]);
                  if (lparameter[i] <= 0)
                  {
                     error_message << "Bad value for parameter " << i+1;
                     tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
                  }
               }
            }
            if(boost::python::len(vprobability) > 0)
            {
                error_message << "Multivariate Poisson distributions are not defined "
                              << "by probabilities but by parameters";
                tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            break;
         }
         case MMULTINOMIAL :
         {
            if (boost::python::len(isup_bound) > 0)
            {
               error_message << "Bad number of superior bounds";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }

            if (boost::python::len(vparameter) != 1)
            {
               error_message << "Bad number of parameters";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            else
            {
               lparameter.resize(1);
               lparameter[0] = boost::python::extract<double>(vparameter[0]);
            }
            if (lparameter[0] <= 0.)
            {
               error_message << "Bad value of parameter";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }

            plength = boost::python::len(vprobability);

            if (plength != inb_variable)
            {
               error_message << "Bad number of probabilities";
               tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }
            else
            {
               double rparameter = 0.;

               lprobability.resize(plength);
               for(i = 0; i < plength; i++)
               {
                  lprobability[i] = boost::python::extract<double>(vprobability[i]);
                  if ((lprobability[i] > 1.) || (lprobability[i] < 0.))
                  {
                     error_message << "Bad value for probability " << i+1;
                     tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
                  }
                  else
                     rparameter += lprobability[i];
               }
               // check that probabilities sum to 1.
               if (abs(1-rparameter) > 1-CUMUL_THRESHOLD)
               {
                  error_message << "Bad values for probabilities";
                  tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
               }
            }
            break;
         }
         default:
         {
            error_message << "Bad identifier for Discrete Multivariate Distribution";
            tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
         }
      }

      dist = new DiscreteMultivariateParametric(inb_variable, iident, iinf_bound,
                                                 lsup_bound, lparameter, lprobability,
                                                 cumul_threshold);

      return boost::shared_ptr<DiscreteMultivariateParametric>(dist);
   }

   static boost::shared_ptr<DiscreteMultivariateParametric>
   Dmp_wrapper_build_param0(int inb_variable, int iident, unsigned int iinf_bound,
                            boost::python::list isup_bound,
                            boost::python::list vparameter,
                            boost::python::list vprobability)
   { return Dmp_wrapper_build_param(inb_variable, iident, iinf_bound, isup_bound,
                                    vparameter, vprobability); } // CUMUL_THRESHOLD }

   static str Dmp_wrapper_ascii_write0(const DiscreteMultivariateParametric& dist)
   {
      std::stringstream s;
      str res;

      dist.line_write(s);
      res = str(s.str());
      return res;
   }


   static bool Dmp_wrapper_file_ascii_write2(const DiscreteMultivariateParametric& dist,
                                             const char * path, bool exhaustive)
   {
      StatError error;
      bool res;

      res = dist.ascii_write(error, path, exhaustive);
      if (!res)
         tree_statistic::throw_stat_tree_error(error);
      return res;
   }

   static bool Dmp_wrapper_file_ascii_write1(const DiscreteMultivariateParametric& dist,
                                             const char * path)
   { return Dmp_wrapper_file_ascii_write2(dist, path, false); }

   static MultiPlotSet*
   Dmp_wrapper_get_plotable(DiscreteMultivariateParametric& dist)
   {
     StatError error;
     MultiPlotSet *ret = dist.get_plotable();
     if (ret == NULL)
       stat_tool::wrap_util::throw_error(error);
     return ret;
   }

   static DiscreteMultivariateParametric*
   Dmd_estimate(const Vectors& vec, int ident, unsigned int min_inf_bound, bool min_inf_bound_flag)
   {
      unsigned int i, v, j;
      double loglikelihood;
      ostringstream error_message;
      DiscreteMultivariateReestimation<int> *histo = NULL;
      DiscreteMultivariateParametric *res = NULL;

      if (vec.get_nb_vector() < 2)
      {
         error_message << "Bad number of vectors: " << vec.get_nb_vector();
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      // compute histogram
      histo = new DiscreteMultivariateReestimation<int>(vec);

      if (histo->get_nb_variable() < 2)
      {
         error_message << "Bad number of variables: " << histo->get_nb_variable();
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      res = new DiscreteMultivariateParametric();

      // estimate distribution
      switch(ident)
      {
         case MNEGATIVE_MULTINOMIAL :
         {
            loglikelihood
               = histo->negative_multinomial_estimation(res, min(min_inf_bound, histo->get_offset()), min_inf_bound_flag,
                                                        CUMUL_THRESHOLD);
            break;
         }
         case MPOISSON :
         {
            loglikelihood
               = histo->mpoisson_estimation_ML(res, min(min_inf_bound, histo->get_offset()),
                                               min_inf_bound_flag, CUMUL_THRESHOLD);
                        break;
         }
         case MMULTINOMIAL :
         {
             loglikelihood
                = histo->mmultinomial_estimation(res, min(min_inf_bound, histo->get_offset()),
                                                 min_inf_bound_flag, CUMUL_THRESHOLD);
             break;
         }
         default:
         {
            error_message << "Unknown type of discrete multivariate distribution: " << ident;
            tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
         }
      }

      delete histo;
      histo = NULL;

      if (loglikelihood <= D_INF)
      {
         error_message << "Estimation of discrete multivariate distribution failed";
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      return res;
   }

   static DiscreteMultivariateParametric*
   Dmd_nmultinomial_estimate(const Vectors& vec, unsigned int min_inf_bound, bool min_inf_bound_flag)
   {
      unsigned int i, v, j;
      double loglikelihood;
      ostringstream error_message;
      DiscreteMultivariateReestimation<int> *histo = NULL;
      DiscreteMultivariateParametric *res = NULL;

      if (vec.get_nb_vector() < 2)
      {
         error_message << "Bad number of vectors: " << vec.get_nb_vector();
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      // compute histogram
      histo = new DiscreteMultivariateReestimation<int>(vec);

      if (histo->get_nb_variable() < 2)
      {
         error_message << "Bad number of variables: " << histo->get_nb_variable();
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      res = new DiscreteMultivariateParametric();

      // estimate distribution
      loglikelihood = histo->negative_multinomial_estimation(res, min(min_inf_bound, histo->get_offset()), min_inf_bound_flag,
                                                             CUMUL_THRESHOLD);

      delete histo;
      histo = NULL;

      if (loglikelihood <= D_INF)
      {
         error_message << "Cannot estimate negative multinomial distribution";
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      return res;
   }

}; // WRAP


void class_dmparametric()
{
  class_< DiscreteMultivariateParametric, bases< DiscreteMultivariateDistribution > >
    ("_DiscreteMultivariateParametric", "Parametric Discrete Multivariate Distributions", init< const DiscreteMultivariateParametric& >())
        .def("__init__", make_constructor(WRAP::Dmp_wrapper_ascii_read))
        .def("__init__", make_constructor(WRAP::Dmp_wrapper_ascii_read1))
        .def("__init__", make_constructor(WRAP::Dmp_wrapper_build_param0))
        .def("__init__", make_constructor(WRAP::Dmp_wrapper_build_param))

        .def("FileAsciiWrite", WRAP::Dmp_wrapper_file_ascii_write1)
        .def("FileAsciiWrite", WRAP::Dmp_wrapper_file_ascii_write2)

        DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::Dmp_wrapper_get_plotable,
                                 "return plotable")

        .def_readonly("get_parameter", &DiscreteMultivariateParametric::get_parameter,
                      "Distribution parameters")
        .def_readonly("get_probability", &DiscreteMultivariateParametric::get_probability,
                      "Distribution probabilities")

        // .def(self_ns::str(self)) //__str__
        .def("__str__", WRAP::Dmp_wrapper_ascii_write0)
    ;

    def("NegativeMultinomialEstimation", WRAP::Dmd_nmultinomial_estimate,
                                         return_value_policy< manage_new_object >(),
                                         "NegativeMultinomialEstimation(self, int, bool) -> _DiscreteMultivariateParametric \n\n"
                                         "Estimate a negative multinomial distribution");

    def("MultivariateParametricEstimation", WRAP::Dmd_estimate,
                                            return_value_policy< manage_new_object >(),
                                            "MultivariateParametricEstimation(self, int, int, bool) -> _DiscreteMultivariateParametric \n\n"
                                            "Estimate a parametric discrete multivariate distribution");
};
#undef WRAP

#define WRAP MultinomialCompoundWrap
class WRAP
{

public :

   static boost::shared_ptr<MultinomialCompoundDiscreteParametric>
   Mcdp_wrapper_build_distrib1(int inb_variable,
                               boost::python::list iprobability,
                               const DiscreteParametric& dist,
                               double cumul_threshold)
   {
      unsigned int i;
      double rparameter = 0.; // sum of probabilities
      std::vector<double> lprobability;
      ostringstream error_message;
      MultinomialCompoundDiscreteParametric *rdist = NULL;

      lprobability.resize(inb_variable);
      for(i = 0; i < inb_variable; i++)
      {
         lprobability[i] = boost::python::extract<double>(iprobability[i]);
         if ((lprobability[i] > 1.) || (lprobability[i] < 0.))
         {
            error_message << "Bad value for probability " << i+1;
            tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
         }
         else
            rparameter += lprobability[i];
      }
      // check that probabilities sum to 1.
      if (abs(1-rparameter) > 1-CUMUL_THRESHOLD)
      {
         error_message << "Bad values for probabilities";
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      rdist = new MultinomialCompoundDiscreteParametric(inb_variable,
                                                        lprobability, dist, cumul_threshold);
       if (rdist == NULL)
       {
          error_message << "Error";
          tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
       }
       return boost::shared_ptr<MultinomialCompoundDiscreteParametric>(rdist);
   }

   static boost::shared_ptr<MultinomialCompoundDiscreteParametric>
   Mcdp_wrapper_build_distrib0(int inb_variable, boost::python::list iprobability,
                               const DiscreteParametric& dist)
   {return Mcdp_wrapper_build_distrib1(inb_variable, iprobability,
                                        dist, CUMUL_THRESHOLD); }

   static MultinomialCompoundDiscreteParametric*
   Mcdp_estimate(const Vectors& vec, int ident, unsigned int min_inf_bound, bool min_inf_bound_flag)
   {
      unsigned int i, v, j;
      double loglikelihood;
      ostringstream error_message;
      DiscreteMultivariateReestimation<int> *histo = NULL;
      MultinomialCompoundDiscreteParametric *res = NULL;

      if (vec.get_nb_vector() < 2)
      {
         error_message << "Bad number of vectors: " << vec.get_nb_vector();
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      // compute histogram
      histo = new DiscreteMultivariateReestimation<int>(vec);

      if (histo->get_nb_variable() < 2)
      {
         error_message << "Bad number of variables: " << histo->get_nb_variable();
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      res = new MultinomialCompoundDiscreteParametric();

            switch(ident){
                              case POISSON :
                                case NEGATIVE_BINOMIAL :
                case BINOMIAL : {
                    loglikelihood = histo->compound_multinomial_estimation(res, ident, min(min_inf_bound, histo->get_offset()), min_inf_bound_flag, CUMUL_THRESHOLD);
            break;
                }
                default :
                  error_message << "Unknown type of compound multinomial discrete distribution: " << ident;
          tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
            }


      delete histo;
      histo = NULL;

      if (loglikelihood <= D_INF)
      {
         error_message << "Cannot estimate negative multinomial distribution";
         tree_statistic::throw_stat_tree_error((error_message.str()).c_str());
      }

      return res;
   }

};

void class_mcdparametric()
{
  class_< MultinomialCompoundDiscreteParametric, bases< DiscreteMultivariateParametric > >
    ("_MultinomialCompoundDiscreteParametric", "Parametric Discrete Compound Multinomial Distributions", init< const MultinomialCompoundDiscreteParametric& >())
       .def("__init__", make_constructor(WRAP::Mcdp_wrapper_build_distrib0))
       .def("__init__", make_constructor(WRAP::Mcdp_wrapper_build_distrib1))
    ;

    def("CompoundMultinomialEstimation", WRAP::Mcdp_estimate,
                                         return_value_policy< manage_new_object >(),
                                         "CompoundMultinomialEstimation(self, int, int, bool) -> _MultinomialCompoundDiscreteParametric\n\n"
                                         "Estimate a compound multinomial distribution");

};
#undef WRAP

