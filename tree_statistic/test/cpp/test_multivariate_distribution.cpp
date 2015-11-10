/****************************************************************
 *
 *  Test of the MultivariateDiscreteDistribution methods
 *  as defined in multivariate_distributions.h
 */

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"

// required by "tree_statistic/multivariate_distribution.h"
#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_label.h"
#include "tree_statistic/tree_labels.h"
#include "statiskit/core/data/marginal/multivariate.h"

#include "tree_statistic/multivariate_distribution.h"
#include "tree_statistic/mv_histogram_tools.h"
#include <iostream>
#include <string>
#include <cstdlib>

#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>

#include "assert.h"

#include <boost/math/special_functions/gamma.hpp>

using namespace stat_tool;
using namespace Stat_trees;

using namespace boost::math;

int main(void)
{
   typedef DiscreteMultivariateReestimation<int>::value value;
   typedef DiscreteMultivariateDistribution::index_array index_array;
   index_array::const_iterator it;

   const bool test_loop = false, test_gamma = false, test_cumul = false;
   // const bool test_loop = true, test_gamma = true, test_cumul = true;
   unsigned int i, nb_variable = 4, var, counts,
                inf_bound = 0, max_value = 0,
                sample_size = 2000;
   double sum, cumul, proba, likelihood, eparameter;
   std::vector<unsigned int> sup_bound;
   std::vector<int>  sums;
   std::vector<double> parameter, probability;
   DiscreteMultivariateReestimation<int> *reestimation;
   std::vector<index_array> sample;
   std::vector<value> elements;
   Statiskit::storages_type scalars(nb_variable, Statiskit::INTEGER);
   StatError error;
   const char * outputpath = "./output.dmd";
   const char * multinomial_distrib = "./multi.dmd";
   const char * vec_path = "./nmvectors.vec";
   MultiPlotSet *plot = NULL;
   Vectors *vec = NULL;
   Statiskit::Marginal::Multivariate::CountTable *mv_histo = NULL;

   index_array indices, min_indices;

   IidDiscreteMultivariateParametric* iid_binomial = NULL;
   DiscreteMultivariateParametric *nmultinomial = NULL, *edist = NULL;
   DiscreteMultivariateDistribution *multivariate_distribution = NULL;

   multivariate_distribution = new DiscreteMultivariateDistribution();
   cout << "Empty DiscreteMultivariateDistribution: " << endl;
   multivariate_distribution->ascii_write(cout);
   cout << endl;
   delete multivariate_distribution;
   multivariate_distribution = NULL;

   iid_binomial = new IidDiscreteMultivariateParametric(BINOMIAL, nb_variable,
                                                        2, 3, D_DEFAULT, 0.3);

   indices = iid_binomial->loop_init();
   min_indices.resize(indices.size());

   if (test_loop)
      for(i = 1; i <= 257; i++)
      {
         cout << "indices : " << indices << endl;
         indices = iid_binomial->loop_next_low();
      }

   iid_binomial->line_write(cout);
   cout << endl;
   cout << "Ascii write (not exhaustive): " << endl;
   iid_binomial->ascii_write(cout, false);
   cout << "Ascii write (exhaustive): " << endl;
   iid_binomial->ascii_write(cout, true);
   // write distribution into a file
   iid_binomial->ascii_write(error, outputpath, true);

   // read distribution from file
   multivariate_distribution = discrete_multivariate_ascii_read(error, outputpath);
   if (error.get_nb_error() > 0)
   {
      cout << error;
      return 0;
   }
   cout << "Read IidDiscreteMultivariateParametric from file: " << endl;
   multivariate_distribution->ascii_write(cout);
   delete multivariate_distribution;
   multivariate_distribution = NULL;


   mv_histo = new Statiskit::Marginal::Multivariate::CountTable(scalars);

   // simulation
   cout << "simulation: " << endl;
   for(i = 0; i < 10; i++)
   {
      indices = iid_binomial->simulation();
      cout << indices << "(" << iid_binomial->get_mass(indices) << ")" << endl;
      mv_histo->add(Stat_histogram_value(indices), 1);
   }

   reestimation = stat_histogram_to_DiscreteMultivariateReestimation(*mv_histo);
   // histogram
   cout << "Histogram: " << endl;
   cout << "Number of different values: " << reestimation->nb_value_computation() << endl;
   reestimation->offset_computation();
   reestimation->max_computation();
   reestimation->mean_computation();
   reestimation->variance_computation();
   reestimation->print(cout);

   // elements of histogram
   cout << "Elements of histogram: " << endl;
   elements = reestimation->get_elements();
   assert(elements.size() > 0);
   for(i = 0; i < elements.size(); i++)
      cout << elements[i] << ": " << reestimation->get_frequency(elements[i]) << endl;

   for(i = 0; i < elements[0].size(); i++)
      elements[0][i] -= 1;

   cout << "Frequency for element " << elements[0] << ": ";
   cout << reestimation->get_frequency(elements[0]) << endl;

   for(i = 0; i < elements[0].size(); i++)
      elements[0][i] += 2;

   cout << "Frequency for element " << elements[0] << ": ";
   cout << reestimation->get_frequency(elements[0]) << endl;

   delete iid_binomial;
   iid_binomial = NULL;
   delete reestimation;
   reestimation = NULL;

   probability.resize(nb_variable);
   indices.resize(nb_variable);
   parameter.resize(1);
   parameter[0] = 13.0;
   sum = 0.;
   for(i = 0; i < nb_variable; i++)
   {
       probability[i] = 1. - 1./(nb_variable+1);
       sum += (1. - probability[i]);
       min_indices[i] = 2;
   }
   assert(sum < 1.);
   // min_indices[nb_variable-1] -= 1;

   // test gamma function
   if (test_gamma)
   {
      cout << "test gamma function: " << endl;

      cout << "tgamma(3): " << tgamma(3) << endl;

      cout << "tgamma(3.5): " << tgamma(3.5) << endl;

      cout << "tgamma_delta_ratio(4.8, -1.2): " << tgamma_delta_ratio(4.8, -1.2)  << endl;

      cout << "... to be compared with: " << tgamma(4.8) / tgamma(3.6) << endl;
   }

   nmultinomial = new DiscreteMultivariateParametric(nb_variable, MNEGATIVE_MULTINOMIAL,
                                                     inf_bound, sup_bound, parameter, probability);

   nmultinomial->line_write(cout);
   cout << endl;
   nmultinomial->ascii_write(cout, true);

   indices = nmultinomial->loop_init(min_indices);
   cumul = 0.;

   for(i = 0; i < NB_VALUE; i++)
   {
      proba = nmultinomial->get_mass(indices);
      cumul += proba;
      if (i < 10)
         cout << indices << " : " << proba << endl;
      indices = nmultinomial->loop_next_low(min_indices);
   }
   cout << "cumul. prob. ( " << NB_VALUE << "): " << cumul << endl;

   if (test_cumul)
   {
      cumul = 0.;
      cout << "nb. necessary values to reach CUMUL_THRESHOLD: " << endl;
      while (cumul < CUMUL_THRESHOLD)
      {
         i++;
         proba = nmultinomial->get_mass(indices);
         cumul += proba;
         indices = nmultinomial->loop_next_low(min_indices);
         if ((i % 500000) == 0)
            cout << i << " iterations (" << indices << ")" << endl;
      }

      cout << "nb. val. ( " << CUMUL_THRESHOLD << "): " << i << endl;

      // f.y.i. 4,100,625
   }

   // simulation
   cout << "simulation: " << endl;
   for(i = 0; i < 10; i++)
   {
      indices = nmultinomial->simulation();
      cout << indices << "(" << nmultinomial->get_mass(indices) << ")" << endl;
   }

   sample.resize(sample_size);

   // simulation and estimation
   cout << "simulation: " << endl;
   for(i = 0; i < sample_size; i++)
   {
      indices = nmultinomial->simulation();
      if (sample_size < 100)
         cout << indices << endl;
      sample[i] = indices;
      it = std::max_element(indices.begin(), indices.end());
      max_value = MAX(max_value, *it);
   }


   reestimation = new DiscreteMultivariateReestimation<int>(nb_variable, inf_bound, max_value-inf_bound+1);
   for(i = 0; i < sample_size; i++)
      reestimation->update(sample[i], 1);


   // log likelihood
   cout << "Log likelihood of true parameter: "
        << reestimation->likelihood_computation(*nmultinomial) << endl;

   // histogram
   cout << "Histogram: " << endl;
   cout << "Number of different values: " << reestimation->nb_value_computation() << endl;
   reestimation->offset_computation();
   reestimation->max_computation();
   reestimation->mean_computation();
   reestimation->variance_computation();
   sums = reestimation->sum_computation();

   cout << "Number of elements: " << reestimation->get_nb_element() << endl;
   cout << "Offset: " << reestimation->get_offset() << endl;
   cout << "Sums - offset: ";
   counts = 0;
   for(var = 0; var < nb_variable; var++)
   {
      sums[var] -= reestimation->get_offset() * reestimation->get_nb_element();
      counts += sums[var];
      cout << sums[var] << " " ;
   }
   cout << endl;
   cout << "Marginals: " << endl;

   for(var = 0; var < nb_variable; var++)
   {
      reestimation->get_marginal_ptr(var)->print(cout);
   }


   if (sample_size < 100)
      reestimation->print(cout);

   edist = new DiscreteMultivariateParametric();

   // likelihood function as a function of parameter
   # ifdef DEBUG
   if (sample_size < 100)
   {
      eparameter = 0.5;
      for(i = 0; i < 800; i++)
      {
      // likelihood =  reestimation->negative_multinomial_marginal_estimation(edist, 0, true, CUMUL_THRESHOLD);
         eparameter += 0.05 / 10;
         cout << eparameter << " "
              //<< reestimation->negative_multinomial_parameter_estimation_score_target(*nmultinomial,
              //                                                                        sums,     counts, eparameter)
              << reestimation->negative_multinomial_parameter_partial_likelihood(sums, 2,
                                                                                 counts,
                                                                                 eparameter)
              << " "
              << reestimation->negative_multinomial_parameter_estimation_likelihood_target(*nmultinomial,
                                                                                           2, sums, counts,
                                                                                           eparameter)
              << " "
              << reestimation->negative_multinomial_parameter_partial_likelihood_derivative(sums, 2,
                                                                                            counts,
                                                                                            eparameter,
                                                                                            1e-7)

              << " "
              << reestimation->negative_multinomial_parameter_estimation_score_value(*nmultinomial,
                                                                                     2, sums, counts, eparameter)
              << endl;

      }
   }
   #endif

   likelihood =  reestimation->negative_multinomial_estimation(edist, 0, true, CUMUL_THRESHOLD);
   // likelihood =  reestimation->negative_multinomial_estimation(edist, 2, true, CUMUL_THRESHOLD);

   cout << "Estimate distribution: " << endl;
   edist->line_write(cout);
   cout << endl;
   edist->ascii_write(cout, true);

   // plotable
   plot = edist->get_plotable(mv_histo);
   delete plot;
   plot = NULL;

   delete mv_histo;
   mv_histo = NULL;

   cout << "Log likelihood: " << likelihood << endl;
   cout << "Likelihood equation value: "
        << reestimation->negative_multinomial_parameter_estimation_score_target(*edist, edist->get_inf_bound(), sums, counts,
                                                                                edist->get_parameter()[0]) << endl;
   cout << "Empirical means: " << endl;
   cout << reestimation->get_means() << endl;

   // write distribution into a file
   edist->ascii_write(error, outputpath, true);
   delete edist;
   edist = NULL;

   // simulate and estimate from vectors
   vec = nmultinomial->simulation(sample_size);
   assert(vec->get_nb_vector() == sample_size);
   vec->ascii_data_write(error, vec_path);
   delete reestimation;
   reestimation = NULL;
   reestimation = new DiscreteMultivariateReestimation<int>(*vec);
   delete vec;
   vec = NULL;

   edist = new DiscreteMultivariateParametric();

   likelihood =  reestimation->negative_multinomial_estimation(edist, 0, true, CUMUL_THRESHOLD);
   delete reestimation;
   reestimation = NULL;
   cout << endl << "Estimate distribution from vectors: " << endl;
   edist->line_write(cout);
   cout << endl;
   edist->ascii_write(cout, true);


   // read distribution from file
   // negative multinomial
   cout << "Read DiscreteMultivariateParametric from file: " << endl;
   multivariate_distribution = discrete_multivariate_ascii_read(error, outputpath);
   if (error.get_nb_error() > 0)
   {
      cout << error;
      return 0;
   }
   multivariate_distribution->ascii_write(cout);
   delete multivariate_distribution;
   multivariate_distribution = NULL;

   // multinomial
   multivariate_distribution = discrete_multivariate_ascii_read(error, multinomial_distrib);
   if (error.get_nb_error() > 0)
   {
      cout << error;
      return 0;
   }
   delete multivariate_distribution;
   multivariate_distribution = NULL;


   delete edist;
   edist = NULL;
   delete nmultinomial;
   nmultinomial = NULL;

   return 0;

}
