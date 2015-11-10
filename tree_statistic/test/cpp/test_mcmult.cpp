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
   typedef DiscreteMultivariateDistribution::index_array index_array;
   index_array::const_iterator it;
   unsigned int i, nb_variable = 4, var, counts, inf_bound = 2, max_value = 0, sample_size = 1000;
   double sum, cumul, proba, likelihood, eparameter;
   std::vector<unsigned int> sup_bound;
   std::vector<int>  sums;
   std::vector<double> probability;
   double parameter;
   DiscreteMultivariateReestimation<int> *reestimation;
   std::vector<index_array> sample;
   StatError error;
   const char * mcompound_distrib = "./mcompound.dmd";
   MultiPlotSet *plot = NULL;
   DiscreteMultivariateDistribution *multivariate_distribution = NULL;
   DiscreteParametric *cdist = NULL;

   index_array indices, min_indices;

   MultinomialCompoundDiscreteParametric *cmultinomial = NULL, *edist = NULL;
   probability.resize(4,0);
   probability[0] = 20;
   probability[1] = 10;
   probability[2] = 5;
   probability[3] = 5;
   cdist = new DiscreteParametric(BINOMIAL, inf_bound, 10, D_DEFAULT, 0.5);

   cmultinomial = new MultinomialCompoundDiscreteParametric(4, probability, *cdist, CUMUL_THRESHOLD);

   multivariate_distribution = cmultinomial;
   assert(multivariate_distribution == cmultinomial);

   cout << "Pointer to a MultinomialCompoundDiscreteParametric" << endl;
   cmultinomial->line_write(cout);
   cout << endl;
   cmultinomial->ascii_write(cout, true);

   cout << "Copy pointer into a MultivariateDistribution" << endl;
   multivariate_distribution->ascii_write(cout, true);

   sample.resize(sample_size);
   //simulation
   cout << "simulation: " << endl;
   for(i = 0; i < sample_size; i++){
       indices = cmultinomial->simulation();
           if(i < 10)
               cout << indices << "(" << cmultinomial->get_mass(indices) << ")" << endl;
           sample[i] = indices;
           it = std::max_element(indices.begin(), indices.end());
           max_value = MAX(max_value, *it);
   }

   reestimation = new DiscreteMultivariateReestimation<int>(nb_variable, 0, max_value+1);

   for(i = 0; i < sample_size; i++){
       reestimation->update(sample[i], 1);
   }

   //log likelihood
   //cout << "Log likelihood of true parameter: " << reestimation->likelihood_computation(*cmultinomial) << endl;

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
   cout << "Sums : ";
   counts = 0;
   for(var = 0; var < nb_variable; var++){
       counts += sums[var];
       cout << sums[var] << " " ;
   }
   cout << endl;
   cout << "Marginals: " << endl;
   for(var = 0; var < nb_variable; var++){
       reestimation->get_marginal_ptr(var)->print(cout);
   }
   edist = new MultinomialCompoundDiscreteParametric();
   reestimation->compound_multinomial_estimation(edist, BINOMIAL, 0, true, CUMUL_THRESHOLD);
   cout << "Estimate distribution: " << endl;
   edist->line_write(cout);
   cout << endl;
   edist->ascii_write(cout, false);
   delete edist;

   // read distribution from file

   cout << "Read compound multinomial from a file: " << endl;
   multivariate_distribution = discrete_multivariate_ascii_read(error, mcompound_distrib);
   if (error.get_nb_error() > 0)
   {
      cout << error;
      return 1;
   }

   multivariate_distribution->ascii_write(cout);

   delete multivariate_distribution;
   multivariate_distribution = NULL;

   return 0;
}
