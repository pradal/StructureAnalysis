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
     unsigned int i, nb_variable = 4, var, counts, inf_bound = 2, max_value = 0, sample_size = 3000;
     double sum, cumul, proba, likelihood, eparameter;
     std::vector<unsigned int> sup_bound;
     std::vector<int>  sums;
     std::vector<double> parameter, probability;
     DiscreteMultivariateReestimation<int> *reestimation;
     std::vector<index_array> sample;
     StatError error;
     MultiPlotSet *plot = NULL;

     index_array indices, min_indices;
     DiscreteMultivariateParametric *multinomial = NULL, *edist = NULL;
     DiscreteMultivariateDistribution *multivariate_distribution = NULL;

     probability.resize(4);
     probability[0] = 20;
     probability[1] = 5;
     probability[2] = 4;
     probability[3] = 2;
     parameter.resize(1,24);

     multinomial = new DiscreteMultivariateParametric();
     delete multinomial;

     multinomial = new DiscreteMultivariateParametric(nb_variable, MMULTINOMIAL,
                                                      inf_bound, sup_bound,
                                                      parameter, probability,
                                                      CUMUL_THRESHOLD);
     multinomial->line_write(cout);
     cout << endl;
     multinomial->ascii_write(cout, true);

     sample.resize(sample_size);
     // simulation
     cout << "simulation: " << endl;
     for(i = 0; i < sample_size; i++){
         indices = multinomial->simulation();
         if(i < 10)
             cout << indices << "(" << multinomial->get_mass(indices) << ")" << endl;
         sample[i] = indices;
         it = std::max_element(indices.begin(), indices.end());
         max_value = MAX(max_value, *it);
   }



     reestimation = new DiscreteMultivariateReestimation<int>(nb_variable, inf_bound, max_value-inf_bound+1);
     for(i = 0; i < sample_size; i++){
         reestimation->update(sample[i], 1);
     }

     // log likelihood
     cout << "Log likelihood of true parameter: " << reestimation->likelihood_computation(*multinomial) << endl;

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
     for(var = 0; var < nb_variable; var++){
         sums[var] -= reestimation->get_offset() * reestimation->get_nb_element();
         counts += sums[var];
         cout << sums[var] << " " ;
     }
     cout << endl;
     cout << "Marginals: " << endl;
     for(var = 0; var < nb_variable; var++){
         reestimation->get_marginal_ptr(var)->print(cout);
     }

     edist = new DiscreteMultivariateParametric();
     reestimation->mmultinomial_estimation(edist, 2, true, CUMUL_THRESHOLD);
     cout << "Estimate distribution: " << endl;
     edist->line_write(cout);
     cout << endl;
     edist->ascii_write(cout, true);

     DiscreteMultivariateReestimation<int> *breestimation = NULL;
     breestimation = reestimation->get_bootstrap_distribution(sample_size);
    // log likelihood
     cout << "Log likelihood of true parameter: " << breestimation->likelihood_computation(*multinomial) << endl;

     // histogram
     cout << "Histogram: " << endl;
     cout << "Number of different values: " << breestimation->nb_value_computation() << endl;
     breestimation->offset_computation();
     breestimation->max_computation();
     breestimation->mean_computation();
     breestimation->variance_computation();
     sums = breestimation->sum_computation();

}
