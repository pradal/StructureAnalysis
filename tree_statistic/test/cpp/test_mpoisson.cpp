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
     std::vector<double> parameter, probability;
     DiscreteMultivariateReestimation<int> *reestimation;
     std::vector<index_array> sample;
     StatError error;
     const char * mpoisson_distrib = "./mpoisson.dmd";
     MultiPlotSet *plot = NULL;

     index_array indices, min_indices;
     IidDiscreteMultivariateParametric* iid_binomial = NULL;
     DiscreteMultivariateParametric *nmultinomial = NULL, *nmpoisson = NULL, *mdist = NULL, *emdist = NULL;
     DiscreteMultivariateDistribution *multivariate_distribution = NULL;

     parameter.resize(5);
     parameter[0] = 0.5;
     parameter[1] = 0.5;
     parameter[2] = 0.4;
     parameter[3] = 0.2;
     parameter[4] = 0.3;
     nmpoisson = new DiscreteMultivariateParametric(nb_variable, MPOISSON, inf_bound, sup_bound, parameter, probability);
     nmpoisson->line_write(cout);
     cout << endl;
     nmpoisson->ascii_write(cout, true);

     sample.resize(sample_size);
     // simulation
     cout << "simulation: " << endl;
     for(i = 0; i < sample_size; i++){
         indices = nmpoisson->simulation();
         if(i < 10)
             cout << indices << "(" << nmpoisson->get_mass(indices) << ")" << endl;
         sample[i] = indices;
         it = std::max_element(indices.begin(), indices.end());
         max_value = MAX(max_value, *it);
   }



     reestimation = new DiscreteMultivariateReestimation<int>(nb_variable, inf_bound, max_value-inf_bound+1);
     for(i = 0; i < sample_size; i++){
         reestimation->update(sample[i], 1);
     }

     // log likelihood
     cout << "Log likelihood of true parameter: " << reestimation->likelihood_computation(*nmpoisson) << endl;

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

     mdist = new DiscreteMultivariateParametric();
     likelihood = reestimation->mpoisson_estimation_M(mdist, 0, true, CUMUL_THRESHOLD);
     if (likelihood > D_INF)
     {
        cout << "Estimate distribution: " << endl;
        mdist->line_write(cout);
        cout << endl;
        mdist->ascii_write(cout, true);
     }
     else
        cout << "Cannot estimate distribution" << endl;
     delete mdist;

     emdist = new DiscreteMultivariateParametric();
     likelihood = reestimation->mpoisson_estimation_ML(emdist, 0, true, CUMUL_THRESHOLD);
     if (likelihood > D_INF)
     {
        cout << "Estimate distribution: " << endl;
        emdist->line_write(cout);
        cout << endl;
        emdist->ascii_write(cout, true);
     }
         else
             cout << "Cannot estimate distribution" << endl;
     delete emdist;

     // read distribution from file

     // multivariate Poisson
     multivariate_distribution = discrete_multivariate_ascii_read(error, mpoisson_distrib);
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
