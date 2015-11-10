/****************************************************************
 *
 *  Test of the data structure and algorithms for Markov
 *  out-trees as defined in markov_out_tree.h
 */

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"

// required by "tree_statistic/multivariate_distribution.h"
#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_label.h"
#include "tree_statistic/tree_labels.h"
#include "statiskit/core/data/marginal/multivariate.h"
#include "tree_statistic/mv_histogram_tools.h"

#include "tree_statistic/multivariate_distribution.h"

using namespace stat_tool;
using namespace Stat_trees;


int main(void)
{
   // typedef GenerationProcess::index_array index_array;

   const int nb_variable = 3, nb_vector = 12;
   unsigned int i;

   const unsigned int inf_bound = 1;
   std::vector<unsigned int> sup_bound(nb_variable);
   std::vector<double> parameter, probability;
   std::vector<unsigned int> value(nb_variable, 0);

   DiscreteMultivariateDistribution *tmp_distr = NULL;
   DiscreteMultivariateParametric *tmp_param = NULL;
   Statiskit::Marginal::Multivariate::CountTable *mv_hist = NULL;
   Statiskit::Marginal::Univariate::Table<int> *sum_mv_hist = NULL;
   FrequencyDistribution *fd = NULL;

   probability.resize(nb_variable);
   probability[0] = 0.5;
   probability[1] = 0.25;
   probability[2] = 0.25;
   sup_bound[0] = 0;
   sup_bound[1] = 0;
   sup_bound[2] = 0;
   parameter.resize(1);
   parameter[0] = 12;

   tmp_param = new DiscreteMultivariateParametric(nb_variable, MMULTINOMIAL,
                                                  inf_bound, sup_bound,
                                                  parameter, probability,
                                                  CUMUL_THRESHOLD);
   tmp_param->ascii_print(cout, false, false, false, false, NULL);
   tmp_distr = tmp_param;
   assert(tmp_distr = tmp_param);
   tmp_distr->ascii_print(cout, false, false, false, false, NULL);
   delete tmp_param;
   tmp_param = NULL;

   mv_hist = new Statiskit::Marginal::Multivariate::CountTable(Statiskit::storages_type(nb_variable, Statiskit::INTEGER));

   for(i=0; i < nb_vector-1; i++)
   {
      value[0] = i;
      value[1] = i / 2;
      mv_hist->add(Stat_histogram_value(value), 1);
   }
   value[0] = nb_vector-2;
   value[1] = (nb_vector-2)/2;
   mv_hist->add(Stat_histogram_value(value), 1);

   cout << "Statiskit::Marginal::Multivariate::CountTable: " << endl;
   print_stat_histogram(*mv_hist, cout);

   sum_mv_hist = Stat_trees::Stat_pseudo_histogram_get_sum(*mv_hist, Stat_histogram_value_sum_int);
   delete mv_hist;
   mv_hist = NULL;

   fd = Stat_trees::Scalar_stat_histogram_to_frequency_distribution(*sum_mv_hist);
   delete sum_mv_hist;
   sum_mv_hist = NULL;

   cout << endl << "FrequencyDistribution for components sum" << endl;
   fd->ascii_write(cout, true, true);
   delete fd;
   fd = NULL;

   return 0;
}
