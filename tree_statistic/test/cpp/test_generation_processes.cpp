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

#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"   // definition of DiscreteParametricModel class
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/sequence_label.h"

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"

#include "tree_statistic/int_fl_containers.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include "tree_statistic/typed_edge_trees.h"


#include "tree_statistic/hidden_markov_tree.h"
#include "tree_statistic/markov_out_tree.h"

using namespace Stat_trees;


int main(void)
{
   typedef GenerationProcess::index_array index_array;

   register int i, t;
   const int nb_state = 3, nb_output_process = 2, nb_trees = 1;
   const unsigned int nb_children_branching = 1, nb_factors = 0,
                      nb_generation = 2, nb_generation_values = nb_state * nb_state;

   unsigned int inf_bound = 1, total_size;
   double sum = 0.0;
   std::vector<unsigned int> sup_bound, depths(1);
   std::vector<double> parameter, probability, entropy;
   int *sel_var = new int[1];

   index_array *factor_values = NULL;
   index_array generation_values(0);
   GenerationProcess *pgeneration_distributions = NULL;
   DiscreteMultivariateDistribution **pgeneration = new DiscreteMultivariateDistribution*[nb_generation_values];
   DiscreteMultivariateParametric **pgeneration_param = new DiscreteMultivariateParametric*[nb_generation_values];
   DiscreteMultivariateDistribution *tmp_distr = NULL;
   DiscreteMultivariateParametric *tmp_param = NULL;
   MultinomialCompoundDiscreteParametric *tmp_comp = NULL;
   DiscreteParametric *cdist = NULL;


   StatError error;

   MultiPlotSet *plot = NULL;
   std::vector<DiscreteDistributionData*> distrib_data(0);

   depths[0] = 5;


   // compound multinomials
   probability.resize(nb_state);
   parameter.resize(1);
   for(i = 0; i < nb_state; i++)
       probability[i] = 1. / nb_state;

   for(i = 0; i < nb_generation_values; i++)
   {
      cdist = new DiscreteParametric(BINOMIAL, inf_bound, i, D_DEFAULT, 0.5);
      tmp_comp = new MultinomialCompoundDiscreteParametric(nb_state, probability, *cdist, CUMUL_THRESHOLD);
      delete cdist;
      cdist = NULL;
      tmp_distr = tmp_comp;
      assert(tmp_comp == tmp_distr);
      cout << tmp_comp << endl;
      cout << tmp_distr << endl;
      tmp_comp->ascii_print(cout, false, false, false, false, NULL);
      tmp_distr->ascii_print(cout, false, false, false, false, NULL);
      delete tmp_comp;
   }

   probability.resize(nb_state);
   parameter.resize(1);
   sum = 0.;
   for(i = 0; i < nb_state; i++)
   {
       probability[i] = 1. - 1./(nb_state + 1);
       sum += (1. - probability[i]);
   }
   assert(sum < 1.);

   for(i = 0; i < nb_generation_values; i++)
   {
      parameter[0] = i + 1;
      tmp_param = new DiscreteMultivariateParametric(nb_state, MNEGATIVE_MULTINOMIAL,
                                                     inf_bound, sup_bound, parameter, probability);
      tmp_distr = tmp_param;
      assert(tmp_distr == tmp_param);
      cout << tmp_distr << endl;
      cout << tmp_param << endl;
      tmp_param->ascii_print(cout, false, false, false, false, NULL);
      tmp_distr->ascii_print(cout, false, false, false, false, NULL);
      pgeneration[i] = tmp_param->copy_ptr();
      pgeneration_param[i] = tmp_param->copy_ptr();
      delete tmp_param;
   }


   pgeneration_distributions = new GenerationProcess(nb_state,
                                                     nb_children_branching,
                                                     nb_factors, factor_values,
                                                     nb_generation, pgeneration);
   // test code and decode
   assert(pgeneration_distributions->code(pgeneration_distributions->decode(nb_generation_values-1))
           == nb_generation_values - 1);

   for(i = 0; i < nb_generation_values; i++)
   {
      delete pgeneration[i];
      pgeneration[i] = NULL;
      delete pgeneration_param[i];
      pgeneration_param[i] = NULL;
   }
   delete [] pgeneration;
   delete [] pgeneration_param;

   // test get_plotable for MarkovOutTreeData

   // test get_distribution_data
   // nb_state = mtd->get_nb_states();


   return 0;
}
