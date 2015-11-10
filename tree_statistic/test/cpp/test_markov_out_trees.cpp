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
#include "tree_statistic/hidden_markov_ind_out_tree.h"
#include "tree_statistic/markov_out_tree.h"

using namespace sequence_analysis;
using namespace Stat_trees;


int main(void)
{
   typedef GenerationProcess::index_array index_array;

   register int i, t;
   const int nb_state = 3, nb_output_process = 2, nb_trees = 1;
   const unsigned int unordered_rank = 2, nb_vomc = 1,
                      nb_children_branching = 1;

   unsigned int inf_bound = 1, total_size, nb_generation_values_fact = 0,
                nb_generation_values = nb_state * nb_state,
                nb_factors = 0, nb_generation = 2;
   double sum = 0.0;
   std::vector<unsigned int> sup_bound, depths(1);
   std::vector<double> parameter, probability, entropy;
   std::vector<GenerationProcess::index_array> comb;
   int *sel_var = new int[1];

   index_array *factor_values = NULL;
   index_array generation_values(0);
   index_array factor_indices(0);
   CategoricalTreeProcess **nptprocess = NULL;
   CategoricalProcess **npprocess = NULL; //, **npprocess_copy = NULL;
   DiscreteParametricProcess **dpprocess = NULL, **dpprocess_copy = NULL;
   ContinuousParametricProcess **cprocess = new ContinuousParametricProcess*[nb_output_process];
   VariableOrderMarkov **pvomc_models = NULL;
   GenerationProcess *pgeneration_distributions = NULL;
   DiscreteMultivariateDistribution **pgeneration = new DiscreteMultivariateDistribution*[nb_generation_values];
   MarkovOutTreeData *mtd = NULL, *mtd_cp = NULL, *sel_mtd = NULL;
   Trees *trees = NULL, *trees_v1 = NULL;

   StatError error;

   const char * vomcpath = "./mot1.vomc";
   const char * hmiotpath = "./hmt.hmt";
   const char * output_path = "./mot_output.mot";
   const char * mot_input_path = "./mot3s_2pop.mot";
   const char * mot_fact_path = "./mot_fact.mot";
   const char * mot_sim_path = "./mot3s_2pop_sim.mot";
   const char * mbp_input_path = "./mbp3s.mot";
   const char * mbp_obs_input_path = "./mbp3s_2pop.mot";
   MarkovOutTree *mot = NULL, *mot_read = NULL, *embp = NULL;
   HiddenMarkovIndOutTree *hmiot = NULL;
   MultiPlotSet *plot = NULL;
   std::vector<DiscreteDistributionData*> distrib_data(0);

   depths[0] = 5;

   pvomc_models = new VariableOrderMarkov*[nb_vomc];
   pvomc_models[nb_vomc-1] = variable_order_markov_ascii_read(error, vomcpath);
   if (pvomc_models[nb_vomc-1] == NULL)
   {
      cout << error;
      return 2;
   }
   hmiot = hidden_markov_ind_out_tree_ascii_read(error, hmiotpath);
   if (hmiot == NULL)
   {
      cout << error;
      return 2;
   }

   nptprocess = hmiot->get_categorical_process();
   dpprocess_copy = hmiot->get_iparametric_process();

   npprocess = new CategoricalProcess*[nb_output_process];
   dpprocess = new DiscreteParametricProcess*[nb_output_process];
   // copy output processes and shift indices
   for(i = 0; i < nb_output_process; i++)
   {
      if (nptprocess[i+1] != NULL)
         npprocess[i] = new CategoricalProcess(*nptprocess[i+1]);
      else
         npprocess[i] = NULL;
      dpprocess[i] = dpprocess_copy[i+1];
   }
   for(i = 0; i < nb_output_process; i++)
      cprocess[i] = NULL;

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
      pgeneration[i] = new DiscreteMultivariateParametric(nb_state, MNEGATIVE_MULTINOMIAL,
                                                     inf_bound, sup_bound, parameter, probability);
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
   }
   delete [] pgeneration;

   mot = new MarkovOutTree(nb_state, nb_output_process, npprocess, dpprocess, cprocess,
                           unordered_rank, nb_vomc, pvomc_models, pgeneration_distributions);

   delete pgeneration_distributions;
   pgeneration_distributions = NULL;

   /* delete hmiot;
   hmiot = NULL;
   // pointers in observation processes have been deallocated
   // by "delete hmiot" (even if the pointers are not NULL)
   for(i = 0; i < nb_output_process; i++)
   {
      if (npprocess[i] != NULL)
         delete npprocess[i];
      npprocess[i] = NULL;
      // nptprocess[i] = NULL;
      dpprocess[i] = NULL;
      cprocess[i] = NULL;
   }
   // delete [] nptprocess;
   delete [] npprocess;
   delete [] dpprocess;
   delete [] cprocess; */
   // delete [] dpprocess_copy;

   cout << "Build MarkovOutTree from pointers: " << endl;
   mot->ascii_write(cout, false);
   mot->ascii_write(error, output_path, true);
   assert(error.get_nb_error() == 0);
   mot_read = markov_out_tree_ascii_read(error, mot_input_path);
   cout << error << endl;
   if (mot_read != NULL)
   {
      cout << "Read MarkovOutTree from a file: " << endl;
      mot_read->ascii_write(cout, false);
   }
   assert(mot_read != NULL);

   delete mot_read;
   mot_read = NULL;

   mot_read = markov_out_tree_ascii_read(error, mot_sim_path);
   cout << error << endl;
   if (mot_read != NULL)
   {
      cout << "Read MarkovOutTree from a file: " << endl;
      mot_read->ascii_write(cout, false);
   }
   assert(mot_read != NULL);

   // test get_plotable
   plot = mot_read->get_plotable();
   assert(plot != NULL);
   delete plot;
   plot = NULL;


   // simulate a Markov out tree
   cout << "Simulate Markov out-tree: " << endl;
   mtd = mot_read->simulation(error, depths, factor_indices, false, false);
   cout << error;
   assert(error.get_nb_error() == 0);

   // print the simulated trees

   cout << endl << "Simulated state trees : " << endl;
   for (t = 0; t < mtd->get_nb_trees(); t++)
      (mtd->get_state_tree_ptr(t))->display(cout, 0);
   cout << endl;

   cout << "Simulated trees : " << endl;
   for (t = 0; t < mtd->get_nb_trees(); t++)
   {
      cout << "tree " << t << ": " << endl;
      (mtd->get_tree_ptr(t))->display(cout, 0);
      cout << endl;
   }

   delete mtd;
   mtd = NULL;

   delete mot_read;
   mot_read = NULL;

   // multitype branching process
   mot_read = markov_out_tree_ascii_read(error, mbp_input_path);
   cout << error << endl;
   if (mot_read != NULL)
   {
      cout << "Read MultitypeBranchingProcess from a file: " << endl;
      mot_read->ascii_write(cout, false);
   }
   assert(mot_read != NULL);

   depths.resize(6);
   for(t=0; t < depths.size(); t++)
      depths[t] = 3;

   cout << "Simulate multitype branching process: " << endl;
   mtd = mot_read->simulation(error, depths, factor_indices, false, false);
   cout << error;
   assert(error.get_nb_error() == 0);

   // print the simulated trees

   cout << endl << "Simulated state trees : " << endl;
   for(t = 0; t < mtd->get_nb_trees(); t++)
   {
      cout << "tree " << t << ": " << endl;
      (mtd->get_state_tree_ptr(t))->display(cout, 0);
      cout << endl;
   }
   cout << endl;

   cout << "Generation probabilities: " << endl;
   factor_values = new index_array(1);
   generation_values.resize(3);
   generation_values[0] = 1;
   generation_values[1] = 1;
   generation_values[2] = 0;
   (*factor_values)[0] = 0;
   cout << "0 -> (1, 1): " << mot_read->get_generation_ptr(error, *factor_values)->get_mass(generation_values, true)
        << endl;
   cout << error;
   generation_values[0] = 0;
   generation_values[1] = 0;
   cout << "0 -> (0, 0): " << mot_read->get_generation_ptr(error, *factor_values)->get_mass(generation_values, true)
        << endl;
   cout << error;
   (*factor_values)[0] = 1;
   generation_values[0] = 1;
   generation_values[1] = 1;
   cout << "1 -> (1, 1): " << mot_read->get_generation_ptr(error, *factor_values)->get_mass(generation_values, true)
        << endl;
   cout << error;

   delete factor_values;
   factor_values = NULL;
   // compute partial likelihood
   cout << "Partial likelihood for tree 1: "
        << mot_read->likelihood_computation(*mtd, true, NULL, 1) << endl;

   // print censored vertices
   cout << endl << "Vertices with censored descendants in tree 1: " << endl;
   for(i = 0; i < (mtd->get_state_tree_ptr(1))->get_size(); i++)
      if (mtd->is_virtual(1, i))
         cout << i << " ;" ;
   cout << endl;

   delete mtd;
   mtd = NULL;

   // test estimation
   cout << "Estimate a multitype branching process. " << endl;
   depths.resize(6);
   for(t=0; t < depths.size(); t++)
      depths[t] = 7;

   // re-simulate multitype branching processes
   mtd = mot_read->simulation(error, depths, factor_indices, false, false);
   cout << error;
   assert(error.get_nb_error() == 0);

   cout << "Trees sizes: ";
   total_size = 0;
   for(i = 0; i < mtd->get_nb_trees(); i++)
   {
      cout << mtd->get_size(i) << "; " ;
      total_size += mtd->get_size(i);
   }
   cout << endl;

   if (total_size < 100)
   {
      cout << endl << "Simulated state trees : " << endl;
      for(t = 0; t < mtd->get_nb_trees(); t++)
      {
         cout << "tree " << t << ": " << endl;
         (mtd->get_state_tree_ptr(t))->display(cout, 0);
         cout << endl;
      }
      cout << endl;
   }

   embp = mtd->markov_out_tree_estimation(error, mot_read->get_nb_vomc(),
                                          mot_read->get_nb_ordered_children(),
                                          mot_read->get_nb_children_branching());

   cout << error;
   if (error.get_nb_error() == 0)
   {
      cout << "Estimated model: " << endl;
      embp->ascii_write(cout, false);
      cout << endl;
      delete embp;
      embp = NULL;
   }

   delete mot_read;
   mot_read = NULL;

   // multitype branching process with output processes
   mot_read = markov_out_tree_ascii_read(error, mbp_obs_input_path);
   cout << error << endl;
   if (mot_read != NULL)
   {
      cout << "Read MultitypeBranchingProcess from a file: " << endl;
      mot_read->ascii_write(cout, false);
   }
   assert(mot_read != NULL);

   // test estimation
   cout << endl << "Estimate a multitype branching process with output processes. " << endl;
   depths.resize(6);
   for(t=0; t < depths.size(); t++)
      depths[t] = 7;

   // re-simulate multitype branching processes
   mtd = mot_read->simulation(error, depths, factor_indices, false, false);
   cout << error;
   assert(error.get_nb_error() == 0);

   cout << "Trees sizes: ";
   total_size = 0;
   for(i = 0; i < mtd->get_nb_trees(); i++)
   {
      cout << mtd->get_size(i) << "; " ;
      total_size += mtd->get_size(i);
   }
   cout << endl;

   if (total_size < 100)
   {
      cout << endl << "Simulated state trees : " << endl;
      for(t = 0; t < mtd->get_nb_trees(); t++)
      {
         cout << "tree " << t << ": " << endl;
         (mtd->get_state_tree_ptr(t))->display(cout, 0);
         cout << endl;
      }
      cout << endl;
   }

   embp = mtd->markov_out_tree_estimation(error, mot_read->get_nb_vomc(),
                                          mot_read->get_nb_ordered_children(),
                                          mot_read->get_nb_children_branching());

   cout << error;
   if (error.get_nb_error() == 0)
   {
      cout << "Estimated model: " << endl;
      embp->ascii_write(cout, false);
      cout << endl;
      delete embp;
      embp = NULL;
   }

   // test get_plotable for MarkovOutTreeData
   /* plot = mtd->get_plotable();
   assert(plot != NULL);
   delete plot;
   plot = NULL; */

   // test get_plotable for MarkovOutTree
   /* plot = mot_read->get_plotable();
   assert(plot != NULL);
   delete plot;
   plot = NULL;*/

   mtd_cp = mtd->get_state_markov_out_tree_data();


   // reestimate model from state_tree
   cout << "Estimate a multitype branching process with output processes "
        <<  "from get_state_markov_out_tree_data()."  << endl;
   embp = mtd_cp->markov_out_tree_estimation(error, mot_read->get_nb_vomc(),
                                             mot_read->get_nb_ordered_children(),
                                             mot_read->get_nb_children_branching());

   cout << error;
   if (error.get_nb_error() == 0)
   {
      cout << "Estimated model: " << endl;
      embp->ascii_write(cout, false);
      cout << endl;
      delete embp;
      embp = NULL;
   }

   // build a MarkovOutTreeData from a trees by selecting some state variable
   cout << "build a MarkovOutTreeData from a Trees "
        << " by selecting some state variable" << endl;

   sel_mtd = new MarkovOutTreeData(*mtd_cp, 0);
   sel_mtd->ascii_write(cout, false);

   // test get_plotable for MarkovOutTreeData
   // without a model
   sel_mtd->update_markov_reestimation(error);
   cout << error;
   assert(error.get_nb_error() == 0);

   /* plot = sel_mtd->get_plotable();
   assert(plot != NULL);
   delete plot;
   plot = NULL;*/

   delete sel_mtd;
   sel_mtd = NULL;

   // build a MarkovOutTreeData from a trees by selecting a unique state variable
   cout << "build a MarkovOutTreeData from a Trees "
        << " by selecting the only variable" << endl;

   trees = new Trees(*mtd_cp);
   sel_var[0] = 1;
   trees_v1 = trees->select_variable(error, 1, sel_var);
   delete [] sel_var;
   cout << error;
   assert(error.get_nb_error() == 0);

   delete trees;
   trees = NULL;
   sel_mtd = new MarkovOutTreeData(*trees_v1, 0);
   sel_mtd->ascii_write(cout, true);

   delete trees_v1;
   trees_v1 = NULL;

   // test get_plotable for MarkovOutTreeData
   // without a model
   sel_mtd->update_markov_reestimation(error);
   cout << error;
   assert(error.get_nb_error() == 0);

   /* plot = sel_mtd->get_plotable();
   assert(plot != NULL);
   delete plot;
   plot = NULL;*/

   delete sel_mtd;
   sel_mtd = NULL;

   delete mtd_cp;
   mtd_cp = NULL;

   // test get_distribution_data
   // nb_state = mtd->get_nb_states();
   generation_values.resize(1);
   generation_values[0] = 1;
   distrib_data = mtd->get_distribution_data(error, generation_values);
   cout << error;
   assert(error.get_nb_error() == 0);
   for(i = 0; i < nb_state; i++)
   {
      if (distrib_data[i]->nb_element > 0)
      {
         cout << "Extract DistributionData " << i
              << " for parent " << generation_values[0] << ": " << endl;
         distrib_data[i]->ascii_write(cout, true);
         cout << endl;
      }
      delete distrib_data[i];
      distrib_data[i] = NULL;
   }

   // test copy constructor
   mtd_cp = new MarkovOutTreeData(*mtd, true, true);
   delete mtd;
   mtd = NULL;
   delete mtd_cp;
   mtd_cp = NULL;


   delete mot_read;
   mot_read = NULL;

   delete mot;
   mot = NULL;

   // include external factors
   nb_factors = 2;
   factor_values = new index_array(nb_factors);
   nb_generation_values_fact = 1;
   for(i = 0; i < factor_values->size(); i++)
   {
      (*factor_values)[i] = i + 2;
      nb_generation_values_fact *= (*factor_values)[i];
   }
   // number of factors including ordered brother and parent
   nb_generation = 1 + nb_factors + nb_vomc;
   // number of factors combinations
   nb_generation_values = nb_generation_values_fact * nb_state * pow(nb_state, nb_vomc);
   pgeneration = new DiscreteMultivariateDistribution*[nb_generation_values];

   for(i = 0; i < nb_generation_values; i++)
   {
      parameter[0] = i + 1;
      pgeneration[i] = new DiscreteMultivariateParametric(nb_state, MNEGATIVE_MULTINOMIAL,
                                                          inf_bound, sup_bound, parameter, probability);
   }



   pgeneration_distributions = new GenerationProcess(nb_state,
                                                     nb_children_branching,
                                                     nb_factors, factor_values,
                                                     nb_generation, pgeneration);
   assert(pgeneration_distributions->code(pgeneration_distributions->decode(nb_generation_values_fact-1))
          == nb_generation_values_fact - 1);

   for(i = 0; i < nb_generation_values_fact; i++)
   {
      delete pgeneration[i];
      pgeneration[i] = NULL;
   }
   delete [] pgeneration;

   mot = new MarkovOutTree(nb_state, nb_output_process, npprocess, dpprocess, cprocess,
                           unordered_rank, nb_vomc, pvomc_models, pgeneration_distributions);

   delete hmiot;
   hmiot = NULL;
   delete factor_values;
   factor_values = NULL;

   delete pgeneration_distributions;
   pgeneration_distributions = NULL;

   // pointers in observation processes have been deallocated
   // by "delete hmiot" (even if the pointers are not NULL)
   for(i = 0; i < nb_output_process; i++)
   {
      if (npprocess[i] != NULL)
         delete npprocess[i];
      npprocess[i] = NULL;
      // nptprocess[i] = NULL;
      dpprocess[i] = NULL;
      cprocess[i] = NULL;
   }
   // delete [] nptprocess;
   delete [] npprocess;
   delete [] dpprocess;
   delete [] cprocess;
   // delete [] dpprocess_copy;

   cout << "Build MarkovOutTree with factors from pointers: " << endl;
   mot->ascii_write(cout, false);
   mot->ascii_write(error, output_path, true);
   assert(error.get_nb_error() == 0);
   mot_read = markov_out_tree_ascii_read(error, mot_fact_path);

   delete mot;
   mot = NULL;

   cout << error << endl;
   if (mot_read != NULL)
   {
      cout << "Read MarkovOutTree from a file: " << endl;
      mot_read->ascii_write(cout, false);
   }
   assert(mot_read != NULL);

   // simulate a Markov out tree
   depths.resize(6);
   for(t=0; t < depths.size(); t++)
      depths[t] = 3;
   cout << "Simulate Markov out-tree." << endl;
   factor_indices.resize(2);
   factor_indices[0] = 0;
   factor_indices[1] = 1;
   mtd = mot_read->simulation(error, depths, factor_indices, false, false);
   cout << error;
   assert(error.get_nb_error() == 0);

   // print the simulated trees
   /*
   cout << endl << "Simulated state trees : " << endl;
   for (t = 0; t < mtd->get_nb_trees(); t++)
      (mtd->get_state_tree_ptr(t))->display(cout, 0);
   cout << endl;

   cout << "Simulated trees : " << endl;
   for (t = 0; t < mtd->get_nb_trees(); t++)
   {
      cout << "tree " << t << ": " << endl;
      (mtd->get_tree_ptr(t))->display(cout, 0);
      cout << endl;
   }
   */
   // combinations of factors
   cout << "Combinations of factors: " << endl;
   comb = mtd->get_factor_combinations();
   for (i = 0; i < comb.size(); i++)
      cout << comb[i] << endl;


   delete mot_read;
   mot_read = NULL;

   delete mtd;
   mtd = NULL;

   return 0;
}
