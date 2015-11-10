/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       Forum for OpenAlea developers: Openalea-devlp@lists.gforge.inria.f
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

#include "statiskit/core/data/marginal/multivariate.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution_reestimation.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"   // definition of DiscreteParametricModel class
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/sequence_label.h"

// required by "multivariate_distribution.h"
#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_label.h"
#include "tree_labels.h"
// must be included before multivariate_distribution.h
// since the latter is template and makes use of the former
#include "mv_histogram_tools.h"
#include "multivariate_distribution.h"

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"

#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"

#include "hidden_markov_tree.h"
#include "markov_out_tree.h"

#include <algorithm>

using namespace stat_tool;
using namespace sequence_analysis;
using namespace Stat_trees;

// for simulation; see stat_tool
extern int cumul_method(int nb_value, const double *cumul, double scale);

/*****************************************************************
 *
 *  Simulation of a MarkovOutTree using a StatError,
 *  a vector of tree depths,
 *  the indiced of simulated observed variables to be used as external factors
 *  a flag on the counting distribution computation
 *  and a flag on the Kullback-Leibler discrepancy computation,
 *
 *
 **/

MarkovOutTreeData* MarkovOutTree::simulation(StatError& error,
                                             const std::vector<unsigned int>& depths,
                                             const std::vector<unsigned int>& factors,
                                             bool counting_flag,
                                             bool divergence_flag) const
{
   typedef Trees::tree_type tree_type;
   typedef generic_visitor<tree_type>::vertex_array vertex_array;
   typedef GenerationProcess::index_array index_array;
   typedef MarkovOutTreeData::virtual_vdic virtual_vdic;

   const unsigned int nb_generation = get_nb_generation(),
                      nb_children_branching = get_nb_children_branching(),
                      nb_factors = factors.size();
   const bool parent_generation = (nb_generation == nb_children_branching + nb_factors + 1);
              // true iif children depend on parent state
   bool status = true, remove_cumul = false;
   register int t, i, var, ivar, dvar, fact;
   unsigned int node, cdepth, // current depth
                nb_depth_vertices, // nb vertices at last depth level
                nb_depth_vertices_new, // nb vertices at current depth level
                cnode, factor_config_index, nb_children, fact_val=0;
   const unsigned int nb_trees = depths.size();
   int parent_state;
       // cumul_depth could be used to limit the cumulated depth
   int nb_ioutput_process = 0, nb_doutput_process = 0;
   ostringstream error_message;
   std::deque<unsigned int> parent_state_list, parent_state_id_list;
   Typed_edge_one_int_tree::value state_v; // state variable for current vertex
   Typed_edge_one_int_tree::key v;
   tree_type::value obs_v;
   tree_type::key parent_key;
   MarkovOutTreeData *res = NULL;
   MarkovOutTree *markov;
   VariableOrderMarkovData *ordered_children = NULL;
   VariableOrderMarkov *ordered_children_model = NULL;
   // model used to simulate ordered children; does not need to be deallocated
   index_array factor_config, unordered_children, nb_unordered_children,
               factor_max_values;
   DiscreteMultivariateDistribution *unordered_children_model = NULL;
   // model used to simulate unordered children; does not need to be deallocated
   Typed_edge_one_int_tree *state_tree;
   tree_type *tree;
   virtual_vdic *vdic;

   error.init();

   if ((nb_trees < 1) || (nb_trees > NB_TREES))
   {
      status = false;
      error.update(STAT_TREES_error[TREESTATR_NB_TREES]);
   }
   if ((unordered_rank - 1 > 0) && (vomc_models == NULL))
   {
      status = false;
      error.update(STAT_TREES_parsing[TREESTATP_NB_VOMC]);
   }
   if (nb_factors > 0)
   {
      if (generation_distributions == NULL)
      {
         status = false;
         error.update(STAT_TREES_parsing[TREESTATP_NB_FACTORS]);
      }
      else
      {
         if (generation_distributions->get_nb_factors() != nb_factors)
         {
            status = false;
            error.update(STAT_TREES_parsing[TREESTATP_NB_FACTORS]);
         }
         factor_max_values = generation_distributions->get_factor_values();
      }
   }

   for(t = 0; t < nb_trees; t++)
   {
      // *std::min_element(depths.begin(), depths.end())
      if (depths[t] < 2)
      {
         status = false;
         error_message << STAT_TREES_label[TREESTATL_TREE] << " " << t + 1 << ": "
                       << STAT_TREES_error[TREESTATR_SMALL_TREE_SIZE];
         error.update((error_message.str()).c_str());
      }
      if (depths[t] > MAX_SIZE)
      {
         status = false;
         error_message << STAT_TREES_label[TREESTATL_TREE] << " " << t + 1 << ": "
                       << STAT_TREES_error[TREESTATR_BIG_TREE_SIZE];
      }
   }

   if (status)
   {
      for(var = 0; var < nb_output_process; var++)
      {
         if (nonparametric_process[var+1] != NULL)
            nb_ioutput_process++;
         else
            {
               if (discrete_parametric_process[var+1] != NULL)
                  nb_ioutput_process++;
               else
                  nb_doutput_process++;
            }
      }

      // initializations
      parent_state_id_list.resize(0);
      parent_state_list.resize(0);
      factor_config.resize(nb_generation);

      res = new MarkovOutTreeData(nb_ioutput_process, nb_doutput_process,
                                  nb_trees);

      for(var= 0; var < nb_ioutput_process; var++)
         res->_type[var] = INT_VALUE;

      for(var= 0; var < nb_doutput_process; var++)
         res->_type[var + nb_ioutput_process] = REAL_VALUE;

      res->markov = new MarkovOutTree(*this, false);

      markov = res->markov; // *markov is a MarkovOutTree
      markov->create_cumul();
      remove_cumul = true;
      markov->cumul_computation();
      // N.B. : Chain::create_cumul(), etc.

      obs_v.reset(nb_ioutput_process, nb_doutput_process);

      if (nb_vomc == 1)
         ordered_children_model = vomc_models[0];

      for(t = 0; t < res->get_nb_trees(); t++)
      {
         parent_state_list.clear();
         assert(parent_state_id_list.empty());
         state_tree = res->state_trees[t];
         tree = res->trees[t];
         vdic = res->virt_vertices[t];

         // simulation of root state
         state_v.Int() = cumul_method(markov->nb_state, markov->cumul_initial);
         node = state_tree->root();
         state_tree->put(node, state_v);

         // simulation of the root output variables
         ivar = 0;
         dvar = 0;
         for(var = 0; var < markov->nb_output_process; var++)
            if (markov->nonparametric_process[var+1] != NULL)
               obs_v.Int(ivar++)= markov->nonparametric_process[var+1]->observation[state_v.Int()]->simulation();
           else
              {
                 if (markov->discrete_parametric_process[var+1] != NULL)
                    obs_v.Int(ivar++)= markov->discrete_parametric_process[var+1]->observation[state_v.Int()]->simulation();
                 else
                    obs_v.Double(dvar++)= markov->continuous_parametric_process[var+1]->observation[state_v.Int()]->simulation();
              }

         if (nb_output_process > 0)
         {
            assert(node == tree->root());
            tree->put(node, obs_v);
         }

         parent_state_id_list.push_back(state_tree->root());
         parent_state_list.push_back(state_v.Int());
         nb_depth_vertices_new = 1;

         for(cdepth = 0; cdepth < depths[t]-1; cdepth++)
         {
            // simulate leaf vertices
            nb_depth_vertices = nb_depth_vertices_new;
            nb_depth_vertices_new = 0;
            for(cnode = 0; cnode < nb_depth_vertices; cnode++)
            {
               if (status)
               {
                  // add children to current vertex

                  // get next vertex in the list of parents
                  node = parent_state_id_list.front();
                  parent_state_id_list.pop_front();
                  parent_state = parent_state_list.front();
                  parent_state_list.pop_front();
                  // current type of factor for unordered children
                  // (parent, ordered brother vertex states, external factors)
                  factor_config_index = 0;
                  // simulate children states:
                  // 1) add ordered vertices

                  // select VariableOrderMarkov model and initial probabilities
                  if (nb_vomc == nb_state)
                     ordered_children_model = vomc_models[parent_state];

                  if (nb_vomc == 1)
                     for(i = 0; i < nb_state; i++)
                        ordered_children_model->initial[i] = double(parent_state == i);

                  factor_config_index = 0;
                  if (parent_generation)
                     factor_config[factor_config_index++] = parent_state;

                  if (unordered_rank - 1 > 0)
                  {
                     ordered_children = ordered_children_model->simulation(error, 1, max((unsigned int)2, unordered_rank - 1), false);
                     if (error.get_nb_error() > 0)
                        cerr << error;
                     assert(error.get_nb_error() == 0);
                     for(i = 0; i < unordered_rank - 1; i++)
                     {
                        state_v.Int() = ordered_children->get_int_sequence(0, 0, i);
                        v = state_tree->add_vertex(state_v);
                        state_tree->add_edge(node, v);
                        parent_state_id_list.push_back(v);
                        parent_state_list.push_back(state_v.Int());

                        // simulation of the output variables
                        ivar = 0;
                        dvar = 0;
                        for(var = 0; var < markov->nb_output_process; var++)
                           if (markov->nonparametric_process[var+1] != NULL)
                              obs_v.Int(ivar++)= markov->nonparametric_process[var+1]->observation[state_v.Int()]->simulation();
                          else
                             {
                                if (markov->discrete_parametric_process[var+1] != NULL)
                                   obs_v.Int(ivar++)= markov->discrete_parametric_process[var+1]->observation[state_v.Int()]->simulation();
                                else
                                   obs_v.Double(dvar++)= markov->continuous_parametric_process[var+1]->observation[state_v.Int()]->simulation();
                             }

                        if (nb_output_process > 0)
                        {
                           assert(v == tree->add_vertex(obs_v));
                           tree->add_edge(node, v, false);
                        }

                        nb_depth_vertices_new++;
                        if (i < nb_children_branching)
                           factor_config[factor_config_index++] = state_v.Int();
                     }
                     delete ordered_children;
                     ordered_children = NULL;
                     // add type < to 1st child
                     state_tree->set_edge_type(node, v, true);
                     if (nb_output_process > 0)
                        tree->set_edge_type(node, v, true);
                  }

                  // 2) add unordered vertices

                  if (generation_distributions != NULL)
                  {
                     // external factors (if any)
                     if (nb_factors > 0)
                        for(fact = 0; fact < nb_factors; fact++)
                        {
                           fact_val = obs_v.Int(factors[fact]);
                           if (fact_val >= factor_max_values[fact])
                           {
                              status = false;
                              error_message << STAT_TREES_label[TREESTATL_TREE] << " " << t + 1 << ": "
                                            << STAT_TREES_parsing[TREESTATP_NB_FACTOR_VALUES] << ": "
                                            << fact_val;
                              error.update((error_message.str()).c_str());
                              fact_val = 0;

                           }
                           factor_config[factor_config_index++] = fact_val;
                        }
                     unordered_children_model = generation_distributions->generation[generation_distributions->code(factor_config)];

                     // number of unordered vertices of each type
                     nb_unordered_children = unordered_children_model->simulation();

                     unordered_children.resize(0);
                     // unordered_children should be permuted at random
                     // (permutation with uniform distribution)
                     nb_children = 0;
                     for(var = 0; var < nb_state; var++)
                     {
                        for(i = 0; i < nb_unordered_children[var]; i++)
                           unordered_children.push_back(var);
                        nb_children += nb_unordered_children[var];
                     }
                     if ((nb_children == 0) && (unordered_rank < 2))
                        // current vertex actually has no children
                        vdic->insert(pair<int, bool>(node, false));

                     std::random_shuffle(unordered_children.begin(), unordered_children.end());

                     for(i = 0; i < nb_children; i++)
                     {
                        state_v.Int() = unordered_children[i];
                        v = state_tree->add_vertex(state_v);
                        state_tree->add_edge(node, v);
                        parent_state_id_list.push_back(v);
                        parent_state_list.push_back(state_v.Int());

                        // simulation of the output variables
                        ivar = 0;
                        dvar = 0;
                        for(var = 0; var < markov->nb_output_process; var++)
                           if (markov->nonparametric_process[var+1] != NULL)
                              obs_v.Int(ivar++)= markov->nonparametric_process[var+1]->observation[state_v.Int()]->simulation();
                          else
                             {
                                if (markov->discrete_parametric_process[var+1] != NULL)
                                   obs_v.Int(ivar++)= markov->discrete_parametric_process[var+1]->observation[state_v.Int()]->simulation();
                                else
                                   obs_v.Double(dvar++)= markov->continuous_parametric_process[var+1]->observation[state_v.Int()]->simulation();
                             }

                        if (nb_output_process > 0)
                        {
                           assert(v == tree->add_vertex(obs_v));
                           tree->add_edge(node, v, false);
                        }
                     }
                     nb_depth_vertices_new += nb_children;
                  }
               } // end if status
            } // end previous list of leaf vertices
         } // end for: current generation

         // add virtual vertices (leaf vertices at last step)
         if (status)
         {
            assert(parent_state_id_list.size() == nb_depth_vertices_new);
            nb_depth_vertices = nb_depth_vertices_new;
            for(cnode = 0; cnode < nb_depth_vertices; cnode++)
            {
               // get next vertex in the list of parents
               node = parent_state_id_list.front();
               parent_state_id_list.pop_front();
               vdic->insert(pair<int, bool>(node, true));
            }
         }
         else
         {
            parent_state_id_list.clear();
            parent_state_list.clear();
         }
      }

      if (status)
      {
         // extraction of the characteristics for the simulated trees

         res->min_max_value_computation();

         res->chain_data = new ChainData(type, nb_state, nb_state);

         if (nb_ioutput_process + nb_doutput_process > 0);
            res->build_characteristics();
         res->build_size_frequency_distribution();
         res->build_nb_children_frequency_distribution();
         res->_nb_states = nb_state;
         if (nb_ioutput_process + nb_doutput_process > 0);
            res->build_observation_frequency_distribution();
         // res->build_state_characteristics(); // called by build_observation_frequency_distribution();
         res->build_state_frequency_distribution(error, nb_children_branching,
                                                 nb_factors, &factor_max_values,
                                                 nb_generation, unordered_rank,
                                                 nb_vomc, &factors);
         if (error.get_nb_error() > 0)
         {
            status = false;
            delete res;
            res = NULL;
         }
         if (status && (!divergence_flag))
         {
            markov->characteristic_computation(*res, counting_flag);
            markov->log_computation();

            // likelihood computation
            res->likelihood = markov->likelihood_computation(*res, &factors);

         }
      }
      else
      {
         delete res;
         res = NULL;
      }
      if ((status) && (remove_cumul))
         markov->remove_cumul();

   }
   return res;
}

/*****************************************************************
 *
 *  Compute the log-likelihood of a MarkovOutTree for a given set
 *  of observed trees, a tree identifier, the indices of factor variables
 *  and the index of the state variable
 *
 **/

double MarkovOutTree::likelihood_computation(const Trees& trees, unsigned int state_variable,
                                             const std::vector<unsigned int>* factor_indices,
                                             int index) const
{
   double loglikelihood = D_INF;

   return loglikelihood;
}

/*****************************************************************
 *
 *  Compute the log-likelihood of a MarkovOutTree
 *  for a given MarkovOutTreeData instance
 *
 **/

double MarkovOutTree::likelihood_computation(const MarkovOutTreeData& trees,
                                             const std::vector<unsigned int>* factor_indices,
                                             int index) const
{
   return likelihood_computation(trees, false, factor_indices, I_DEFAULT);
}


/*****************************************************************
 *
 *  Compute the log-likelihood of a MarkovOutTree
 *  for a given MarkovOutTreeData instance,
 *  the indices of factor variables and a tree identifier,
 *  accounting or not for possible permutations of subtrees
 *
 **/

double MarkovOutTree::likelihood_computation(const MarkovOutTreeData& trees, bool partial,
                                             const std::vector<unsigned int>* factor_indices,
                                             int index) const
{
   bool reestimation_computation = false;
   unsigned int t, inb_children_branching = 0, inb_factors = 0, inb_generation = 0;
   double loglikelihood = 0, inc_loglikelihood = D_INF;
   int *iindividual = NULL;
   StatError error;
   MarkovOutTreeData *cp_trees = NULL;
   MarkovOutTreeData::index_array *ifactor_values = NULL;
   const MarkovOutTreeReestimation< int > *reestimation = trees.mot_reestimation;

   // determine whether the multivariate histograms
   // and sequences of ordered children need to be computed
   if (generation_distributions != NULL)
   {
      inb_children_branching = generation_distributions->nb_children_branching;
      inb_factors = generation_distributions->nb_factors;
      inb_generation = generation_distributions->nb_generation;
      ifactor_values = generation_distributions->factor_values;
   }
   if (reestimation == NULL)
      reestimation_computation = true;
   else
   {
      if ((reestimation->get_nb_state() != nb_state) ||
          (reestimation->get_nb_children_branching() != inb_children_branching) ||
          (reestimation->get_nb_factors() != inb_factors) ||
          (reestimation->get_nb_generation() != inb_generation) ||
          (reestimation->get_unordered_rank() != unordered_rank) ||
          (reestimation->get_nb_vomc() != nb_vomc))
         reestimation_computation = true;
   }
   if ((index != I_DEFAULT) && (trees.get_nb_trees() != 1))
      // compute likelihood on a subsample
      reestimation_computation = true;

   if (reestimation_computation)
   {
      if ((index != I_DEFAULT) && (trees.get_nb_trees() != 1))
      // compute likelihood on a subsample
      {
         iindividual = new int[1];
         iindividual[0] = index;
         cp_trees = trees.select_individual(error, 1, iindividual, true);
         if (cp_trees == NULL)
            return D_INF;
      }
      else
         cp_trees = new MarkovOutTreeData(trees, true, true);

      // cp_trees->build_size_frequency_distribution();
      cp_trees->min_max_value_computation();
      // cp_trees->build_nb_children_frequency_distribution();
      cp_trees->build_chain_reestimation(type, nb_state);
      cp_trees->build_markov_reestimation(inb_children_branching, inb_factors, ifactor_values,
                                          inb_generation, unordered_rank, nb_vomc, factor_indices);
   }

   // compute part of loglikelihood related to initial probabilities
   if (cp_trees == NULL)
      inc_loglikelihood = initial_state_likelihood_computation(trees);
   else
      inc_loglikelihood = initial_state_likelihood_computation(*cp_trees);
   if (inc_loglikelihood <= D_INF)
      return D_INF;
   else
      loglikelihood += inc_loglikelihood;

   // compute part of loglikelihood related to output processes

   if (nb_output_process > 0)
   {
      if (cp_trees == NULL)
         inc_loglikelihood = output_likelihood_computation(trees);
      else
         inc_loglikelihood = output_likelihood_computation(*cp_trees);

      if (inc_loglikelihood <= D_INF)
         return D_INF;
      else
         loglikelihood += inc_loglikelihood;
   }

   if (nb_vomc > 0)
   {
      // compute part of loglikelihood related to ordered children
      if (cp_trees == NULL)
         inc_loglikelihood = ordered_model_likelihood_computation(trees);
      else
         inc_loglikelihood = ordered_model_likelihood_computation(*cp_trees);

      if (inc_loglikelihood <= D_INF)
         return D_INF;
      else
         loglikelihood += inc_loglikelihood;
   }

   // compute part of loglikelihood related to unordered children
   if (generation_distributions != NULL)
      if (partial)
      {
         if (cp_trees == NULL)
            loglikelihood += unordered_model_partial_likelihood_computation(trees);
         else
            loglikelihood += unordered_model_partial_likelihood_computation(*cp_trees);

      }
      else
      {
         // should be:
         // loglikelihood += unordered_model_likelihood_computation(*cp_trees);
         if (cp_trees == NULL)
            loglikelihood = unordered_model_likelihood_computation(trees);
         else
            loglikelihood = unordered_model_likelihood_computation(*cp_trees);
      }

   if (reestimation_computation)
   {
      delete cp_trees;
      cp_trees = NULL;
   }
   if (iindividual != NULL)
   {
      delete [] iindividual;
      iindividual = NULL;
   }

   return loglikelihood;
}

/* Simulation of trees using an empirical distribution for
   tree depths */

MarkovOutTreeData* MarkovOutTree::simulation(StatError& error, const FrequencyDistribution& idepth,
                                             bool depth_flag,
                                             bool counting_flag,
                                             bool divergence_flag) const
{
   MarkovOutTreeData *res = NULL;

   return res;
}

/* Simulation of trees using maximal values for tree depth or size */
MarkovOutTreeData* MarkovOutTree::simulation(StatError& error, int inb_trees,
                                             int imax_depth, int imax_size,
                                             bool depth_flag,
                                             bool counting_flag) const
{
   MarkovOutTreeData *res = NULL;

   return res;
}

/*****************************************************************
 *
 *  Compute output conditional distributions for MarkovOutTree class
 *  using a MarkovOutTreeData object,
 *  the stored conditional probabilities,
 *  a flag on the computation of the log probabilities
 *  and the index of the considered tree
 *
 **/

void MarkovOutTree::output_conditional_distribution(const MarkovOutTreeData& trees,
                                                    double_array_3d& output_cond,
                                                    bool log_computation,
                                                    int index) const
{
   typedef MarkovOutTreeData::tree_type tree_type;
   typedef tree_type::vertex_iterator vertex_iterator;
   typedef tree_type::value value;

   int nb_trees = trees._nb_trees, t, current_size;
   register int j, var;
   double proba;
   Typed_edge_int_fl_tree<Int_fl_container> *current_tree;
   vertex_iterator it, end;
   value val;

   // output_cond[t][j][u] corresponds to conditional distribution, given the state variable
   // is equal to j, taken at the value of node u of tree t

   assert(nb_output_process == trees._nb_integral + trees._nb_float);

   if (output_cond == NULL)
   {
      output_cond = new double_array_2d[nb_trees];
      for(t = 0; t < nb_trees; t++)
      output_cond[t] = NULL;
   }

   for(t = 0; t < nb_trees; t++)
      if ((index == I_DEFAULT) || (index == t))
      {
         if (output_cond[t] == NULL)
         {
            output_cond[t] = new double*[nb_state];
            for(j = 0; j < nb_state; j++)
            output_cond[t][j] = NULL;
         }

         current_tree = trees.trees[t];
         current_size = current_tree->get_size();
         for(j = 0; j < nb_state; j++)
            if (output_cond[t][j] == NULL)
               output_cond[t][j]= new double[current_size];

         Tree_tie::tie(it, end) = current_tree->vertices();
         while (it < end)
         {
            val = current_tree->get(*it);
            for(j= 0; j < nb_state; j++)
            {
               if (log_computation)
                  output_cond[t][j][*it] = 0.;
               else
                  output_cond[t][j][*it] = 1.;
               for(var = 0; var < nb_output_process; var++)
               {
                  if (nonparametric_process[var+1] != NULL)
                  {
                     if (log_computation)
                     {
#                       ifdef DEBUG
                        if (((nonparametric_process[var+1]->observation[j]->mass[val.Int(var)] > 0)
                             && (abs(log(nonparametric_process[var+1]->observation[j]->mass[val.Int(var)])-nonparametric_process[var+1]->observation[j]->cumul[val.Int(var)]) > DOUBLE_ERROR))
                           || ((nonparametric_process[var+1]->observation[j]->mass[val.Int(var)] == 0) && (abs(nonparametric_process[var+1]->observation[j]->cumul[val.Int(var)]-D_INF) > DOUBLE_ERROR)))
                           cout << "Warning: computation error at observation[" << var+1 << "]["
                                << j << "]." << endl;
#                       endif
                        if (nonparametric_process[var+1]->observation[j]->mass[val.Int(var)] != 0)
                           output_cond[t][j][*it] += log(nonparametric_process[var+1]->observation[j]->mass[val.Int(var)]);
                        else
                           output_cond[t][j][*it] = D_INF;
                     }
                     else
                        output_cond[t][j][*it] *= nonparametric_process[var+1]->observation[j]->mass[val.Int(var)];
                  }
                  if (discrete_parametric_process[var+1] != NULL)
                  {
                     if (log_computation)
                     {
#                       ifdef DEBUG
                        if (((discrete_parametric_process[var+1]->observation[j]->mass[val.Int(var)] > 0)
                             && (abs(log(discrete_parametric_process[var+1]->observation[j]->mass[val.Int(var)])-discrete_parametric_process[var+1]->observation[j]->cumul[val.Int(var)]) > DOUBLE_ERROR))
                           || ((discrete_parametric_process[var+1]->observation[j]->mass[val.Int(var)] == 0) && (abs(discrete_parametric_process[var+1]->observation[j]->cumul[val.Int(var)]-D_INF) > DOUBLE_ERROR)))
                           cout << "Warning: computation error at transition[" << var+1 << "]["
                                << j << "]." << endl;
#                       endif
                           if (discrete_parametric_process[var+1]->observation[j]->mass[val.Int(var)] != 0)
                              output_cond[t][j][*it] += log(discrete_parametric_process[var+1]->observation[j]->mass[val.Int(var)]);
                           else
                              output_cond[t][j][*it] = D_INF;
                     }
                     else
                        output_cond[t][j][*it] *= discrete_parametric_process[var+1]->observation[j]->mass[val.Int(var)];
                  }
                  if (continuous_parametric_process[var+1] != NULL)
                  {
                     assert(false);
                     //proba = continuous_parametric_process[var+1]->observation[j]->mass_computation(val.Double(var) - seq.min_interval[var+1] / 2 , val.Double(var) + seq.min_interval[var+1] / 2);
                     if (log_computation)
                     {
                        if (proba > 0)
                           output_cond[t][j][*it] += log(proba);
                        else
                           output_cond[t][j][*it] = D_INF;
                     }
                     else
                        output_cond[t][j][*it] *= proba;
                  }
               // case of floating observed processes, not implemented for the moment
               }
            }
            it++;
         }
   }
}

/*****************************************************************
 *
 *  Compute part of the loglikelihood related to ordered children
 *  in a MarkovOutTree, for a given MarkovOutTreeData instance.
 *  It is assumed that this instance contains updated mot_reestimation
 *
 **/

double MarkovOutTree::ordered_model_likelihood_computation(const MarkovOutTreeData& trees) const
{
   double loglikelihood = 0, inc_loglikelihood = 0;
   unsigned int dist;

   // this method should be updated to take into account potential missing values
   return loglikelihood;
}

/*****************************************************************
 *
 *  Compute part of the loglikelihood related to initial states
 *  in a MarkovOutTree, for a given MarkovOutTreeData instance.
 *  It is assumed that this instance contains updated chain_data
 *
 **/

double MarkovOutTree::initial_state_likelihood_computation(const MarkovOutTreeData& trees) const
{
   double loglikelihood = 0;
   register int j;
   ChainData *chain_data = trees.chain_data;

   if (chain_data != NULL)
   {
      for(j = 0; j < nb_state; j++)
      {
         if (chain_data->initial[j] > 0)
            if (initial[j] > 0)
               loglikelihood += chain_data->initial[j] * log(initial[j]);
            else
               return D_INF;
      }
   }
   else
      return D_INF;
   return loglikelihood;
}

/*****************************************************************
 *
 *  Compute part of the loglikelihood related to output processes
 *  in a MarkovOutTree for a given MarkovOutTreeData instance.
 *
 **/

double MarkovOutTree::output_likelihood_computation(const MarkovOutTreeData& trees) const
{
   typedef MarkovOutTreeData::tree_type tree_type;
   typedef MarkovOutTreeData::state_tree_type state_tree_type;
   typedef MarkovOutTreeData::state_value value;
   typedef MarkovOutTreeData::vertex_iterator vertex_iterator;

   unsigned int t, j, state;
   double loglikelihood = 0, inc_loglikelihood = 0;
   double_array_3d output_cond = NULL;
   vertex_iterator it, end;
   tree_type *current_tree;
   state_tree_type *current_state_tree;

   output_conditional_distribution(trees, output_cond, true);

   for(t = 0; t < trees.get_nb_trees(); t++)
   {
      current_tree = trees.trees[t];
      current_state_tree = trees.state_trees[t];
      Tree_tie::tie(it, end) = current_tree->vertices();
      while (it < end)
      {
          state = (current_state_tree->get(*it)).Int();
          inc_loglikelihood = output_cond[t][state][*it];
          if (inc_loglikelihood <= D_INF)
             return D_INF;
          else
             loglikelihood += inc_loglikelihood;
          it++;
      }
   }

   for(t = 0; t < trees.get_nb_trees(); t++)
   {
      for(j = 0; j < nb_state; j++)
      {
         delete [] output_cond[t][j];
         output_cond[t][j] = NULL;
      }

      delete [] output_cond[t];
      output_cond[t] = NULL;
   }

   delete [] output_cond;
   output_cond = NULL;

   return loglikelihood;
}

/*****************************************************************
 *
 *  Compute part of the loglikelihood related to
 *  unordered children in a MarkovOutTree
 *  for a given MarkovOutTreeData instance.
 *  Possible permutations of subtrees are not accounted for.
 *  It is assumed that this instance contains updated mot_reestimation
 *
 **/

double MarkovOutTree::unordered_model_partial_likelihood_computation(const MarkovOutTreeData& trees) const
{
   double loglikelihood = 0, inc_loglikelihood = 0;
   unsigned int dist;
   MarkovOutTreeData::mvHistogram * const *generation_reestim = trees.mot_reestimation->get_generation_reestim_ptr();

   // this method should be updated to take into account potential missing values
   if ((generation_reestim != NULL) && (generation_distributions != NULL)
       && (generation_distributions->generation != NULL))
   {
      for(dist = 0; dist < generation_distributions->nb_generation_values; dist++)
      {
         inc_loglikelihood = generation_distributions->generation[dist]->likelihood_computation(*generation_reestim[dist]);
         if (inc_loglikelihood != D_INF)
            loglikelihood += inc_loglikelihood;
         else
            return D_INF;

      }
   }
   return loglikelihood;
}


/*****************************************************************
 *
 *  Compute part of the loglikelihood related to
 *  unordered children in a MarkovOutTree
 *  for a given MarkovOutTreeData instance.
 *  Possible permutations of subtrees are accounted for.
 *  It is assumed that this instance contains updated mot_reestimation
 *
 **/

double MarkovOutTree::unordered_model_likelihood_computation(const MarkovOutTreeData& trees) const
{ return D_INF; }

/*****************************************************************
 *
 *  Estimate the variable order Markov models of a MarkovOutTree
 *  from a MarkovOutTreeData instance
 *
 **/

void MarkovOutTree::vomc_estimation(StatError& error, const MarkovOutTreeData& trees)
{ }

/*****************************************************************
 *
 *  Estimate the generation distributions of a MarkovOutTree
 *  from a MarkovOutTreeData instance, the vector of distribution families,
 *  minimal possible inferior bound and a flag on estimating this bound
 *
 **/

void MarkovOutTree::generation_process_estimation(StatError& error,
                                                  const MarkovOutTreeData& trees,
                                                  const std::vector<int>& generation_types,
                                                  unsigned int generation_min_inf_bound,
                                                  bool generation_min_inf_bound_flag)
{
   // prerequisite: trees.mot_reestimation != NULL and is compatible with *this
   // if generation_types is empty, the family must be estimated

   const bool estimate_types = (generation_types.size() == 0);
   const unsigned int inb_children_branching = trees.mot_reestimation->get_nb_children_branching(),
                      inb_factors = trees.mot_reestimation->get_nb_factors(),
                      inb_generation = trees.mot_reestimation->get_nb_generation(),
                      inb_generation_values = trees.mot_reestimation->get_nb_generation_values();
   unsigned int g;
   double likelihood;
   MarkovOutTreeData::index_array ifactor_values = trees.mot_reestimation->get_factor_values();
   ostringstream error_message;
   DiscreteMultivariateDistribution **pgeneration = new DiscreteMultivariateDistribution*[inb_generation_values];
   MarkovOutTreeData::mvHistogram *generation_reestim = NULL;
   DiscreteMultivariateReestimation<int>* dmr = NULL;

   error.init();

   if (generation_distributions != NULL)
   {
      delete generation_distributions;
      generation_distributions = NULL;
   }

   if (estimate_types)
      for(g = 0; g < inb_generation_values; g++)
      {
         generation_reestim = trees.mot_reestimation->get_generation_reestim_ptr()[g];
         if (generation_reestim->get_total() <= nb_state)
         {
            error_message << "To few elements (" << generation_reestim->get_total()
                          << ") for states and factors configuration " << trees.mot_reestimation->decode(g);
#           ifdef MESSAGE
            print_stat_histogram(*generation_reestim, cout);
#           endif
            pgeneration[g] = NULL;
         }
         else
         {
            dmr = stat_histogram_to_DiscreteMultivariateReestimation(*generation_reestim);
            pgeneration[g] = dmr->type_parametric_estimation(error,
                                                             likelihood,
                                                             generation_min_inf_bound,
                                                             generation_min_inf_bound_flag);
            delete dmr;
            dmr = NULL;
            if ((pgeneration[g] == NULL) || (likelihood == D_INF))
            {
               error_message << " - at states and factors configuration " << trees.mot_reestimation->decode(g);
#              ifdef MESSAGE
               print_stat_histogram(*generation_reestim, cout);
#              endif
            }
         }
      }
   else
      for(g = 0; g < inb_generation_values; g++)
      {
         generation_reestim = trees.mot_reestimation->get_generation_reestim_ptr()[g];
         if (generation_reestim->get_total() <= nb_state)
         {
            error_message << "To few elements (" << generation_reestim->get_total()
                          << ") for states and factors configuration " << trees.mot_reestimation->decode(g);
#           ifdef MESSAGE
            print_stat_histogram(*generation_reestim, cout);
#           endif
            pgeneration[g] = NULL;
         }
         else
         {
            dmr = stat_histogram_to_DiscreteMultivariateReestimation(*generation_reestim);
            pgeneration[g] = dmr->parametric_estimation(error,
                                                        generation_types[g],
                                                        likelihood,
                                                        generation_min_inf_bound,
                                                        generation_min_inf_bound_flag);
            delete dmr;
            dmr = NULL;
            if ((pgeneration[g] == NULL) || (likelihood == D_INF))
            {
               error_message << " - at states and factors configuration " << trees.mot_reestimation->decode(g);
#              ifdef MESSAGE
               print_stat_histogram(*generation_reestim, cout);print_stat_histogram(*generation_reestim, cout);
#              endif
            }
         }
      }

   if (error.get_nb_error() == 0)
      generation_distributions = new GenerationProcess(nb_state, inb_children_branching,
                                                       inb_factors, &ifactor_values,
                                                       inb_generation, pgeneration);
   else
      error.update((error_message.str()).c_str());

   for(g = 0; g < inb_generation_values; g++)
      if (pgeneration[g] != NULL)
      {
         delete pgeneration[g];
         pgeneration[g] = NULL;
      }
   delete [] pgeneration;
   pgeneration = NULL;
}

/*****************************************************************
 *
 *  Estimate the output processes of a MarkovOutTree
 *  from a MarkovOutTreeData instance using a StatError object,
 *  the vector of output variables to ignore
 *  and a flag on common dispersion parameters (continuous output processes only)
 *
 **/

void MarkovOutTree::output_process_estimation(StatError &error,
                                              const MarkovOutTreeData& trees,
                                              const std::vector<bool> &skipped_variables,
                                              bool common_dispersion)
{
   bool status = true;
   register int var, j, val;
   const int nb_variable = trees.get_nb_int() + trees.get_nb_float();
   int type, max_nb_value = 0;
   int *cp_frequency = NULL;
   int nb_value[nb_output_process];
   double observation_likelihood, min_likelihood, buff;
   double *variance = NULL, **mean_direction = NULL;
   ostringstream error_message, correction_message;
   FrequencyDistribution *hobservation = NULL;
   FrequencyDistribution ** const * observation_distribution = trees.observation_distribution;

   if (nb_output_process != nb_variable)
   {
      status = false;
      error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
   }
   else
   {
      // check variable types
      for(var = 0; var < nb_variable; var++)
      {
         if (!skipped_variables[var])
         {
            type = trees.get_type(var);
            if ((type != INT_VALUE) && (type != REAL_VALUE)
                 && (type != STATE))
            {
               status = false;
               error_message << STAT_label[STATL_VARIABLE] << " " << var + 1 << ": "
                             << STAT_error[STATR_VARIABLE_TYPE];
               correction_message << STAT_variable_word[INT_VALUE] << " or "
                                  << STAT_variable_word[REAL_VALUE];
               error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
            }
         }
      }

      if (status)
      {
         for(var = 0; var < nb_variable; var++)
            if (!skipped_variables[var])
               {

               if ((trees.characteristics[var] != NULL)
                   && (trees.characteristics[var]->marginal_distribution != NULL))
                  nb_value[var] = trees.characteristics[var]->marginal_distribution->nb_value;
               else
                  nb_value[var] = I_DEFAULT;
            }
         // initialization and allocation of output processes
         output_process_init(nb_value);

         // compute maximal observed value for any integral variable
         // modelled by a discrete parametric process
         max_nb_value = 0;
         for(var = 0; var < nb_output_process; var++)
            if ((discrete_parametric_process[var+1] != NULL) &&
               (max_nb_value < trees.characteristics[var]->marginal_distribution->nb_value))
                max_nb_value = trees.characteristics[var]->marginal_distribution->nb_value;

         if (max_nb_value > 0)
           hobservation = new FrequencyDistribution(max_nb_value);

      }
   }

   if (status)
      for(var = 0; var < nb_output_process; var++)
         if (!skipped_variables[var])
         {

         if (nonparametric_process[var+1] != NULL)
         {
            if (observation_distribution[var] != NULL)
               for(j = 0; j < nb_state; j++)
               {
                  // make Distribution and FrequencyDistribution coincide
                  if (observation_distribution[var][j]->nb_value < nonparametric_process[var+1]->observation[j]->nb_value)
                  {
                     cp_frequency = observation_distribution[var][j]->frequency;
                     observation_distribution[var][j]->frequency = new int[nonparametric_process[var+1]->observation[j]->nb_value];
                     for(val = 0; val < observation_distribution[var][j]->nb_value; val++)
                        observation_distribution[var][j]->frequency[val] = cp_frequency[val];
                     delete [] cp_frequency;
                     cp_frequency = NULL;
                     observation_distribution[var][j]->nb_value = trees.characteristics[var]->marginal_distribution->nb_value;
                     for(val = observation_distribution[var][j]->nb_value; val < nonparametric_process[var+1]->observation[j]->nb_value; val++)
                        observation_distribution[var][j]->frequency[val] = 0;

                     observation_distribution[var][j]->nb_value = nonparametric_process[var+1]->observation[j]->nb_value;
                     observation_distribution[var][j]->nb_value_computation();
                     observation_distribution[var][j]->offset_computation();
                     observation_distribution[var][j]->nb_element_computation();
                     observation_distribution[var][j]->max_computation();
                     observation_distribution[var][j]->mean_computation();
                     observation_distribution[var][j]->variance_computation(true);
                  }
                  reestimation(trees.characteristics[var]->marginal_distribution->nb_value,
                               observation_distribution[var][j]->frequency ,
                               nonparametric_process[var+1]->observation[j]->mass,
                               MIN_PROBABILITY, false);
               }
         }
         else // (nonparametric_process[var+1] == NULL)
         {
            if (observation_distribution[var] != NULL)
            {
               for(j = 0; j < nb_state; j++)
               {
                  observation_distribution[var][j]->nb_value_computation();
                  observation_distribution[var][j]->offset_computation();
                  observation_distribution[var][j]->nb_element_computation();
                  observation_distribution[var][j]->max_computation();
                  observation_distribution[var][j]->mean_computation();
                  observation_distribution[var][j]->variance_computation(true);
               }
               if (discrete_parametric_process[var+1] != NULL)
               {
                  for(j = 0; j < nb_state; j++)
                  {
                     observation_likelihood =
                        // observation_distribution[var][j]->distribution_estimation(discrete_parametric_process[var+1]->observation[j]);
                        observation_distribution[var][j]->Reestimation<int>::type_parametric_estimation(discrete_parametric_process[var+1]->observation[j],
                                                                                                        0, true, OBSERVATION_THRESHOLD);
                        // observation_distribution[var][j]->Reestimation<int>::parametric_estimation(discrete_parametric_process[var+1]->observation[j],
                        //                                                                            0, true, OBSERVATION_THRESHOLD);
                     /*
                     hobservation->update(observation_distribution[var],
                                          MAX((int)(observation_distribution[var][j]->nb_element *
                                                    MAX(sqrt(observation_distribution[var][j]->variance), 1.) * OBSERVATION_COEFF),
                                                        MIN_NB_ELEMENT));
                     observation_likelihood = hobservation->Reestimation<int>::type_parametric_estimation(discrete_parametric_process[var+1]->observation[j],
                                                                                                                            0, true, OBSERVATION_THRESHOLD);
                     */
                     if (observation_likelihood == D_INF)
                        min_likelihood = D_INF;
                     else
                     {
                        discrete_parametric_process[var+1]->observation[j]->computation(trees.characteristics[var]->marginal_distribution->nb_value,
                                                                                        OBSERVATION_THRESHOLD);

                        if (discrete_parametric_process[var+1]->observation[j]->ident == BINOMIAL)
                           for(val = discrete_parametric_process[var+1]->observation[j]->nb_value; val < trees.characteristics[var]->marginal_distribution->nb_value; val++)
                              discrete_parametric_process[var+1]->observation[j]->mass[val] = 0.;
                     }
                  }
               }
               else // (discrete_parametric_process[var+1] == NULL)
               {
                  /*
                  switch (continuous_parametric_process[var+1]->ident)
                  {
                     case GAUSSIAN :
                     {
                        for(j = 0; j < nb_state; j++)
                          continuous_parametric_process[var+1]->observation[j]->location = observation_distribution[var][j]->mean;

                        switch (common_dispersion)
                        {
                           case false :
                           {
                              for(j = 0; j < nb_state; j++)
                                 continuous_parametric_process[var+1]->observation[j]->dispersion = sqrt(observation_distribution[var][j]->variance);
                              break;
                           }

                           case true :
                           {
                              variance[0] = 0.;
                              buff = 0.;

                              for(j = 0; j < nb_state; j++)
                              {
                                 for(val = observation_distribution[var][j]->offset; val < observation_distribution[var][j]->nb_value; val++)
                                 {
                                    diff = val - observation_distribution[var][j]->mean;
                                    variance[0] += observation_distribution[var][j]->frequency[val] * diff * diff;
                                 }
                                 buff += observation_distribution[var][j]->nb_element;
                              }

                              variance[0] /= buff;

                              for(j = 0; j < nb_state; j++)
                                 continuous_parametric_process[var+1]->observation[j]->dispersion = sqrt(variance[0]);
                              break;
                           }
                        } // end switch common_dispersion
                        break;
                     } // end case GAUSSIAN

                     case VON_MISES:
                     {
                        for(j = 0; j < nb_state; j++)
                        {
                           observation_distribution[var][j]->mean_direction_computation(mean_direction[j]);
                           continuous_parametric_process[var+1]->observation[j]->location = mean_direction[j][3];
                        }

                        switch (common_dispersion)
                        {
                           case false :
                           {
                              for(j = 0; j < nb_state; j++)
                                 continuous_parametric_process[var+1]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                              break;
                           }

                           case true :
                           {
                              global_mean_direction = 0.;
                              buff = 0.;

                              for(j = 0; j < nb_state; j++)
                              {
                                 global_mean_direction += observation_distribution[var][j]->nb_element * mean_direction[j][2];
                                 buff += observation_distribution[var][j]->nb_element;
                              }
                              concentration = von_mises_concentration_computation(global_mean_direction / buff);

                              for(j = 0; j < nb_state; j++)
                                 continuous_parametric_process[var+1]->observation[j]->dispersion = concentration;
                              break;
                           }
                        }
                        break;
                     } // end case VON_MISES
                  } // end switch (continuous_parametric_process[var+1]->ident)
                  */
               } // end if (nonparametric_process[var+1] != NULL)
            } // end if (observation_distribution[var] != NULL)
            else // (observation_distribution[var] == NULL)
            {
               /*
               switch (continuous_parametric_process[var+1]->ident)
               {
                   case GAUSSIAN:
                   {
                      for(j = 0; j < nb_state; j++)
                      {
                         mean[j] = 0.;
                         state_frequency[j] = 0.;
                      }

                      switch (tree.get_type(var))
                      {
                         case INT_VALUE :
                         {
                            for(j = 0; j < trees.get_nb_trees; j++)
                               for(pos = 0; pos < length[j]; pos++)
                                 for(m = 0; m < nb_state; m++)
                                 {
                                    mean[m] += state_sequence_count[j][k][m] * int_sequence[j][i][k];
                                    state_frequency[m] += state_sequence_count[j][k][m];
                                 }
                            break;
                         }

                         case REAL_VALUE :
                         {
                            for(j = 0; j < nb_sequence; j++)
                               for(pos = 0; pos < length[j]; pos++)
                                  for(m = 0; m < nb_state; m++)
                                  {
                                    mean[m] += state_sequence_count[j][k][m] * real_sequence[j][i][k];
                                    state_frequency[m] += state_sequence_count[j][k][m];
                                  }
                            break;
                      }

                      for (j = 0;j < nb_state;j++)
                      {
                         mean[j] /= state_frequency[j];
                         continuous_parametric_process[var+1]->observation[j]->location = mean[j];
                      }

                      switch (common_dispersion)
                      {
                         case false :
                         {
                            for(j = 0; j < nb_state; j++)
                               variance[j] = 0.;

                            switch (tree.get_type(var))
                            {
                                case INT_VALUE :
                                {
                                  for (j = 0;j < nb_sequence;j++) {
                                    for (k = 0;k < length[j];k++) {
                                      for (m = 0;m < hsmarkov->nb_state;m++) {
                                        diff = int_sequence[j][i][k] - mean[m];
                                        variance[m] += state_sequence_count[j][k][m] * diff * diff;
                                      }
                                    }
                                  }
                                  break;
                                }

                                case REAL_VALUE :
                                {
                                  for (j = 0;j < nb_sequence;j++) {
                                    for (k = 0;k < length[j];k++) {
                                      for (m = 0;m < hsmarkov->nb_state;m++) {
                                        diff = real_sequence[j][i][k] - mean[m];
                                        variance[m] += state_sequence_count[j][k][m] * diff * diff;
                                      }
                                    }
                                  }
                                  break;
                                }
                            } // end switch (tree.get_type(i))

                            for(j = 0; j < nb_state; j++)
                               continuous_parametric_process[var+1]->observation[j]->dispersion = sqrt(variance[j] / state_frequency[j]);
                            break;
                         }

                         case true :
                         {
                            for(j = 1; j < nb_state; j++)
                               state_frequency[0] += state_frequency[j];

                            variance[0] = 0;

                            switch (trees.get_type(var))
                            {
                               case INT_VALUE : {
                                 for (j = 0;j < nb_sequence;j++) {
                                   for (k = 0;k < length[j];k++) {
                                     for (m = 0;m < hsmarkov->nb_state;m++) {
                                       diff = int_sequence[j][i][k] - mean[m];
                                       variance[0] += state_sequence_count[j][k][m] * diff * diff;
                                     }
                                   }
                                 }
                                 break;
                               }

                               case REAL_VALUE : {
                                 for (j = 0;j < nb_sequence;j++) {
                                   for (k = 0;k < length[j];k++) {
                                     for (m = 0;m < hsmarkov->nb_state;m++) {
                                       diff = real_sequence[j][i][k] - mean[m];
                                       variance[0] += state_sequence_count[j][k][m] * diff * diff;
                                     }
                                   }
                                 }
                                 break;
                               }
                            }

                            for(j = 0; j < nb_state; j++)
                               continuous_parametric_process[var+1]->observation[j]->dispersion = sqrt(variance[0] / state_frequency[0]);
                            break;
                         }
                      } // end switch (common_dispersion)
                      break;
                   } // end case GAUSSIAN

                   case VON_MISES :
                   {
                      for(j = 0; j < nb_state; j++)
                      {
                         mean_direction[j][0] = 0.;
                         mean_direction[j][1] = 0.;

                         state_frequency[j] = 0.;
                      }

                      switch (trees.get_type(var)) {
                      recopier du cas discret ?

                      case INT_VALUE : {
                        for (j = 0;j < nb_sequence;j++) {
                           for (k = 0;k < length[j];k++) {
                             for (m = 0;m < hsmarkov->nb_state;m++) {
                               mean_direction[m][0] += state_sequence_count[j][k][m] * cos(int_sequence[j][i][k] * M_PI / 180);
                               mean_direction[m][1] += state_sequence_count[j][k][m] * sin(int_sequence[j][i][k] * M_PI / 180);
                               state_frequency[m] += state_sequence_count[j][k][m];
                             }
                           }
                         }
                         break;
                       }

                       case REAL_VALUE : {
                         for (j = 0;j < nb_sequence;j++) {
                           for (k = 0;k < length[j];k++) {
                             for (m = 0;m < hsmarkov->nb_state;m++) {
                               switch (hsmarkov->continuous_parametric_process[i + 1]->unit) {
                               case DEGREE :
                                 mean_direction[m][0] += state_sequence_count[j][k][m] * cos(real_sequence[j][i][k] * M_PI / 180);
                                 mean_direction[m][1] += state_sequence_count[j][k][m] * sin(real_sequence[j][i][k] * M_PI / 180);
                                 break;
                               case RADIAN :
                                 mean_direction[m][0] += state_sequence_count[j][k][m] * cos(real_sequence[j][i][k]);
                                 mean_direction[m][1] += state_sequence_count[j][k][m] * sin(real_sequence[j][i][k]);
                                 break;
                               }

                               state_frequency[m] += state_sequence_count[j][k][m];
                             }
                           }
                         }
                         break;
                       }
                      }

                      for (j = 0;j < hsmarkov->nb_state;j++) {
                        mean_direction[j][0] /= state_frequency[j];
                        mean_direction[j][1] /= state_frequency[j];

                        mean_direction[j][2] = sqrt(mean_direction[j][0] * mean_direction[j][0] +
                                                    mean_direction[j][1] * mean_direction[j][1]);


                        if (mean_direction[j][2] > 0.) {
                          mean_direction[j][3] = atan(mean_direction[j][1] / mean_direction[j][0]);

                          if (mean_direction[j][0] < 0.) {
                            mean_direction[j][3] += M_PI;
                          }
                          else if (mean_direction[j][1] < 0.) {
                            mean_direction[j][3] += 2 * M_PI;
                          }

                          if (hsmarkov->continuous_parametric_process[i + 1]->unit == DEGREE) {
                            mean_direction[j][3] *= (180 / M_PI);
                          }
                        }

                        hsmarkov->continuous_parametric_process[i + 1]->observation[j]->location = mean_direction[j][3];
                      }

                      switch (common_dispersion) {

                      case false : {
                        for (j = 0;j < hsmarkov->nb_state;j++) {
                          hsmarkov->continuous_parametric_process[i + 1]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                        }
                        break;
                      }

                      case true : {
                        global_mean_direction = 0.;
                        buff = 0.;

                        for (j = 0;j < hsmarkov->nb_state;j++) {
                          global_mean_direction += state_frequency[j] * mean_direction[j][2];
                          buff += state_frequency[j];
                        }
                        concentration = von_mises_concentration_computation(global_mean_direction / buff);

                        for (j = 0;j < hsmarkov->nb_state;j++) {
                          hsmarkov->continuous_parametric_process[i + 1]->observation[j]->dispersion = concentration;
                        }
                        break;
                      }
                     }
                     break;
               } // end switch (continuous_parametric_process[var+1]->ident)
               */
            } // end if (observation_distribution[var] != NULL)
         } // end if (nonparametric_process[var+1] != NULL)
         } // end  if (!skipped_variables[var])
   // end for var

   if (hobservation != NULL)
   {
      delete hobservation;
      hobservation = NULL;
   }
}

/*****************************************************************
 *
 *  Estimate the initial distribution of a MarkovOutTree
 *  from a MarkovOutTreeData instance.
 *  Update the variable order Markov chains if necessary
 *  (should be called after vomc_estimation)
 *
 **/

void MarkovOutTree::initial_distribution_estimation(const MarkovOutTreeData& trees)
{
   // prerequisites: trees.chain_data and trees.mot_reestimation are non-null
   // and compatible with *this
   // unordered_rank, nb_vomc and vomc_models are known at this step

   register int s;

   if (initial == NULL)
      initial = new double[nb_state];

   for(s = 0; s < nb_state; s++)
      initial[s] = 1. / nb_state;

   // reestimation of the initial distribution
   // does not seem to work if initial[s] is 0
   reestimation(nb_state, trees.chain_data->initial,
                initial, MIN_PROBABILITY, false);

   if (nb_vomc == 1)
      for(s = 0; s < nb_state; s++)
         vomc_models[0]->initial[s] = initial[s];
}

// methods below should be used for testing purposes only...
/*
void MarkovOutTree::get_state_marginal_distribution(const MarkovOutTreeData& trees,
                                     double_array_3d& state_marginal) const;
void MarkovOutTree::get_output_conditional_distribution(const MarkovOutTreeData& trees,
                                         double_array_3d& output_cond) const;

*/

/*****************************************************************
 *
 *  Check whether a MarkovOutTree can be estimated using specified
 *  number of ordered children, number of ordered children
 *  involved in generation processes, integral variables used
 *  as factors in generation processes, a flag on using parent state
 *  as factor and families of generation processes.
 *  A StatError object is updated otherwise.
 *  The number of generation factors, their values, the types of
 *  distributions to be estimated, the list of variables acting
 *  as input factors and not as outputs, and the list of variables
 *  to discard in the analysis for this very reason.
 *
 **/

bool MarkovOutTreeData::markov_out_tree_check_estimation(StatError& error,
                                                         unsigned int inb_vomc,
                                                         unsigned int inb_ordered_children,
                                                         unsigned int inb_children_branching,
                                                         const std::vector<unsigned int>& factor_variables,
                                                         bool parent_dependent,
                                                         const std::vector<int>& generation_types,
                                                         unsigned int& inb_generation,
                                                         unsigned int& inb_generation_values,
                                                         std::vector<int>& generation_types_estimation,
                                                         const GenerationProcess::index_array *ifactor_values,
                                                         std::vector<bool>& skipped_variables) const

{
   // typedef GenerationProcess::index_array gen_index_array;
   // prerequisites:
   //    if inb_factors > 0, characteristics must be computed and have its marginals
   //    if inb_ordered_children > 0, hnb_children must be computed

   bool status = true; //, reestimation_computation = false;
   const unsigned int iunordered_rank = inb_ordered_children+1,
                      inb_factors = factor_variables.size();
   unsigned int fval;
   register int i;
   ostringstream error_message;

   error.init();

   if (state_trees == NULL)
   {
      status = false;
      error_message << STAT_TREES_error[TREESTATR_STATE_TREES];
      error.update((error_message.str()).c_str());
   }

   if (((inb_vomc > 0) && (iunordered_rank < 2))
       || ((inb_vomc < 1) && (iunordered_rank > 1)))
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_NB_VOMC] << " (: "
                    << inb_vomc << ") or " << STAT_TREES_parsing[TREESTATP_ORDERED_CHILDREN]  << " (: "
                    << iunordered_rank << ")";
      error.update((error_message.str()).c_str());
   }
   if (inb_ordered_children > hnb_children->offset)
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_ORDERED_CHILDREN]  << " (: "
                    << iunordered_rank << "): should be at least "
                    << inb_ordered_children << endl;
      error.update((error_message.str()).c_str());
   }

   if ((inb_vomc > 0) && (inb_vomc != 1) && (inb_vomc != _nb_states))
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_NB_VOMC] << "s : "
                    << inb_vomc;
      error.update((error_message.str()).c_str());
   }
   // compute the number of generation distributions
   if (parent_dependent)
   {
      inb_generation_values *= _nb_states;
      inb_generation++;
   }
   else
      if (inb_children_branching + inb_factors == 0)
         inb_generation_values = 0;

   for(i = 0; i < inb_factors; i++)
   {
      // check that factor_variables[i] is integral, starts with value 0
      // and that each value between 0 and max_value is present
      if (factor_variables[i] > _nb_integral)
      {
         status = false;
         error_message << "bad index for variable associated with factor " << i+1 << ": "
                       << factor_variables[i];
         error.update((error_message.str()).c_str());
      }
      else
      {
         if ((_type[factor_variables[i]] != STATE)
             && (_type[factor_variables[i]] != INT_VALUE))
         {
            status = false;
            error_message << STAT_TREES_error[TREESTATR_VARIABLE_TYPE]
                          << " for variable associated with factor " << i+1 << "("
                          << factor_variables[i] << ")";
            error.update((error_message.str()).c_str());
         }
         else
            skipped_variables[factor_variables[i]] = true;
         if ((*ifactor_values)[i] <= _max_value.Int(factor_variables[i]))
         {
            status = false;
            error_message << "Bad maximal value for factor variable: "
                          << factor_variables[i]
                          << "; should be < " << _nb_integral << endl;
            error.update((error_message.str()).c_str());
         }

      }
      if (status)
      {
         assert((characteristics != NULL) && (characteristics[factor_variables[i]] != NULL)
                && (characteristics[factor_variables[i]]->marginal_distribution != NULL));
         for(fval = 0; fval < (*ifactor_values)[i]; fval++)
            if (characteristics[factor_variables[i]]->marginal_distribution->frequency[fval] == 0)
            {
               status = false;
               error_message << "value " << fval << " " << STAT_TREES_error[TREESTATR_NOT_PRESENT]
                             << " in factor associated with variable  " << factor_variables[i] << endl;
               error.update((error_message.str()).c_str());
            }
         if ((*ifactor_values)[i] > 0)
         {
            inb_generation_values *= (*ifactor_values)[i];
            inb_generation++;
         }
         else
         {
            status = false;
            error_message << STAT_TREES_error[TREESTATR_BAD_FACTORS] << "; "
                          << "bad value for factor " << factor_variables[i] << ": "
                          << (*ifactor_values)[i] << endl;
            error.update((error_message.str()).c_str());
         }
      }
   } // end for nb_factor

   if (generation_types.size() != 0)
      if ((generation_types.size() != 1) && (generation_types.size() != inb_generation_values))
      {
         status = false;
         error_message << STAT_TREES_parsing[TREESTATP_NB_GENERATION_PROCESS] << ": "
                       << "bad number of types ( " << generation_types.size() << ");"
                       << "must be 0, 1 or " << inb_generation_values;
         error.update((error_message.str()).c_str());
      }
      else
      {
         if (generation_types.size() == 1)
         {
            generation_types_estimation.resize(inb_generation_values);
            for(i = 0; i < inb_generation_values; i++)
               generation_types_estimation[i] = generation_types[0];
         }
         else
            generation_types_estimation = generation_types;
      }

   inb_generation_values *= pow(_nb_states, inb_children_branching);
   inb_generation += inb_children_branching;

   if ((inb_generation - inb_children_branching - inb_factors < 0)
       || (inb_generation - inb_children_branching - inb_factors > 1))
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_NB_CHILDREN_BRANCHING] << "(: "
                    << inb_generation << ") or " << STAT_TREES_parsing[TREESTATP_NB_CHILDREN_BRANCHING] << "(: "
                    << inb_children_branching << ") or " << STAT_TREES_parsing[TREESTATP_NB_FACTORS]
                    << "(" << inb_factors << ")";
      error.update((error_message.str()).c_str());
   }
   if (inb_children_branching > iunordered_rank - 1)
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_ORDERED_CHILDREN]  << " (: "
                    << iunordered_rank - 1 << ") or " << STAT_TREES_parsing[TREESTATP_NB_CHILDREN_BRANCHING] << "(: "
                    << inb_children_branching << ")";
      error.update((error_message.str()).c_str());
   }

   // ensure that a model exists for the ordered children
   if ((inb_generation == 0) && (inb_vomc == 0))
   {
      status = false;
      error_message << "At least one model for ordered children or one model "
                    << "for unordered children must be defined";
      error.update((error_message.str()).c_str());
   }

   return status;
}


/*****************************************************************
 *
 *  Estimate a MarkovOutTree using a StatError object,
 *  the number of ordered children, the number of ordered children
 *  involved in generation processes, the integral variables used
 *  as factors in generation processes, a flag on using parent state
 *  as factor, the families of generation processes,
 *  minimal possible inferior bound and a flag on estimating this bound
 *
 **/

MarkovOutTree* MarkovOutTreeData::markov_out_tree_estimation(StatError& error, unsigned int inb_vomc,
                                                             unsigned int inb_ordered_children,
                                                             unsigned int inb_children_branching,
                                                             const std::vector<unsigned int>& factor_variables,
                                                             bool parent_dependent,
                                                             const std::vector<int>& generation_types,
                                                             unsigned int generation_min_inf_bound,
                                                             bool generation_min_inf_bound_flag) const
{
   typedef GenerationProcess::index_array gen_index_array;
   // prerequisites:
   //    if inb_factors > 0, characteristics must be computed and have its marginals
   //    if inb_ordered_children > 0, hnb_children must be computed

   bool status = true, reestimation_computation = false;
   const unsigned int iunordered_rank = inb_ordered_children+1,
                      inb_factors = factor_variables.size();
   unsigned int inb_generation = 0, fval, inb_generation_values = 1;
   std::vector<int> generation_types_estimation(0);
   std::vector<bool> skipped_variables(_nb_integral+_nb_float, false);
   gen_index_array *ifactor_values = new gen_index_array(0, inb_factors);
   register int i;
   ostringstream error_message;
   MarkovOutTree *res = NULL;
   MarkovOutTreeData *markov_data = NULL;

   status = markov_out_tree_check_estimation(error, inb_vomc, inb_ordered_children,
                                             inb_children_branching, factor_variables,
                                             parent_dependent, generation_types,
                                             inb_generation,
                                             inb_generation_values,
                                             generation_types_estimation,
                                             ifactor_values,
                                             skipped_variables);

   if (status)
   {
      markov_data = new MarkovOutTreeData(*this, false, false);
      if (markov_data->mot_reestimation == NULL)
         reestimation_computation = true;
      else
      {
         if ((markov_data->mot_reestimation->get_nb_children_branching() != inb_children_branching) ||
             (markov_data->mot_reestimation->get_nb_factors() != inb_factors) ||
             (markov_data->mot_reestimation->get_nb_generation() != inb_generation) ||
             (markov_data->mot_reestimation->get_unordered_rank() != iunordered_rank) ||
             (markov_data->mot_reestimation->get_nb_vomc() != inb_vomc))
            reestimation_computation = true;
      }
      if (reestimation_computation)
      {
         markov_data->min_max_value_computation();
         markov_data->build_chain_reestimation('o', _nb_states);
         markov_data->build_markov_reestimation(inb_children_branching, inb_factors, ifactor_values,
                                                inb_generation, iunordered_rank, inb_vomc);
         markov_data->build_observation_frequency_distribution();
      }
      res = new MarkovOutTree(_nb_states, _nb_integral+_nb_float-inb_factors);

      /* if (markov_data->hsize != NULL)
         res->nonparametric_process[0]->size = new FrequencyDistribution(*markov_data->hsize); */
      // markov_data, nonparametric_process, discrete_parametric_process,
      // continuous_parametric_process, unordered_rank, nb_vomc, vomc_models
      // and generation_distributions remain to be built at this step

      if (inb_vomc > 0)
         res->vomc_estimation(error, *markov_data);
      else
         res->unordered_rank = 1;
      if ((error.get_nb_error() == 0) && (inb_generation > 0))
         res->generation_process_estimation(error, *markov_data, generation_types_estimation,
                                            generation_min_inf_bound, generation_min_inf_bound_flag);
      if ((error.get_nb_error() == 0) && (res->get_nb_output_process() > 0))
         res->output_process_estimation(error, *markov_data, skipped_variables);
      if (error.get_nb_error() == 0)
         // estimate initial distribution after vomcs, and update if necessary
         res->initial_distribution_estimation(*markov_data);

      if (error.get_nb_error() > 0)
      {
         delete markov_data;
         markov_data = NULL;
         delete res;
         res = NULL;
         error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
      }
      else
      {
         markov_data->likelihood = res->likelihood_computation(*markov_data);
         markov_data->build_size_frequency_distribution();
         res->markov_data = markov_data;
         // etc.
      }
   }

   delete ifactor_values;
   ifactor_values = NULL;

   return res;
}
