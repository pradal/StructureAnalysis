/****************************************************************
 *
 *  Test of the estimation algorithms for hidden Markov
 *  out-trees as defined in hidden_markov_out_tree.h:
 *  checking the syntax and the algorithm basic properties
 */

#include "tree_simple.h"
#include "tree_traits.h"
#include "basic_visitors.h"
#include "stat_tools.h"
#include "curves.h"
#include "markovian.h"
#include "sequences.h"
#include "semi_markov.h"
#include "hidden_semi_markov.h"
#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"
#include "int_trees.h"
#include "hidden_markov_tree.h"
#include "hidden_markov_out_tree.h"

using namespace Stat_trees;

int main(void)
{
   typedef Hidden_markov_tree_data::tree_type hmtd_tree_type;
   typedef Hidden_markov_tree_data::pt_observed_trees_array pt_observed_trees_array;
   typedef Hidden_markov_tree_data::observed_trees observed_trees;
   typedef Hidden_markov_tree_data::value value;
   typedef Hidden_markov_tree_data::vertex_iterator vertex_iterator;
   typedef Stat_trees::Hidden_markov_tree::double_array_3d double_array_3d;
   typedef Stat_trees::Hidden_markov_tree::double_array_4d double_array_4d;
   typedef generic_visitor<hmtd_tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   register int t, j, i, nb_states, nb_trees= 5;
   bool plot_ok= false;
   unsigned int u;
   double likelihood, entropy1, entropy2;
   value default_value;
   visitor v;
   vertex_iterator it, end;
   vertex_array va;
   Format_error error;
   Hidden_markov_out_tree *hmot= NULL, *hmot2= NULL, *hmot_init= NULL;
   Hidden_markov_tree_data *hmtd= NULL, *hmtdman= NULL,
                           *hmtdtmp= NULL;
   Trees *ctrees= NULL;
   hmtd_tree_type **ptrees= NULL;
   observed_trees *merged_trees= NULL;
   std::vector<Hidden_markov_tree_data*> tree_list;
   hmtd_tree_type ctree;
   // const char * hsmcpath= "./laricio_3.hsc";
   const char * hmotpath= "./hmot_np_2s.hmt";
   const char * hmotinitpath= "./hmot_np_init_2s.hmt";
   const int size= 15, nb_children_max= 2,
             nb_iterations= 5;
   double_array_3d state_marginal= NULL, output_cond= NULL,
          upward_prob= NULL, upward_parent_prob= NULL,
          downward_prob= NULL,
          state_entropy= NULL;
   double_array_4d downward_pair_prob= NULL;
   double **sum_state_marginal= NULL;
   // int *rejected_variables= NULL;

   // reading and printing a hidden Markov out tree
   hmot= hidden_markov_out_tree_ascii_read(error, hmotpath);
   cout << error;

   if (hmot != NULL)
      hmot->ascii_write(cout, false);

   cout << endl;

   if (hmot != NULL)
   {
      // simulation of hidden Markov out trees
      hmtdtmp= hmot->simulation(error, nb_trees, size, nb_children_max);
      cout << error;
      // add a tree of size 2
      hmtdman= hmot->simulation(error, 1, 2, 1);
      cout << error;
      tree_list.push_back(hmtdman);
      hmtd= hmtdtmp->merge(error, tree_list);
      delete hmtdtmp;
      hmtdtmp= NULL;
      delete hmtdman;
      hmtdman= NULL;

      if (hmtd != NULL)
      {
         nb_trees= hmtd->get_nb_trees();
         // print simulated trees

         cout << "Simulated trees : " << endl;
         for(t= 0; t < nb_trees; t++)
         {
            ctree= *(hmtd->get_tree(t));
            ctree.display(cout, 0);
         }

         cout << endl << "Simulated hidden trees : " << endl;
         for (t= 0; t < nb_trees; t++)
            (hmtd->get_state_tree(t))->display(cout, 0);
         cout << endl;

         // compute the state marginal distributions
         hmot->get_state_marginal_distribution(*hmtd, state_marginal);
         sum_state_marginal= new double*[nb_trees];

         hmot->get_output_conditional_distribution(*hmtd, output_cond);

         nb_states= hmot->get_nb_state();


         ptrees= new hmtd_tree_type*[nb_trees];
         default_value.reset(0, 1);
         // default_value.Double(0)= .0;

         for(t= 0; t < nb_trees; t++)
         {
            if (t == 0)
               cout << "Marginal distributions for tree number " << t+1  << " : " << endl;
            ptrees[t]= new hmtd_tree_type(*((hmtd->get_state_tree(t))->get_structure()),
                                     default_value);
            va= v.get_breadthorder(*(ptrees[t]));
            sum_state_marginal[t]= new double[va.size()];
            for(j= 0; j < hmot->get_nb_state(); j++)
            {
               if (t == 0)
                  cout << "state " << j << " : " << endl;
               for(u= 0; u < va.size(); u++)
               {
                  default_value.Double(0)= state_marginal[t][j][va[u]];
                  ptrees[t]->put(va[u], default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());
            }
            if (t == 0)
               cout << "Output conditional distributions for tree number " << t+1  << " : " << endl;

            for(j= 0; j < hmot->get_nb_state(); j++)
            {
               if (t == 0)
                  cout << "state " << j << " : " << endl;
               for(u= 0; u < ptrees[t]->get_size(); u++)
               {
                  default_value.Double(0)= output_cond[t][j][u];
                  ptrees[t]->put(u, default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());

            }
            for(u= 0; u < va.size(); u++)
            {
               sum_state_marginal[t][u]= .0;
               for(j= 0; j < hmot->get_nb_state(); j++)
                  sum_state_marginal[t][u]+= state_marginal[t][j][u];
            }
         }

         likelihood= hmot->get_upward_step(*hmtd, upward_prob, upward_parent_prob,
                                           state_entropy, state_marginal,
                                           output_cond, entropy1);
         hmot->get_downward_step(*hmtd, downward_prob, downward_pair_prob,
                                 upward_prob, upward_parent_prob,
                                 state_marginal, output_cond, entropy2);
         cout << "Log-likelihood of the parameter : " << likelihood << endl;

         for(t= 0; t < nb_trees; t++)
         {
            if (t == 0)
               cout << "State conditional distributions for tree number " << t+1  << " : " << endl;
            for(j= 0; j < hmot->get_nb_state(); j++)
            {
               if (t == 0)
                  cout << "state " << j << " : " << endl;
               tie(it, end)= ptrees[t]->vertices();
               while (it < end)
               {
                  default_value.Double(0)= downward_prob[t][j][*it];
                  ptrees[t]->put(*it++, default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());
            }
         }

         default_value.reset(0, nb_states);
         for(t= 0; t < nb_trees; t++)
         // should be done for t < nb_trees
         {
            if (t == 0)
               cout << "Conditional distributions for pairs of states for tree number " << t+1  << " : " << endl;
            for(j= 0; j < nb_states; j++)
            {
               if (t == 0)
                  cout << "child state " << j << " : " << endl;
               tie(it, end)= ptrees[t]->vertices();
               while (it < end)
               {
                  for(i= 0; i < nb_states; i++)
                     default_value.Double(i)= downward_pair_prob[t][i][j][*it];
                  ptrees[t]->put(*it++, default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());
            }
         }

         delete hmot2;
         cout << "Testing estimation with initial hmt specification." << endl;
         cout << "Initialization HMT " << endl;
         hmot_init= hidden_markov_out_tree_ascii_read(error, hmotinitpath);
         cout << error;
         if (hmot_init != NULL)
            hmot_init->ascii_write(cout, false);

         hmot2= hmtd->hidden_markov_out_tree_estimation(error, cout, *hmot_init,
                                                        true, VITERBI, 1, true);
         cout << error;

         if (hmot2 != NULL)
         {
            likelihood= hmot2->likelihood_computation(*hmtd);
            delete hmot_init;
            hmot_init= NULL;
            cout << "estimated HMT (one iteration)" << endl;
            // cout << *hmot2 << endl;
            hmot2->ascii_write(cout, false); //true
            for(i= 0; i < nb_iterations; i++)
            {
               hmot_init= hmot2;
               hmot2= hmtd->hidden_markov_out_tree_estimation(error, cout, *hmot_init,
                                                              true, VITERBI, 1, true);
               delete hmot_init;
               hmot_init= NULL;
               cout << "relative likelihood increase at iteration " << i << ": "
                    << (hmot2->likelihood_computation(*hmtd) - likelihood)/abs(likelihood) << endl;
               likelihood= hmot2->likelihood_computation(*hmtd);
            }
            cout << endl << "estimated HMT (" << nb_iterations << " iterations)" << endl;
            cout << *hmot2;
   #if 1
            cout << endl << "Testing state profile plot" << endl;
            plot_ok= hmot2->state_profile_plot_write(error, "/home/durand/ftmp",
                                                     4, 8, NULL);
            if (!plot_ok)
               cout << error;
            delete hmot2;
            hmot2= NULL;
   #endif
         }
         cout << endl;

         cout << "Testing estimation with initial self-transition probability "
              << "specification." << endl;

         hmot2= hmtd->hidden_markov_out_tree_estimation(error, cout, 'o', 2, false,
                                                        true, VITERBI, 0.99999, 1);
         cout << error;

         if (hmot2 != NULL)
         {
            likelihood= hmot2->likelihood_computation(*hmtd);
            hmot_init= NULL;
            cout << "estimated HMT (one iteration)" << endl;
            cout << *hmot2;
            for(i= 0; i < nb_iterations; i++)
            {
               hmot_init= hmot2;
               hmot2= hmtd->hidden_markov_out_tree_estimation(error, cout, *hmot_init,
                                                              true, VITERBI, 1, true);
               delete hmot_init;
               hmot_init= NULL;
               cout << "relative likelihood increase at iteration " << i << ": "
                    << (hmot2->likelihood_computation(*hmtd) - likelihood)/abs(likelihood)<< endl;
               likelihood= hmot2->likelihood_computation(*hmtd);
            }
            cout << endl << "estimated HMT (" << nb_iterations << " iterations)" << endl;
            cout << *hmot2;
            delete hmot2;
         }

         cout << endl;

         for(t= 0; t < nb_trees; t++)
         {
            for(j= 0; j < nb_states; j++)
            {
               delete [] state_marginal[t][j];
               delete [] upward_prob[t][j];
               delete [] upward_parent_prob[t][j];
               delete [] downward_prob[t][j];
               for(i= 0; i < nb_states; i++)
                  delete [] downward_pair_prob[t][j][i];
               delete [] downward_pair_prob[t][j];
               delete [] state_entropy[t][j];
            }
            for(u= 0; u < ptrees[t]->get_size(); u++)
               if (fabs(sum_state_marginal[t][u] - 1.) > DOUBLE_ERROR)
                  cout << "State marginal distribution for tree " << t <<
                       " and node " << u << " sums to " << sum_state_marginal[t][u] << endl;
            delete [] state_marginal[t];
            delete [] output_cond[t];
            delete [] upward_prob[t];
            delete [] upward_parent_prob[t];
            delete [] downward_prob[t];
            delete [] downward_pair_prob[t];
            delete [] sum_state_marginal[t];
            delete [] state_entropy[t];
            delete ptrees[t];
         }

         delete [] ptrees;
         delete [] state_marginal;
         delete [] sum_state_marginal;
         delete [] output_cond;
         delete [] upward_prob;
         delete [] upward_parent_prob;
         delete [] downward_prob;
         delete [] downward_pair_prob;
         delete [] state_entropy;
         delete hmot2;
         delete hmtd;
      }
      delete hmot;
   }
   return 0;
}
