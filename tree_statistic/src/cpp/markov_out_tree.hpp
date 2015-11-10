/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       $Source$
 *       $Id: markov_out_tree.hpp 9554 2011-02-03 17:10:30Z jbdurand $
 *
 *       Forum for V-Plants developers:
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

#ifndef MARKOV_OUT_TREE_HPP
#define MARKOV_OUT_TREE_HPP

#include <assert.h>

/*****************************************************************
 *
 *  Constructor of MarkovOutTreeReestimation class using the number
 *  of states, the number of children states, of factors
 *  and parent children determining the number inb_generation
 *  of branching processes, a collection of generation distributions
 *  and the sequences of ordered children
 *
 **/

template<typename Type>
MarkovOutTreeReestimation<Type>::MarkovOutTreeReestimation(int inb_state,
                                                           unsigned int inb_children_branching,
                                                           unsigned int inb_factors,
                                                           const index_array *ifactor_values,
                                                           unsigned int inb_generation,
                                                           unsigned int iunordered_rank,
                                                           unsigned int inb_vomc,
                                                           mvHistogram * const * pgeneration,
                                                           const sequence_analysis::MarkovianSequences * const * pseq)

 : nb_state(inb_state)
 , nb_children_branching(inb_children_branching)
 , nb_factors(inb_factors)
 , factor_values(NULL)
 , nb_generation_values(1)
 , nb_generation(inb_generation)
 , factor_base(0)
 , unordered_rank(iunordered_rank)
 , nb_vomc(inb_vomc)
 , seq(NULL)
 , generation_reestim(NULL)
{
   unsigned int j, s;

   if (nb_state > 0)
   {
      if (nb_generation == nb_factors + nb_children_branching + 1)
      // parent state is a factor
         nb_generation_values *= nb_state;
      for(s = 0; s < nb_children_branching; s++)
         nb_generation_values *= nb_state;

      if (nb_factors > 0)
      {
         factor_values = new index_array(nb_factors);
         for(s = 0; s < nb_factors; s++)
         {
            (*factor_values)[s] = (*ifactor_values)[s];
            nb_generation_values *= (*factor_values)[s];
         }
      }

      if (pgeneration != NULL)
      {
         generation_reestim = new mvHistogram*[nb_generation_values];
         for(s = 0; s < nb_generation_values; s++)
         {
            if (pgeneration[s] != NULL)
               generation_reestim[s] = new mvHistogram(*pgeneration[s]);
            else
               generation_reestim[s] = NULL;
         }
      }
      else
         generation_reestim = NULL;

      factor_base.resize(nb_generation, 0);
      s = 0;
      if (nb_generation == nb_children_branching + nb_factors + 1)
      // parent state is a factor
      {
         factor_base[s] = nb_state;
         s++;
      }
      j = 0;
      while (j < nb_children_branching)
      {
         factor_base[s] = nb_state;
         s++;
         j++;
      }
      j = 0;
      while (j < nb_factors)
      {
         factor_base[s] = (*factor_values)[j];
         s++;
         j++;
      }

      if ((pseq != NULL) && (nb_vomc >0))
      {
         seq = new sequence_analysis::MarkovianSequences*[nb_vomc];
         for(s = 0; s < nb_vomc; s++)
            seq[s] = new sequence_analysis::MarkovianSequences(*pseq[s]);
      }
   }
}

/*****************************************************************
 *
 *  Copy constructor for MarkovOutTreeReestimation class
 *
 **/

template<typename Type>
MarkovOutTreeReestimation<Type>::MarkovOutTreeReestimation(const MarkovOutTreeReestimation<Type>& process)
 : nb_state(0)
 , nb_children_branching(0)
 , nb_factors(0)
 , factor_values(NULL)
 , nb_generation_values(0)
 , nb_generation(0)
 , factor_base(0)
 , unordered_rank(0)
 , nb_vomc(0)
 , seq(NULL)
 , generation_reestim(NULL)
{ copy(process); }


/*****************************************************************
 *
 *  Destructor for MarkovOutTreeReestimation class
 *
 **/

template<typename Type>
MarkovOutTreeReestimation<Type>::~MarkovOutTreeReestimation()
{ remove(); }

/*****************************************************************
 *
 *  Code a configuration of factors of MarkovOutTreeReestimation class
 *  with an unsigned int
 *
 **/

template<typename Type>
unsigned int MarkovOutTreeReestimation<Type>::code(const index_array& fact) const
{
   unsigned int code = 0, base = 1;
   register int j;

   assert(fact.size() == nb_generation);

   /* for(j = nb_generation; j > 0; j--)
   {
      assert(fact[j-1] < factor_base[j-1]);
      code += fact[j-1] * base;
      base *= factor_base[j-1];
   }*/
   for(j = 0; j < nb_generation; j++)
   {
      assert(fact[j] < factor_base[j]);
      code += fact[j] * base;
      base *= factor_base[j];
   }

   return code;
}

/*****************************************************************
 *
 *  Decode a configuration of factors of MarkovOutTreeReestimation class
 *  coded with an unsigned int
 *
 **/

template<typename Type> typename
MarkovOutTreeReestimation<Type>::index_array
MarkovOutTreeReestimation<Type>::decode(unsigned int code) const
{
   unsigned int base;
   register int j;
   index_array decodev(nb_generation, 0);

   assert(code < nb_generation_values);

   /* base = factor_base[nb_generation-1];
   for(j = nb_generation; j > 0; j--)
   {
      decodev[j-1] = code % base;
      code = code / base;
      base *= factor_base[j-1];
   }*/

   base = nb_generation_values;

   for(j = nb_generation; j > 0; j--)
   {
      base = base / factor_base[j-1];
      decodev[j-1] = code / base;
      code = code % base;
   }

   return decodev;
}

/*****************************************************************
 *
 *  Code a configuration of factors of MarkovOutTreeReestimation class
 *  with an unsigned int
 *
 **/

template<typename Type>
unsigned int MarkovOutTreeReestimation<Type>::code_error(StatError& error,
                                                         const index_array& fact) const
{
   bool status = true;
   register int j;

   error.init();

   if (fact.size() != nb_generation)
   {
      status = false;
      error.update(STAT_TREES_error[TREESTATR_BAD_FACTORS]);
   }

   for(j = nb_generation; j > 0; j--)
      if (fact[j-1] >= factor_base[j-1])
      {
         status = false;
         error.update(STAT_TREES_error[TREESTATR_BAD_FACTORS]);
      }

   if (status)
      return code(fact);
   else
      return 0;
}

/*****************************************************************
 *
 *  Return the number of states in MarkovOutTreeReestimation
 *
 **/

template<typename Type>
int MarkovOutTreeReestimation<Type>::get_nb_state() const
{ return nb_state; }

/*****************************************************************
 *
 *  Return the rank (strictly) beyond which children are unordered
 *  in MarkovOutTreeReestimation
 *
 **/

template<typename Type>
int MarkovOutTreeReestimation<Type>::get_unordered_rank() const
{ return unordered_rank; }

/*****************************************************************
 *
 *  Return the number of variable order Markov chains
 *  in MarkovOutTreeReestimation
 *
 **/

template<typename Type>
int MarkovOutTreeReestimation<Type>::get_nb_vomc() const
{ return nb_vomc; }

/*****************************************************************
 *
 *  Return the number of ordered children involved
 *  in the branching processes as factors (due to dependencies
 *  regarding these states, excluding parent state)
 *  in MarkovOutTreeReestimation
 *
 **/

template<typename Type>
unsigned int MarkovOutTreeReestimation<Type>::get_nb_children_branching() const
{ return nb_children_branching; }

/*****************************************************************
 *
 *  Return the number of other factors involved
 *  in the branching process in MarkovOutTreeReestimation
 *
 **/

template<typename Type> unsigned int
MarkovOutTreeReestimation<Type>::get_nb_factors() const
{ return nb_factors; }


/*****************************************************************
 *
 *  Return the number of possible values for each factor
 *  in MarkovOutTreeReestimation
 *
 **/

template<typename Type>
typename MarkovOutTreeReestimation<Type>::index_array
MarkovOutTreeReestimation<Type>::get_factor_values() const
{
   index_array res(0);

   if (factor_values != NULL)
      res = *factor_values;

   return res;
}

/*****************************************************************
 *
 *  Return the number of generation processes in MarkovOutTreeReestimation
 *  (combining all possible values of every factor)
 *
 **/

template<typename Type> unsigned int
MarkovOutTreeReestimation<Type>::get_nb_generation_values() const
{ return nb_generation_values; }


/*****************************************************************
 *
 *  Return the number of factors for generation processes
 *  (including parent state) in MarkovOutTreeReestimation
 *
 **/

template<typename Type> unsigned int
MarkovOutTreeReestimation<Type>::get_nb_generation() const
{ return nb_generation; }

/*****************************************************************
 *
 *  Get every possible combination of factors for generation processes
 *  (including parent state). Parent state comes first (?) if included.
 *
 **/

template<typename Type>
std::vector<GenerationProcess::index_array>
MarkovOutTreeReestimation<Type>::get_factor_combinations() const
{
   int icode;
   std::vector<GenerationProcess::index_array> res(0);

   if (generation_reestim != NULL)
      for (icode = 0; icode < get_nb_generation_values (); icode++)
         // if (generation_reestim[icode]->get_nb_different_values() > 0)
         if (generation_reestim[icode]->get_total() > 0)
         {
            res.push_back(decode(icode));
            if (nb_generation > nb_factors + nb_children_branching)
               assert(res[res.size()-1][0] < nb_state);
         }
   return res;
}



/*****************************************************************
 *
 *  State permutation for MarkovOutTreeReestimation class
 *  (checking validity of permutation)
 *
 **/

template<typename Type> void
MarkovOutTreeReestimation<Type>::state_permutation(StatError& error, int* perm) const
{}

/*****************************************************************
 *
 *  Get marginal of class MarkovOutTreeReestimation
 *  for a given configuration of factors (referred to by a code)
 *  and for a given state
 *  (return a pointer; object should not be deallocated)
 *
 **/

/* template<typename Type> const Reestimation< Type >*
MarkovOutTreeReestimation<Type>::get_marginal_ptr(unsigned int code, int istate) const
{
   assert(istate < nb_state);
   assert(generation_reestim != NULL);
   assert(code < nb_generation_values);
   assert(generation_reestim[code]->marginals != NULL);

   return generation_reestim[code]->get_marginal_ptr(istate);
}  */

/*****************************************************************
 *
 *  Get marginal of class MarkovOutTreeReestimation
 *  for a given configuration of factors (referred to by a code)
 *  and for a given state
 *  (a new instance is allocated)
 *
 **/
template<typename Type> Statiskit::Marginal::Univariate::CountTable*
MarkovOutTreeReestimation<Type>::get_marginal(unsigned int code, int istate) const
{
   uHistogram *res;

   assert(generation_reestim[code] != NULL);

   res = Statiskit::marginalizing(*generation_reestim[code], istate);

   return res;
}

/*****************************************************************
 *
 *  Get marginal FrequencyDistribution for class MarkovOutTreeReestimation
 *  for a given configuration of factors (referred to by a code)
 *  and for a given state
 *  (a new instance is allocated)
 *
 **/

/* template<typename Type> FrequencyDistribution*
MarkovOutTreeReestimation<Type>::get_marginal_frequency(unsigned int code, int istate) const

{
   assert(istate < nb_state);
   assert(generation_reestim != NULL);
   assert(code < nb_generation_values);
   assert(generation_reestim[code]->marginals != NULL);

   return generation_reestim[code]->get_marginal_frequency(istate);
}*/

/*****************************************************************
 *
 *  Update frequency for a given value in MarkovOutTreeReestimation
 *  for a given configuration of factors
 *  and for a given configuration of states
 *  Offset is substracted in code, not in marginals
 *
 **/

template<typename Type> void
MarkovOutTreeReestimation<Type>::update(const index_array& fact,
                                        const index_array& states, Type freq)
{
   int icode = code(fact);
   Statiskit::Marginal::Multivariate::CountTable::key_type value;

   assert(generation_reestim != NULL);

   value = Stat_histogram_value(states);

   generation_reestim[icode]->add(value, freq);
}

/*****************************************************************
 *
 *  Get frequency for a given value in MarkovOutTreeReestimation
 *  for a given configuration of factors and for a given state
 *
 **/

template<typename Type> Type
MarkovOutTreeReestimation<Type>::get_frequency(const index_array& fact,
                                               const index_array& states) const
{
   register int i;
   int icode = code(fact);
   mvHistogram::const_iterator it;
   mvHistogram::key_type value;

   assert(generation_reestim != NULL);

   for(i=0; i < states.size(); i++)
      value[i] = states[i];

   it = generation_reestim[icode]->cfind(value);

   return (*it).second;
}

/*****************************************************************
 *
 *  Get reestimation for generation process of a MarkovOutTreeReestimation
 *  (return a pointer; object should not be deallocated)
 *
 **/

template<typename Type> Statiskit::Marginal::Multivariate::CountTable**
MarkovOutTreeReestimation<Type>::get_generation_reestim_ptr() const
{ return generation_reestim; }

/*****************************************************************
 *
 *  Get reestimation for generation process of a MarkovOutTreeReestimation
 *  using a configuration of factors
 *  (return a pointer; object should not be deallocated)
 *
 **/

template<typename Type> Statiskit::Marginal::Multivariate::CountTable*
MarkovOutTreeReestimation<Type>::get_generation_reestim_ptr(StatError& error,
                                                            const index_array& fact) const
{
   bool status = true;
   unsigned int icode = 0;

   DiscreteMultivariateReestimation<Type> *res = NULL;

   if (error.get_nb_error() > 0)
      status = false;

   if (status)
   {
      icode = code_error(error, fact);
      if (generation_reestim[icode] == NULL)
      {
         status = false;
         error.update(STAT_TREES_error[TREESTATR_CHARACTERISTICS_NOT_COMPUTED]);
      }
      else
         return generation_reestim[icode];
   }
}

/*****************************************************************
 *
 *  Initialize empirical generation distributions
 *  of a MarkovOutTreeReestimation using the minimal and maximal values
 *
 **/

template<typename Type> void
MarkovOutTreeReestimation<Type>::init(const std::vector<unsigned int>& imin_offset,
                                      const std::vector<unsigned int>& imax_nb_values)
{
   unsigned int s;
   const bool const_min_offset = (imin_offset.size() == 1);
   const bool const_max_values = (imax_nb_values.size() == 1);
   // std::vector<Type> scalars(nb_state, Statiskit::INTEGER);
   Statiskit::storages_type scalars(nb_state, Statiskit::INTEGER);

   assert(((imin_offset.size() == nb_generation_values) && (imax_nb_values.size() == nb_generation_values))
           || (const_min_offset && const_max_values));

   if (nb_generation_values > 0)
   {
      generation_reestim = new mvHistogram*[nb_generation_values];
      if ((const_min_offset) && (const_max_values))
      {
         for(s = 0; s < nb_generation_values; s++)
            generation_reestim[s] = new mvHistogram(scalars); // nb_state, imin_offset[0],
                                                    // imax_nb_values[0]);
      }
      else
      {
         for(s = 0; s < nb_generation_values; s++)
            generation_reestim[s] = new mvHistogram(scalars); //nb_state, imin_offset[s],
                                                    //imax_nb_values[s]);
      }

   }

}

/*****************************************************************
 *
 *  Compute / update frequencies for every configuration of states
 *  in MarkovOutTreeReestimation using a given dataset,
 *  the index of the considered tree, and the indices of variables
 *  associated with factors
 *
 **/

template<typename Type> void
MarkovOutTreeReestimation<Type>::build_characteristics(const MarkovOutTreeData& data,
                                                       int index,
                                                       const std::vector<unsigned int>* factor_indices)
{
   // nb_state, nb_children_branching, nb_factors, factor_values,
   // nb_generation_values, nb_generation, and factor_base
   // are supposed known and compliant with the data at this point

   // factor_indices is supposed to be compliant with nb_factors

   // cdeque[pstate] contains every sequences of ordered children
   // with parent pstate
   typedef MarkovOutTreeData::state_tree_type state_tree_type;
   typedef state_tree_type::vertex_descriptor vid;
   typedef state_tree_type::value state_value;
   typedef MarkovOutTreeData::tree_type::value value;
   typedef generic_visitor<state_tree_type> visitor;
   typedef visitor::vertex_array vertex_array;
   typedef state_tree_type::children_iterator children_iterator;

   // is_pstate_factor is true iif parent state is a factor in generation process
   const bool is_pstate_factor = (nb_generation == nb_factors + nb_children_branching + 1);
   const unsigned int nb_trees = data.get_nb_trees();
   unsigned int t, current_size, u, current_nb_children_branching,
            current_rank, code, current_nb_seq, current_factor, current_factor_pos;
   register int i, seq_id, pos, pstate, state;
   const int index_parameter_type = sequence_analysis::IMPLICIT_TYPE;
   int *seq_identifier = NULL, // seq_identifier[id]
       *seq_length = NULL; // seq_length[id]
   int ***int_sequence; // int_sequence[id][var][pos]
   std::deque<std::vector<int> > *cdeque = NULL; // cdeque[pstate][id][pos]
   std::vector<int> children_states_seq(0);
   visitor v;
   vertex_array va;
   state_value val;
   value fval;

   // configuration of preceding ordered states, including parent and factors
   index_array pstates(nb_generation, 0),
   // configuration of unordered states
               state_config(nb_state, 0),
   // configuration of factors, which is (last) part of pstates
               factor_config(nb_factors, 0);
   vid cnode;
   children_iterator ch_it, ch_end;
   Typed_edge_one_int_tree *current_tree;
   sequence_analysis::Sequences *cseq = NULL;


   if ((unordered_rank > 1) && (seq == NULL))
      seq = new sequence_analysis::MarkovianSequences*[nb_vomc];

   if (unordered_rank > 1)
      cdeque = new std::deque<std::vector<int> >[nb_state];

   if (nb_factors > 0)
      assert((factor_indices != NULL) && (factor_indices->size() == nb_factors));

   for(t = 0; t < nb_trees; t++)
      if ((index == I_DEFAULT) || (index == t))
      {
         current_tree = data.state_trees[t];
         current_size = current_tree->get_size();
         va = v.get_breadthorder(*current_tree); // downward recursion
         for(u = 0; u < current_size; u++)
         {
            cnode = va[u];
            val = current_tree->get(cnode);
            pstate = val.Int();
            Tree_tie::tie(ch_it, ch_end) = current_tree->children(cnode);
            // current number of ordered children taken into account in factors
            current_nb_children_branching = 0;
            // current position of factor in list of factors
            current_factor_pos = 0;
            current_rank = 0;
            children_states_seq.clear();
            if (is_pstate_factor)
               // add parent state to the list of preceding states
               pstates[current_factor_pos++] = pstate;
            if (nb_factors > 0)
            {
               fval = data.trees[t]->get(cnode);
               for(current_factor = 0; current_factor < nb_factors; current_factor++)
                  factor_config[current_factor] = fval.Int((*factor_indices)[current_factor]);
            }
            while(ch_it < ch_end)
            {
               val = current_tree->get(*ch_it++);
               state = val.Int();
               if (current_nb_children_branching < nb_children_branching)
               {
                  // I_DEFAULT states should be handled otherwise ?
                  current_nb_children_branching++;
                  pstates[current_factor_pos++] = state;
               }
               if (current_rank < (unordered_rank - 1)) // ordered vertex
               {
                  // create list of ordered children to build
                  // sequence_analysis::MarkovianSequences later
                  children_states_seq.push_back(state);
                  current_rank++;
               }
               else // unordered ordered vertex
                  state_config[state]++;
            }
            if (unordered_rank > 1) // ordered children exist
            {
               // create list of ordered children to build
               // sequence_analysis::MarkovianSequences later
               while (current_rank < (unordered_rank - 1))
               {
                  // complete missing ordered states with default values
                  children_states_seq.push_back(I_DEFAULT);
                  current_rank++;
               }
               cdeque[pstate].push_back(children_states_seq);
            }
            // add factors
            for(current_factor = 0; current_factor < nb_factors; current_factor++)
               pstates[current_factor_pos++] = factor_config[current_factor];
            if (!data.is_virtual(t, cnode) && (nb_generation > 0))
            {
               update(pstates, state_config, 1);
               state_config.assign(nb_state, 0); // set every value to 0
            }
         }
      }

   for(i = 0; i < nb_vomc; i++)
   {
      // create sequences
      current_nb_seq = cdeque[pstate].size();
      seq_identifier = new int[current_nb_seq];
      seq_length = new int[current_nb_seq];
      int_sequence = new int**[current_nb_seq];
      for(seq_id = 0; seq_id < current_nb_seq; seq_id++)
      {
         seq_identifier[seq_id] = seq_id + 1; // cdeque[pstate][id][pos]
         seq_length[seq_id] = cdeque[pstate][seq_id].size();
         int_sequence[seq_id] = new int*[1];
         int_sequence[seq_id][0] = new int[seq_length[seq_id]];
         for(pos = 0; pos < seq_length[seq_id]; pos++)
            int_sequence[seq_id][0][pos] = cdeque[pstate][seq_id][pos];
      }
      cseq = new sequence_analysis::Sequences(current_nb_seq, seq_identifier, seq_length, index_parameter_type, 1, STATE, int_sequence);
      seq[i] = new sequence_analysis::MarkovianSequences(*cseq);
      delete cseq;
      cseq = NULL;
      delete [] seq_identifier;
      delete [] seq_length;
      for(seq_id = 0; seq_id < current_nb_seq; seq_id++)
      {
         delete [] int_sequence[seq_id][0];
         delete [] int_sequence[seq_id];
      }
      delete [] int_sequence;
   }

   if (unordered_rank > 1)
   {
      delete [] cdeque;
      cdeque = NULL;
   }
   // frequency_distribution_computation();
}


/*****************************************************************
 *
 *  Copy operator for MarkovOutTreeReestimation class
 *
 **/

template<typename Type> void
MarkovOutTreeReestimation<Type>::copy(const MarkovOutTreeReestimation<Type>& process)
{
   unsigned int s;
   // remove() must be used before

   nb_state = process.nb_state;
   nb_children_branching = process.nb_children_branching;
   nb_factors = process.nb_factors;
   if (process.factor_values != NULL)
      factor_values = new index_array(*process.factor_values);
   else
      factor_values = NULL;

   nb_generation_values = process.nb_generation_values;
   nb_generation = process.nb_generation;
   factor_base = process.factor_base;
   unordered_rank = process.unordered_rank;
   nb_vomc = process.nb_vomc;

   if (process.seq != NULL)
   {
      seq = new sequence_analysis::MarkovianSequences*[nb_vomc];
      for(s = 0; s < nb_vomc; s++)
         if (process.seq[s] != NULL)
            seq[s] = new sequence_analysis::MarkovianSequences(*process.seq[s]);
         else
            seq[s] = NULL;
   }
   else
      seq = NULL;

   if (process.generation_reestim != NULL)
   {
      generation_reestim = new mvHistogram*[nb_generation_values];
      for(s = 0; s < nb_generation_values; s++)
         if (process.generation_reestim[s] != NULL)
            generation_reestim[s]
               = new mvHistogram(*process.generation_reestim[s]);
         else
            generation_reestim[s] = NULL;
   }
   else
      generation_reestim = NULL;

}

/*****************************************************************
 *
 *  Destructor for MarkovOutTreeReestimation class
 *
 **/

template<typename Type> void
MarkovOutTreeReestimation<Type>::remove()
{
   unsigned int s;

   if (factor_values != NULL)
   {
      delete factor_values;
      factor_values = NULL;
   }

   if (generation_reestim != NULL)
   {
      for(s = 0; s < nb_generation_values; s++)
         if (generation_reestim[s] != NULL)
         {
            delete generation_reestim[s];
            generation_reestim[s] = NULL;
         }
      delete [] generation_reestim;
      generation_reestim = NULL;
   }

   if (seq != NULL)
   {
      for(s = 0; s < nb_vomc; s++)
         if (seq[s] != NULL)
         {
            delete seq[s];
            seq[s] = NULL;
         }
      delete [] seq;
      seq = NULL;
   }
}

/*****************************************************************
 *
 *  Compute FrequencyDistribution nb_element, mean, etc.
 *  in empirical generation distributions
 *  of MarkovOutTreeReestimation class
 *
 **/

/* template<typename Type> void
MarkovOutTreeReestimation<Type>::frequency_distribution_computation() const
{
   unsigned int s;

   if (generation_reestim != NULL)
      for(s = 0; s < nb_generation_values; s++)
         if (generation_reestim[s] != NULL)
         {
            generation_reestim[s]->frequency_distribution_computation();
            generation_reestim[s]->sum_computation();
         }
}*/


#endif
