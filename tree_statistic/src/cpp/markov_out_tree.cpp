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
 *       $Id: markov_out_tree.cpp 2722 2010-12-14 14:17:56Z jbdurand $
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
 *  ---------------------------------------------------------------------------
 */

#include <iomanip>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "statiskit/core/data/marginal/multivariate.h"


#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"   // definition of DiscreteParametricModel class
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/sequence_label.h"

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

using namespace stat_tool;
using namespace sequence_analysis;
using namespace Stat_trees;

// extern int column_width(int value);
extern int stat_tool::column_width(int nb_value , const double *value , double scale);
// extern char* label(const char *file_name);

/*****************************************************************
 *
 *
 * Some notes on how to use statistic::Histograms
 *
 **/

/* Frequency of a given value
   Statiskit::Marginal::Univariate::CountTable::const_iterator it;
   it = hist.cfind(value);
   (*it).second;
*/
/* Also use the functions in mv_histogram_tools to find or list std::vector values */

/*****************************************************************
 *
 *  Constructor of class GenerationProcess using the number
 *  of states, the number of children states, of factors
 *  and parent children determining the number inb_generation
 *  of branching processes, and a collection of generation distributions.
 *
 **/

GenerationProcess::GenerationProcess(int inb_state,
                                     unsigned int inb_children_branching,
                                     unsigned int inb_factors, const index_array *ifactor_values,
                                     unsigned int inb_generation,
                                     DiscreteMultivariateDistribution * const * pgeneration)
 : nb_children_branching(inb_children_branching)
 , nb_factors(inb_factors)
 , factor_values(NULL)
 , nb_generation(inb_generation)
 , factor_base(0)
 , nb_generation_values(1)

{
   unsigned int s, j;

   nb_state = inb_state;

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
         generation = new DiscreteMultivariateDistribution*[nb_generation_values];
         for(s = 0; s < nb_generation_values; s++)
         {
            if (pgeneration[s] != NULL)
               generation[s] = pgeneration[s]->copy_ptr();
            else
               generation[s] = NULL;
         }
      }
      else
         generation = NULL;

      factor_base.resize(nb_generation, 0);
      s = 0;
      if (nb_generation == nb_children_branching + nb_factors + 1)
      {
         // parent state is a factor
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
   }
   else
      generation = NULL;
}

/*****************************************************************
 *
 *  Copy constructor for GenerationProcess class
 *
 **/

GenerationProcess::GenerationProcess(const GenerationProcess& process)
 : nb_state(0)
 , nb_children_branching(0)
 , nb_factors(0)
 , factor_values(NULL)
 , nb_generation(0)
 , nb_generation_values(0)
 , factor_base(0)
 , generation(NULL)
{ copy(process); }

/*****************************************************************
 *
 *  Destructor for GenerationProcess class
 *
 **/

GenerationProcess::~GenerationProcess()
{ remove(); }

/*****************************************************************
 *
 *  Assigment operator for GenerationProcess class
 *
 **/

GenerationProcess& GenerationProcess::operator=(const GenerationProcess& process)
{
  if (&process != this)
  {
     remove();
     copy(process);
  }

  return *this;
}

/*****************************************************************
 *
 *  Code a configuration of factors of GenerationProcess class
 *  with an unsigned int
 *
 **/

unsigned int GenerationProcess::code(const index_array& fact) const
{
   unsigned int code = 0, base = 1;
   register int j;

   assert(fact.size() == nb_generation);

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
 *  Decode a configuration of factors of GenerationProcess class
 *  coded with an unsigned int
 *
 **/

GenerationProcess::index_array GenerationProcess::decode(unsigned int code) const
{
   unsigned int base;
   register int j;
   index_array decodev(nb_generation, 0);

   assert(code < nb_generation_values);

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
 *  Return the number of ordered children involved
 *  in the branching processes as factors (due to dependencies
 *  regarding these states, excluding parent state)
 *  in GenerationProcess
 *
 **/

unsigned int GenerationProcess::get_nb_children_branching() const
{ return nb_children_branching; }

/*****************************************************************
 *
 *  Return the number of other factors involved
 *  in the branching process in GenerationProcess
 *
 **/

unsigned int GenerationProcess::get_nb_factors() const
{ return nb_factors; }

/*****************************************************************
 *
 *  Return the number of possible values for each factor
 *  in GenerationProcess
 *
 **/

GenerationProcess::index_array GenerationProcess::get_factor_values() const
{
   index_array res(0);

   if (factor_values != NULL)
      res = *factor_values;

   return res;
}

/*****************************************************************
 *
 *  Return the number of generation processes in GenerationProcess
 *  (combining all possible values of every factor)
 *
 **/

unsigned int GenerationProcess::get_nb_generation_values() const
{ return nb_generation_values; }

/*****************************************************************
 *
 *  Return the number of factors for generation processes
 *  (including parent state) in GenerationProcess
 *
 **/

unsigned int GenerationProcess::get_nb_generation() const
{ return nb_generation; }

/*****************************************************************
 *
 *  Matplotlib output for GenerationProcess class
 *  using potential data
 *
 **/

MultiPlotSet* GenerationProcess::get_plotable(Statiskit::Marginal::Multivariate::CountTable* const * histo) const
{
   // plot_set[state][view]
   int index = 0;
   const int nb_plot_set = nb_plot_set_computation();

   MultiPlotSet *plot_set = NULL;

   plot_set = new MultiPlotSet(nb_plot_set, 1);
   plot_set->border = "15 lw 0";
   plotable_write(*plot_set, index, histo);

   return plot_set;
}

/*
MultiPlotSet* GenerationProcess::get_plotable(Statiskit::Marginal::Multivariate::CountTable* const * histo) const
{
   // plot_set[state][view]
   register int j, k, p;
   unsigned int i, var;
   const int nb_plot_set = nb_generation_values * nb_state;
   ostringstream title, legend;
   index_array *decodev = NULL; // configuration of factors
   MultiPlotSet *generation_plot = NULL, *plot_set = NULL;

   plot_set = new MultiPlotSet(nb_plot_set);
   plot_set->border = "15 lw 0";

   title.str("");
   title << STAT_TREES_label[TREESTATL_GENERATION_PROCESS] << ": ";
   if (histo != NULL)
      title << " " << STAT_label[STATL_FIT];

   plot_set->title = title.str();

   decodev = new index_array(nb_generation, 0);

   i = 0;
   while (i < nb_plot_set)
   {
      title.str("");
      k = i / nb_state; // generation_value (i.e. configuration of factors)
      // (*plot_set)[i].resize(nb_state);
      j = i % nb_state; // unordered children state

      // should be called every nb_state only
      if (j == 0)
         if ((histo != NULL) && (histo[k] != NULL))
            generation_plot = generation[k]->get_plotable(histo[k]);
         else
            generation_plot = generation[k]->get_plotable();
      (*plot_set)[i].resize((*generation_plot)[j].size());
      for(p = 0; p < (*generation_plot)[j].size(); p++)
         (*plot_set)[i][p] = (*generation_plot)[j][p];
      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << j << " "
             << STAT_TREES_label[TREESTATL_GENERATION]
             << " " << STAT_label[STATL_DISTRIBUTION] << ": ";
      title << legend.str();
      // (*plot_set)[i].legend = legend.str();

      // recode k using factor_base
      *decodev = decode(k);
      j = 0;
      if (nb_generation == nb_children_branching + nb_factors + 1)
      {
         // parent state is a factor
         title << STAT_TREES_label[TREESTATL_PARENT] << " " << STAT_label[STATL_STATE] << " "
               << (*decodev)[j] << " ";
          j++;
       }
       k = 0;
       while (k < nb_children_branching)
       {
          title << STAT_TREES_word[TREESTATW_CHILD] << " " << k + 1 << " "
                << STAT_label[STATL_STATE] << " " << (*decodev)[j] << " ";
          k++;
          j++;
       }
       k = 0;
       while (k < nb_factors)
       {
          title << STAT_TREES_word[TREESTATW_FACTOR] << " " << k << " "
                << STAT_trees_word[STATW_VALUE] << " " << (*decodev)[j] << " ";
          k++;
          j++;
       }
       // title << STAT_label[STATL_DISTRIBUTION];
       (*plot_set)[i].title = title.str();
       if (j == 0)
          if ((histo != NULL) && (histo[k] != NULL))
          {
             delete generation_plot;
             generation_plot = NULL;
          }
       i++;
    } // end while i


   delete decodev;
   decodev = NULL;

   return plot_set;
}
*/

/*****************************************************************
 *
 *  Return the number of views (i.e. the size) in Matplotlib output
 *  for GenerationProcess class using potential data
 *
 **/

 unsigned int GenerationProcess::nb_plot_set_computation(mvHistogram* const * histo) const
{ return nb_generation_values * nb_state; }

/*****************************************************************
 *
 *  Write views into a MultiPlotSet object
 *  for GenerationProcess class using instance to update,
 *  current index in this instance and potential data
 *
 **/

void GenerationProcess::plotable_write(MultiPlotSet &plot_set, int &index, mvHistogram* const * histo) const
{
   // plot_set[state][view]
   register int j, k, p;
   unsigned int i, var;
   const int nb_plot_set = nb_plot_set_computation();
   ostringstream title, legend;
   index_array *decodev = NULL; // configuration of factors
   MultiPlotSet *generation_plot = NULL;

   title.str("");
   title << STAT_TREES_label[TREESTATL_GENERATION_PROCESS] << ": ";
   if (histo != NULL)
      title << " " << STAT_label[STATL_FIT];

   plot_set.title = title.str();

   decodev = new index_array(nb_generation, 0);

   i = 0;
   while (i < nb_plot_set)
   {
      plot_set.variable[index] = 0; // I_DEFAULT ?
      plot_set.viewpoint[index] = 0; // I_DEFAULT ? STATE_PROCESS ?

      title.str("");
      k = i / nb_state; // generation_value (i.e. configuration of factors)
      // (*plot_set)[i].resize(nb_state);
      j = i % nb_state; // unordered children state

      // should be called every nb_state only
      if (j == 0)
         if ((histo != NULL) && (histo[k] != NULL))
            generation_plot = generation[k]->get_plotable(histo[k]);
         else
            generation_plot = generation[k]->get_plotable();
      plot_set[index].resize((*generation_plot)[j].size());
      for(p = 0; p < (*generation_plot)[j].size(); p++)
         plot_set[index][p] = (*generation_plot)[j][p];
      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << j << " "
             << STAT_TREES_label[TREESTATL_GENERATION]
             << " " << STAT_label[STATL_DISTRIBUTION] << ": ";
      title << legend.str();
      // (*plot_set)[i].legend = legend.str();

      // recode k using factor_base
      *decodev = decode(k);
      j = 0;
      if (nb_generation == nb_children_branching + nb_factors + 1)
      {
         // parent state is a factor
         title << STAT_TREES_label[TREESTATL_PARENT] << " " << STAT_label[STATL_STATE] << " "
               << (*decodev)[j] << " ";
          j++;
       }
       k = 0;
       while (k < nb_children_branching)
       {
          title << STAT_TREES_word[TREESTATW_CHILD] << " " << k + 1 << " "
                << STAT_label[STATL_STATE] << " " << (*decodev)[j] << " ";
          k++;
          j++;
       }
       k = 0;
       while (k < nb_factors)
       {
          title << STAT_TREES_word[TREESTATW_FACTOR] << " " << k << " "
                << STAT_TREES_word[TREESTATW_VALUE] << " " << (*decodev)[j] << " ";
          k++;
          j++;
       }
       // title << STAT_label[STATL_DISTRIBUTION];
       plot_set[index].title = title.str();
       if (j == 0)
          if ((histo != NULL) && (histo[k] != NULL))
          {
             delete generation_plot;
             generation_plot = NULL;
          }
       i++;
       index++;
    } // end while i


   delete decodev;
   decodev = NULL;

}

/*****************************************************************
 *
 *  Copy operator for GenerationProcess class
 *
 **/

void GenerationProcess::copy(const GenerationProcess& process)
{
   unsigned int s;

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

   if (nb_generation_values > 0)
   {
      generation = new DiscreteMultivariateDistribution*[nb_generation_values];

      if (process.generation != NULL)
         for(s = 0; s < nb_generation_values; s++)
         {
            if (process.generation[s] != NULL)
               generation[s] = process.generation[s]->copy_ptr();
            else
               generation[s] = NULL;
         }
   }
   else
      generation = NULL;
}

/*****************************************************************
 *
 *  Destructor for GenerationProcess class
 *
 **/

void GenerationProcess::remove()
{
   unsigned int s;

   if (factor_values != NULL)
   {
      delete factor_values;
      factor_values = NULL;
   }

   if (generation != NULL)
   {
      for(s = 0; s < nb_generation_values; s++)
         if (generation[s] != NULL)
         {
            delete generation[s];
            generation[s] = NULL;
         }
      delete [] generation;
      generation = NULL;
   }
}

/*****************************************************************
 *
 *  Print a GenerationProcess
 *  using an output stream,
 *  empirical observation distributions,
 *  and flags on the level of detail and on the file use
 *
 **/

std::ostream& GenerationProcess::ascii_print(std::ostream &os,
                                             mvHistogram **empirical_observation,
                                             bool exhaustive, bool file_flag) const
{
   register int i, j, k, code, base;
   index_array *decodev = NULL; // configuration of factors
   // int buff , width[2];
   double *pmass;

   if (generation != NULL)
   {
      decodev = new index_array(nb_generation, 0);
      // print generation processes
      os << nb_generation << " "
         << STAT_TREES_word[TREESTATW_FACTOR] << " ";
      if (nb_generation > 1)
         os << STAT_TREES_word[TREESTATW_GENERATION_PROCESSES];
      else
         os << STAT_TREES_word[TREESTATW_GENERATION_PROCESS];
      os << endl << endl;

      if (nb_generation > 0)
      {
         os << nb_children_branching << " ";
         if (nb_children_branching > 1)
            os << STAT_TREES_word[TREESTATW_CHILDREN] << " ";
         else
            os << STAT_TREES_word[TREESTATW_CHILD] << " ";

         if (nb_children_branching > 1)
            os << STAT_TREES_word[TREESTATW_GENERATION_PROCESSES];
         else
            os << STAT_TREES_word[TREESTATW_GENERATION_PROCESS];
         os << endl << endl;
         os << nb_factors << " ";
         if (nb_factors > 1)
            os << STAT_TREES_word[TREESTATW_EXTERNAL_FACTORS] << " ";
         else
            os << STAT_TREES_word[TREESTATW_EXTERNAL_FACTOR] << " ";
         if (nb_factors > 1)
            os << STAT_TREES_word[TREESTATW_GENERATION_PROCESSES];
         else
            os << STAT_TREES_word[TREESTATW_GENERATION_PROCESS];
         os << endl << endl;

         // print nb_values of factors
         for(i = 0; i < nb_factors; i++)
         {
            os << STAT_TREES_word[TREESTATW_FACTOR] << " " << i+1
               << " : " << (*factor_values)[i] << " "
               << STAT_label[STATL_VALUES] << endl;
         }
         // compute the number of values for each generation process
         // parent state first, then children states, then external factors
         for(i = 0; i < nb_generation_values; i++)
         {
            // recode i using factor_base
            *decodev = decode(i);
            j = 0;
            if (nb_generation == nb_children_branching + nb_factors + 1)
            {
               // parent state is a factor
               os << STAT_TREES_word[TREESTATW_PARENT_STATE] << " " << (*decodev)[j] << " ";
               j++;
            }
            k = 0;
            while (k < nb_children_branching)
            {
               os << STAT_TREES_word[TREESTATW_CHILD] << " " << k + 1 << " "
                  << STAT_word[STATW_STATE] << " " << (*decodev)[j] << " ";
               k++;
               j++;
            }
            k = 0;
            while (k < nb_factors)
            {
               os << STAT_TREES_word[TREESTATW_FACTOR] << " " << k << " "
                  << STAT_TREES_word[TREESTATW_VALUE] << " " << (*decodev)[j] << " ";
               k++;
               j++;
            }
            os << STAT_TREES_word[TREESTATW_GENERATION_DISTRIBUTION] << endl;
            if (generation[i] != NULL)
            {
               if (empirical_observation == NULL)
                  generation[i]->ascii_print(os, file_flag, exhaustive, false, NULL);
               else
                  generation[i]->ascii_print(os, file_flag, exhaustive, false, empirical_observation[i]);
            }
            os << endl;
         } // end for
      } // end if (nb_generation > 0)
      delete decodev;
      decodev = NULL;
   }

   return os;
}


std::ostream& GenerationProcess::spreadsheet_print(std::ostream &os, int process,
                                                   mvHistogram **empirical_observation) const
{}

/*****************************************************************
 *
 *  Gnuplot output of a GenerationProcess, using
 *  a file prefix, the title of the output figures,
 *  the identifier of the considered process
 *  and the empirical observation distributions
 *
 **/

bool GenerationProcess::plot_print(const char * prefix, const char * title, int process,
                                   mvHistogram **empirical_observation) const
{
   bool status = false;

   return status;
}

/*****************************************************************
 *
 *  State permutation for GenerationProcess class
 *  (checking validity of permutation)
 *
 **/

void GenerationProcess::state_permutation(StatError& error, int* perm) const
{
   register int i, j;
   bool status = true;
   // each element of check_perm must be used exactly once
   bool * const check_perm = new bool[nb_state];

   // check permutation
   error.init();

   for (i= 0; i < nb_state; i++)
      check_perm[i] = false; // indicates if ith element was used in perm

   for (i= 0; i < nb_state; i++)
   {
      if (check_perm[perm[i]])
         status = false;
      else
         check_perm[perm[i]] = true;
   }

   if (!status)
      error.update(STAT_TREES_error[TREESTATR_NO_PERMUTATION]);
   else
   {
      // permutation of the labels of the generation distributions
      DiscreteMultivariateDistribution **pgeneration = new DiscreteMultivariateDistribution*[nb_state];

      for(i = 0; i < nb_state; i++)
         pgeneration[perm[i]] = generation[i];

      for(i = 0; i < nb_state; i++)
         generation[i] = pgeneration[i];

      delete [] pgeneration;
      pgeneration = NULL;

      // permutation of states withing generation distributions

      for(i = 0; i < nb_state; i++)
         generation[i]->variable_permutation(error, perm, false);

   }
}

GenerationProcess* Stat_trees::generation_process_parsing(StatError &error,
                                                          std::ifstream &in_file,
                                                          int &line,
                                                          int nb_state,
                                                          double cumul_threshold)
{
   RWLocaleSnapshot locale("en");
   RWCString buffer, token;
   size_t position;
   bool status = true, lstatus;
   register int i, j, k, l, m, n, o;
   // int line_save;
   int nb_children_branching_read = -1,
                nb_factors_read = -1, nb_generation_read = -1,
                nb_expected_i, fact;
   unsigned int nb_children_branching = 0,
                nb_factors = 0, nb_generation = 0, nb_generation_values = 1;
   long value;
   GenerationProcess::index_array decodev;
   GenerationProcess *gprocess = NULL, *gprocess_copy = NULL;
   MarkovOutTree *markov = NULL;
   GenerationProcess::index_array *factor_values = NULL;
   DiscreteMultivariateDistribution **pgeneration = NULL;

   while (buffer.readLine(in_file, false))
   {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS)
         buffer.remove(position);

      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull()))
      {
         switch (i)
         {
            // read the number of factors involved in generation process
            case 0 :
            {
               lstatus = locale.stringToNum(token, &value);
               if (lstatus)
               {
                  if (value < 0)
                     lstatus = false;
                  else
                  {
                      nb_generation_read = value;
                      nb_generation = (unsigned int)(nb_generation_read);
                  }
               }
               if (!lstatus)
               {
                  status = false;
                  error.update(STAT_TREES_parsing[TREESTATP_NB_GENERATION_PROCESS], line, i + 1);
               }
               break;
            }

            // test for keyword FACTOR

            case 1 :
            {
               if (token != STAT_TREES_word[TREESTATW_FACTOR])
               {
                  status = false;
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_FACTOR], line, i + 1);
               }
               break;
            }

            case 2 :
            {
               if (((nb_generation > 1) && (token != STAT_TREES_word[TREESTATW_GENERATION_PROCESSES]))
                   || ((nb_generation <= 1) && (token != STAT_TREES_word[TREESTATW_GENERATION_PROCESS])))
               {
                  status = false;
                  if (nb_generation > 1)
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_GENERATION_PROCESSES], line, i + 1);
                  else
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_GENERATION_PROCESS], line, i + 1);
               }
               break;
            }
         }
         i++;
      }

      if (i > 0)
      {
         if (i != 3)
         {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
         }
         break;
      }
   }

   if (nb_generation_read == -1)
   {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT], line);
   }

   // read CHILD(REN) GENERATION_PROCESS
   while (buffer.readLine(in_file, false))
   {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS)
         buffer.remove(position);

      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull()))
      {
         switch (i)
         {
            // read the number of factors associated with ordered children
            case 0 :
            {
               lstatus = locale.stringToNum(token, &value);
               if (lstatus)
               {
                  if (value < 0)
                     lstatus = false;
                  else
                  {
                      nb_children_branching_read = value;
                      nb_children_branching = (unsigned int)(nb_children_branching_read);
                  }
               }
               if (!lstatus)
               {
                  status = false;
                  error.update(STAT_TREES_parsing[TREESTATP_NB_CHILDREN_BRANCHING], line, i + 1);
               }
               break;
            }

            // test for keyword CHILD || CHILDREN

            case 1 :
            {
               if (((nb_children_branching > 1) && (token != STAT_TREES_word[TREESTATW_CHILDREN]))
                   || ((nb_children_branching <= 1) && (token != STAT_TREES_word[TREESTATW_CHILD])))
               {
                  status = false;
                  if (nb_children_branching > 1)
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_CHILDREN], line, i + 1);
                  else
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_CHILD], line, i + 1);
               }
               break;
            }

            case 2 :
            {
               if (((nb_children_branching > 1) && (token != STAT_TREES_word[TREESTATW_GENERATION_PROCESSES]))
                   || ((nb_children_branching <= 1) && (token != STAT_TREES_word[TREESTATW_GENERATION_PROCESS])))
               {
                  status = false;
                  if (nb_children_branching > 1)
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_GENERATION_PROCESSES], line, i + 1);
                  else
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_GENERATION_PROCESSES], line, i + 1);
               }
               break;
            }
         }
         i++;
      }

      if (i > 0)
      {
         if (i != 3)
         {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
         }
         break;
      }
   }

   if (nb_children_branching == -1)
   {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT], line);
   }

   // read EXTERNAL_FACTOR(S) GENERATION_PROCESS
   while (buffer.readLine(in_file, false))
   {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS)
         buffer.remove(position);

      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull()))
      {
         switch (i)
         {
            // read the number of external factors
            case 0 :
            {
               lstatus = locale.stringToNum(token, &value);
               if (lstatus)
               {
                  if (value < 0)
                     lstatus = false;
                  else
                  {
                      nb_factors_read = value;
                      nb_factors = (unsigned int)(nb_factors_read);
                  }
               }
               if (!lstatus)
               {
                  status = false;
                  error.update(STAT_TREES_parsing[TREESTATP_NB_FACTORS], line, i + 1);
               }
               break;
            }

            // test for keyword EXTERNAL_FACTOR(S)

            case 1 :
            {
               if (((nb_factors > 1) && (token != STAT_TREES_word[TREESTATW_EXTERNAL_FACTORS]))
                   || ((nb_factors <= 1) && (token != STAT_TREES_word[TREESTATW_EXTERNAL_FACTOR])))
               {
                  status = false;
                  if (nb_factors > 1)
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_EXTERNAL_FACTORS], line, i + 1);
                  else
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_EXTERNAL_FACTOR], line, i + 1);
               }
               break;
            }

            case 2 :
            {
               if (((nb_factors > 1) && (token != STAT_TREES_word[TREESTATW_GENERATION_PROCESSES]))
                   || ((nb_factors <= 1) && (token != STAT_TREES_word[TREESTATW_GENERATION_PROCESS])))
               {
                  status = false;
                  if (nb_factors > 1)
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_GENERATION_PROCESSES], line, i + 1);
                  else
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_GENERATION_PROCESSES], line, i + 1);
               }
               break;
            }
         }
         i++;
      }

      if (i > 0)
      {
         if (i != 3)
         {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
         }
         break;
      }
   }

   if (nb_factors == -1)
   {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT], line);
   }

   if (nb_factors > 0)
   // parse values of factors
   {
      factor_values = new GenerationProcess::index_array(nb_factors);
      // read number of values of each factor
      for(fact=0; fact < nb_factors; fact++)
      {
         // read FACTOR fact : nb values
         while (buffer.readLine(in_file, false))
         {
            line++;

#           ifdef DEBUG
            cout << line << "  " << buffer << endl;
#           endif

            position = buffer.first('#');
            if (position != RW_NPOS)
               buffer.remove(position);

            i = 0;

            RWCTokenizer next(buffer);

            while (!((token = next()).isNull()))
            {
               switch (i)
               {
                  // read the number of external factors
                  case 0 :
                  {
                     // read FACTOR
                     if (token != STAT_TREES_word[TREESTATW_FACTOR])
                     {
                        status = false;
                        error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_FACTOR], line, i + 1);
                     }
                     break;
                  }
                  case 1 :
                  {
                     // read factor number
                     lstatus = locale.stringToNum(token, &value);
                     if ((lstatus) && (value != fact + 1))
                        lstatus = false;
                     if (!lstatus)
                     {
                        status = false;
                        error.update(STAT_TREES_parsing[TREESTATP_NB_GENERATION_PROCESS], line, i + 1);
                     }
                     break;
                  }

                  case 2 :
                  {
                     // test separator ":"
                     if (token != ":")
                     {
                        status = false;
                        error.update(STAT_parsing[STATP_SEPARATOR], line, i + 1);
                     }
                     break;
                  }

                  case 3 :
                  {
                     // read number of values for factor
                     lstatus = locale.stringToNum(token, &value);
                     if ((lstatus) && (value <= 1))
                        lstatus = false;
                     else
                        (*factor_values)[fact] = value;
                     if (!lstatus)
                     {
                        status = false;
                        error.update(STAT_TREES_parsing[TREESTATP_NB_FACTOR_VALUES], line, i + 1);
                     }
                     break;
                  }

                  case 4 :
                  {
                     // read values
                     if (token != STAT_label[STATL_VALUES])
                     {
                        status = false;
                        error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_label[STATL_VALUES], line, i + 1);
                     }
                  }
               }
               i++;
            }

            if (i > 0)
            {
               if (i != 5)
               {
                  status = false;
                  error.update(STAT_parsing[STATP_FORMAT] , line);
               }
               break;
            }
         }
      } // end for fact
   }

   if (status)
   {
      // read GENERATION_DISTRIBUTIONS

      gprocess_copy = new GenerationProcess(nb_state, nb_children_branching,
                                            nb_factors, factor_values, nb_generation,
                                            pgeneration);

      nb_generation_values = gprocess_copy->nb_generation_values;
      pgeneration = new DiscreteMultivariateDistribution*[nb_generation_values];
      decodev.resize(nb_generation, 0);

      for(j = 0; j < nb_generation_values; j++)
      {
         // recode j using factor_base
         decodev = gprocess_copy->decode(j);
         k = 0; // index in decode v
         l = 0; // index in child
         m = 0; // index in state (or factor ?)

         if (nb_generation == nb_children_branching + nb_factors + 1)
            // parent state is a factor
            n = 0; // index in parent
         else
            n = 1; // skip reading parent state values

         o = 0; // index in external factor

         while (buffer.readLine(in_file, false))
         {
            line++;

   #        ifdef DEBUG
            cout << line << "  " << buffer << endl;
   #        endif

            position = buffer.first('#');
            if (position != RW_NPOS)
               buffer.remove(position);

            i = 0;

            RWCTokenizer next(buffer);

            // read words in current line
            while (!((token = next()).isNull()))
            {
               switch (i % 2)
               {
                  // read keyword
                  case 0 :
                  {
                     if (n == 0) // parent state is a factor
                     {
                        if (token != STAT_TREES_word[TREESTATW_PARENT_STATE])
                        {
                           lstatus = false;
                           error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_PARENT_STATE], line, i + 1);
                        }
                     }
                     else
                     {
                        // read states of ordered brother vertices
                        if (l < nb_children_branching)
                        {
                           if (m == 0)
                           {
                              if (token != STAT_TREES_word[TREESTATW_CHILD])
                              {
                                 lstatus = false;
                                 error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_CHILD], line, i + 1);
                              }
                           }
                           else
                           {
                              if (token != STAT_word[STATW_STATE])
                              {
                                 lstatus = false;
                                 error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_word[STATW_STATE], line, i + 1);
                              }
                           }
                        }
                        else
                        {
                           // read Id of factor
                           if (o < nb_factors)
                           {
                              if (m == 0)
                              {
                                 if (token != STAT_TREES_word[TREESTATW_FACTOR])
                                 {
                                    lstatus = false;
                                    error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_FACTOR], line, i + 1);
                                 }
                              }
                              else
                              {
                                 if (token != STAT_TREES_word[TREESTATW_VALUE])
                                 {
                                    lstatus = false;
                                    error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_VALUE], line, i + 1);
                                 }
                              }
                           }
                           else
                           {
                              if (token != STAT_TREES_word[TREESTATW_GENERATION_DISTRIBUTION])
                              {
                                 lstatus = false;
                                 error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_GENERATION_DISTRIBUTION], line, i + 1);
                              }
                           }
                        }
                     }
                     break;
                  }

                  // read value
                  case 1 :
                  {
                     lstatus = locale.stringToNum(token, &value);
                     if (lstatus)
                     {
                        if (n == 0)
                        {
                           // read parent state
                           if (value != decodev[k])
                              lstatus = false;
                           k++;
                           n++;
                        }
                        else
                        {
                           // read child state
                           if (l < nb_children_branching)
                           {
                              if (m == 0)
                              {
                                 // read child number
                                 if (value != l + 1)
                                    lstatus = false;
                              }
                              if (m == 1)
                              {
                                 // read child state value
                                 if (value != decodev[k])
                                    lstatus = false;
                                 k++;
                                 l++; // read next CHILD
                              }
                              m = 1-m; // read STATE if CHILD was read and vice-versa
                           }
                           else
                              if (o < nb_factors)
                              {
                                 if (m == 0)
                                 {
                                    // read factor ID
                                    if (value != o)
                                       lstatus = false;
                                 }
                                 if (m == 1)
                                 {
                                    // read factor value
                                    if (value != decodev[k])
                                       lstatus = false;
                                    k++;
                                    o++; // read next FACTOR
                                 }
                                 m = 1-m; // read VALUE if FACTOR was read and vice-versa
                              }

                        }
                     }
                     if (!lstatus)
                     {
                        status = false;
                        if (n == 0)
                           error.update(STAT_TREES_word[TREESTATW_PARENT_STATE], line, i + 1);
                        if (l <= nb_children_branching)
                        {
                           if (m == 0)
                              error.update(STAT_TREES_word[TREESTATW_CHILD], line, i + 1);
                           else
                              error.update(STAT_word[STATW_STATE], line, i + 1);
                        }
                        if (o <= nb_factors)
                        {
                           if (m == 0)
                              error.update(STAT_TREES_word[TREESTATW_FACTOR], line, i + 1);
                           else
                              error.update(STAT_TREES_word[TREESTATW_VALUE], line, i + 1);
                        }
                     }
                     break;
                  }
               }
               i++;
            }

            if (i > 0)
            {
               nb_expected_i = 2 * ((nb_generation - nb_children_branching - nb_factors)
                                    + 2 * nb_children_branching
                                    + 2 * nb_factors) + 1;
               if (i != nb_expected_i)
               {
                  status = false;
                  error.update(STAT_parsing[STATP_FORMAT] , line);
               }
               break;
            }
         }
         pgeneration[j] = NULL;
         pgeneration[j] = discrete_multivariate_parsing(error, in_file, line,
                                                        MCOMPOUND_MULTINOMIAL,
                                                        cumul_threshold);
         if (pgeneration[j] == NULL)
            status = false;
         else
            // check number of states
            if (pgeneration[j]->nb_variable != nb_state)
            {
               status = false;
               error.update(STAT_parsing[STATP_NB_STATE], line);
            }
      } // end for
   }

   if (status)
   {
      gprocess_copy->generation = pgeneration;
      gprocess = new GenerationProcess(*gprocess_copy);
      delete gprocess_copy;
      gprocess_copy = NULL;
   }

   if (factor_values != NULL)
   {
      delete factor_values;
      factor_values = NULL;
   }

   return gprocess;
}

/*****************************************************************
 *
 *  Default constructor of MarkovOutTree class
 *
 **/

MarkovOutTree::MarkovOutTree()
 : Chain('o', 0, false)
 , markov_data(NULL)
 , nb_output_process(0)
 , nonparametric_process(NULL)
 , discrete_parametric_process(NULL)
 , continuous_parametric_process(NULL)
 , unordered_rank(0)
 , nb_vomc(0)
 , vomc_models(NULL)
 , generation_distributions(NULL)
{}

/*****************************************************************
 *
 *  Constructor of MarkovOutTree class using the number of states,
 *  the number and values of observation processes, the number of
 *  ordered children, the number and values of Markov models
 *  for ordered children, and the branching processes
 *  The CategoricalTreeProcess npprocess should not include
 *  the state process, and all observation processes are indexed
 *  starting from (observed variable) 0
 *
 **/


MarkovOutTree::MarkovOutTree(int inb_state, int inb_output_process, CategoricalProcess * const * npprocess,
                             DiscreteParametricProcess * const * dpprocess, ContinuousParametricProcess *const * cprocess,
                             unsigned int iunordered_rank, unsigned int inb_vomc,
                             VariableOrderMarkov * const *  pvomc_models,
                             const GenerationProcess *pgeneration)
 : Chain('o', inb_state, true)
 , markov_data(NULL)
 , nb_output_process(inb_output_process)
 , nonparametric_process(NULL)
 , discrete_parametric_process(NULL)
 , continuous_parametric_process(NULL)
 , unordered_rank(iunordered_rank)
 , nb_vomc(inb_vomc)
 , vomc_models(NULL)
 , generation_distributions(NULL)
{
   unsigned int i;

   observation_copy(inb_output_process, npprocess, dpprocess, cprocess);

   if (nb_vomc > 0)
   {
      vomc_models = new VariableOrderMarkov*[nb_vomc];
      for(i = 0; i < nb_vomc; i++)
         vomc_models[i] = new VariableOrderMarkov(*pvomc_models[i]);
   }

   if (nb_vomc == 1)
      for(i = 0; i < nb_state; i++)
         initial[i] = vomc_models[0]->initial[i];

   if (pgeneration != NULL)
      generation_distributions = new GenerationProcess(*pgeneration);
}

/*****************************************************************
 *
 *  Constructor of MarkovOutTree class using a Markov chain
 *  (only the number of states and initial probabilities may be used),
 *  the number and values of observation processes, the number of
 *  ordered children, the number and values of Markov models
 *  for ordered children, and the branching processes
 *  The CategoricalTreeProcess npprocess should not include
 *  the state process, and all observation processes are indexed
 *  starting from (observed variable) 0
 *
 **/


MarkovOutTree::MarkovOutTree(const Chain& markov, int inb_output_process, CategoricalProcess * const * npprocess,
                             DiscreteParametricProcess * const * dpprocess, ContinuousParametricProcess *const * cprocess,
                             unsigned int iunordered_rank, unsigned int inb_vomc,
                             VariableOrderMarkov * const *  pvomc_models,
                             const GenerationProcess *pgeneration)
 : Chain(markov)
 , markov_data(NULL)
 , nb_output_process(inb_output_process)
 , nonparametric_process(NULL)
 , discrete_parametric_process(NULL)
 , continuous_parametric_process(NULL)
 , unordered_rank(iunordered_rank)
 , nb_vomc(inb_vomc)
 , vomc_models(NULL)
 , generation_distributions(NULL)
{
   unsigned int i;

   observation_copy(inb_output_process, npprocess, dpprocess, cprocess);

   if (nb_vomc > 0)
   {
      vomc_models = new VariableOrderMarkov*[nb_vomc];
      for(i = 0; i < nb_vomc; i++)
         vomc_models[i] = new VariableOrderMarkov(*pvomc_models[i]);
   }

   if (pgeneration != NULL)
      generation_distributions = new GenerationProcess(*pgeneration);
}

/*****************************************************************
 *
 *  Create a MarkovOutTree copying every field from a given MarkovOutTree
 *  except for the number and values of observation processes
 *  (given as argument)
 *
 **/

MarkovOutTree::MarkovOutTree(const MarkovOutTree * markov, int inb_output_process,
                             CategoricalProcess * const * nonparametric_observation,
                             DiscreteParametricProcess * const * discrete_parametric_observation,
                             ContinuousParametricProcess * const * continuous_parametric_observation,
                             int size)
 : StatInterface()
 , Chain(*markov)
 , markov_data(NULL)
 , nb_output_process(0)
 , nonparametric_process(NULL)
 , discrete_parametric_process(NULL)
 , continuous_parametric_process(NULL)
 , unordered_rank(0)
 , nb_vomc(0)
 , vomc_models(NULL)
 , generation_distributions(NULL)
{
   if (markov != NULL)
      state_copy(*markov);

   observation_copy(inb_output_process, nonparametric_observation,
                    discrete_parametric_observation, continuous_parametric_observation);

   characteristic_computation(size, true);
}

/*****************************************************************
 *
 *  Copy constructor of MarkovOutTree class using a flag
 *  on copying markov_data
 *
 **/

MarkovOutTree::MarkovOutTree(const MarkovOutTree& markov, bool data_flag)
 : StatInterface()
 , Chain(markov)
 , markov_data(NULL)
 , nb_output_process(0)
 , nonparametric_process(NULL)
 , discrete_parametric_process(NULL)
 , continuous_parametric_process(NULL)
 , unordered_rank(0)
 , nb_vomc(0)
 , vomc_models(NULL)
 , generation_distributions(NULL)
{ copy(markov, data_flag); }

/*****************************************************************
 *
 *  Destructor for MarkovOutTree class
 *
 **/

MarkovOutTree::~MarkovOutTree()
{ remove(); }

/*****************************************************************
 *
 *  Return the data part of a MarkovOutTree, keeping a reference
 *  on current object instance
 *
 **/

MarkovOutTreeData* MarkovOutTree::extract_data(StatError& error) const
{
   bool status= true;
   MarkovOutTreeData *res = NULL;

   error.init();

   if (markov_data == NULL)
   {
      status= false;
      error.update(STAT_error[STATR_NO_DATA]);
   }

   if (nb_output_process != (markov_data->_nb_integral
        + markov_data->_nb_float))
   {
      status= false;
      error.update(STAT_TREES_error[TREESTATR_STATE_TREES]);
   }

   if (status)
   {
      res = new MarkovOutTreeData(*markov_data, true, true);
      res->markov = new MarkovOutTree(*this, false);
   }

   return res;
}

/*****************************************************************
 *
 *  Return the number of states of a MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_state() const
{ return nb_state; }

/*****************************************************************
 *
 *  Return the number of output processes of a MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_output_process() const
{ return nb_output_process; }

/*****************************************************************
 *
 *  Return the number of nonparametric output processes
 *  of a MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_nonparametric_output_process() const
{
   unsigned int var, res = 0;

   if (nonparametric_process == NULL)
      return 0;
   else
      for(var = 1; var <= nb_output_process; var++)
         if (nonparametric_process[var] != NULL)
            res++;
   return res;
}

/*****************************************************************
 *
 *  Return the number of discrete parametric output processes
 *  of a MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_discrete_parametric_output_process() const
{
   unsigned int var, res = 0;

   if (discrete_parametric_process == NULL)
      return 0;
   else
      for(var = 1; var <= nb_output_process; var++)
         if (discrete_parametric_process[var] != NULL)
            res++;
   return res;
}

/*****************************************************************
 *
 *  Return the number of continuous (parametric) output processes
 *  of a MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_continuous_output_process() const
{
   unsigned int var, res = 0;

   if (continuous_parametric_process == NULL)
      return 0;
   else
      for(var = 1; var <= nb_output_process; var++)
         if (continuous_parametric_process[var] != NULL)
            res++;
   return res;
}

/*****************************************************************
 *
 *  Print a MarkovOutTree on a single line
 *  using an output stream
 *
 **/

ostream& MarkovOutTree::line_write(ostream& os) const
{
   os << nb_state << " " << STAT_word[STATW_STATES];

   return os;
}

/*****************************************************************
 *
 *  Print a MarkovOutTree
 *  using an output stream and a flag on the level of detail
 *
 **/

ostream& MarkovOutTree::ascii_write(ostream& os, bool exhaustive) const
{ return ascii_write(os, markov_data, exhaustive, false); }

/*****************************************************************
 *
 *  Print a MarkovOutTree into a file
 *  using a StatError object, the path
 *  and a flag on the level of detail
 *
 **/

bool MarkovOutTree::ascii_write(StatError& error, const char * path,
                                bool exhaustive) const
{
   bool status;
   ofstream out_file(path);

   error.init();

   if (!out_file)
   {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
   }
   else
   {
      status = true;
      ascii_write(out_file, markov_data, exhaustive, true);
   }
   return status;
}

/*****************************************************************
 *
 *  Print a MarkovOutTree in a spreadsheet fashion
 *  using a StatError object and the path
 *
 **/

bool MarkovOutTree::spreadsheet_write(StatError& error,
                                      const char * path) const
{
   bool status;
   ofstream out_file(path);

   error.init();

   if (!out_file)
   {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
   }

   else
   {
      status = true;
      spreadsheet_write(out_file, markov_data);
   }
   return status;
}

/*****************************************************************
 *
 *  Gnuplot output for MarkovOutTree class
 *  using a StatError object, a prefix for the files
 *  and the title of output figures
 *
 **/

bool MarkovOutTree::plot_write(StatError& error,
                               const char * prefix,
                               const char * title) const
{
   bool status = plot_write(prefix, title, markov_data);

   error.init();

   if (!status)
     error.update(STAT_error[STATR_FILE_PREFIX]);

   return status;
}

/*****************************************************************
 *
 *  Return the number of views (i.e. the size) in Matplotlib output
 *  due to generation process for MarkovOutTree class
 *
 **/

unsigned int MarkovOutTree::nb_generation_plot_set_computation() const
{
   unsigned int nb_plot_set = 0;
   MarkovOutTreeData::mvHistogram **histo = NULL;
   // DiscreteMultivariateReestimation<int> **histo = NULL;

   // compute the number of views

   if (generation_distributions != NULL)
   {
      if ((markov_data != NULL) && (markov_data->mot_reestimation != NULL))
         histo = markov_data->mot_reestimation->get_generation_reestim_ptr();
      nb_plot_set += generation_distributions->nb_plot_set_computation(histo);
   }

   return nb_plot_set;
}

/*****************************************************************
 *
 *  Matplotlib output for MarkovOutTree class
 *
 **/

 MultiPlotSet* MarkovOutTree::get_plotable() const
{ return get_plotable(markov_data); }

/*****************************************************************
 *
 *  Computation of the characteristic distributions
 *  for MarkovOutTree class using the considered variable,
 *  the tree depth (i.e. number of generations) and a flag
 *  on the computation of counting distributions
 *
 **/

void MarkovOutTree::characteristic_computation(int depth, bool counting_flag,
                                               int variable)
{
   bool computation[NB_OUTPUT_PROCESS+1];
   register int i, j, k;
   double sum, *memory = NULL;
   DiscreteParametric ddepth(UNIFORM, depth, depth, D_DEFAULT, D_DEFAULT);

   if (nb_component > 0)
   {

     // computation of the intensity and interval distributions
     // at state level

     if (((variable == I_DEFAULT) || (variable == 0)) &&
         ((!(nonparametric_process[0]->size != NULL))
             || (ddepth != *(nonparametric_process[0]->depth))))
     {
        computation[0] = true;
        nonparametric_process[0]->create_characteristic((Distribution)ddepth, false, counting_flag);

        if (memory == NULL)
        {
           switch (type)
           {
              case 'o' :
              {
                 memory = memory_computation();
                 break;
              }

              case 'e' :
              {
                 memory= new double[nb_state];
                 for(j = 0; j < nb_state; j++)
                    memory[j] = initial[j];
                 break;
              }
           }
        }

        for(i = 0; i < nb_state; i++)
           state_leave_probability(memory, i);

       // computation of the stationnary distribution in the case of
       // an equilibrium state process

       if (type == 'e')
       {
          sum = 0.;
          for(i = 0; i < nb_state; i++)
          {
             initial[i] = 1. / nb_state;
             sum += initial[i];
          }
          // the true formula seems to imply the computation of the
          // transposed transition matrix eigenvectors

          for(i = 0; i < nb_state; i++)
             initial[i] /= sum;
       }

       // index_state_distribution();

       for(i = 0; i < nb_state; i++)
       {
          state_no_occurrence_probability(i);
          state_first_occurrence_root_distribution(i);
          state_first_occurrence_leaves_distribution(i);

          state_sojourn_size_distribution(memory, i);
       }
   }
   else // nb_component < 1
      computation[0] = false;

     // computation of the intensity and interval distributions
     // at observation level

   for(i = 1; i <= nb_output_process; i++)
   {
      if ((nonparametric_process[i] != NULL) && ((variable == I_DEFAULT) || (i == variable)) &&
          ((nonparametric_process[i]->depth == NULL) || ((nonparametric_process[i]->size != NULL) &&
            (ddepth != *(nonparametric_process[i]->depth)))))
      {
         computation[i] = true;
         nonparametric_process[i]->create_characteristic(ddepth, true, counting_flag);

         // index_output_distribution(i);

         for(j = 0; j < nonparametric_process[i]->nb_value; j++)
         {
            output_no_occurrence_probability(i, j);
            output_first_occurrence_root_distribution(i, j);
            output_first_occurrence_leaves_distribution(i, j);

            output_leave_probability(memory, i, j);

            for(k = 0; k < nb_state; k++)
            {
               if ((nonparametric_process[i]->observation[k]->mass[j] > 0.) &&
                  ((state_type[k] != 'a') || (nonparametric_process[i]->observation[k]->mass[j] < 1.)))
                break;
            }

            if (k < nb_state)
               output_sojourn_size_distribution(memory, i, j);
            else
            {
               nonparametric_process[i]->absorption[j]= 1.;
               delete nonparametric_process[i]->sojourn_size[j];
               nonparametric_process[i]->sojourn_size[j]= NULL;
            }
         }
      }
      else
        computation[i] = false;
   }

   delete [] memory;

   if (counting_flag)
   {
     // computation of the counting distributions at state level

     if (computation[0])
        for (i = 0; i < nb_state; i++)
        {
           state_nb_pattern_mixture(i, 'r');
           state_nb_pattern_mixture(i, 'o');
        }

     // computation of the counting distributions at observation level

     for(i = 1; i <= nb_output_process; i++)
        if (computation[i])
           for (j = 0; j < nonparametric_process[i]->nb_value; j++)
           {
             output_nb_zones_mixture(i, j);
             output_nb_occurrences_mixture(i, j);
           }
      }
   }
}

/*****************************************************************
 *
 *  Computation of the characteristic distributions
 *  for MarkovOutTree class
 *  using a MarkovOutTreeData object,
 *  a flag on the computation of counting distributions,
 *  the considered variable and a flag indicating whether the tree
 *  depth (i.e. number of generations)
 *
 **/

void MarkovOutTree::characteristic_computation(const MarkovOutTreeData& tree,
                                               bool counting_flag,
                                               int variable,
                                               bool depth_flag)

{
   if (nb_component > 0)
   {
      register int i, j, k;
      int nb_value;
      double *memory= NULL;
      Distribution ddepth;

      if (tree.hdepth != NULL)
         ddepth = *(tree.hdepth);

      // computation of the characteristic distributions at state level

      if (((variable == I_DEFAULT) || (variable == 0)) &&
          // the 4 conditions below are nested
          ((!depth_flag) || ((depth_flag) && ((nonparametric_process[0]->depth == NULL)
           || (ddepth != *(nonparametric_process[0]->depth))))))
      {
         nonparametric_process[0]->create_characteristic(ddepth, true, counting_flag);

         memory = memory_computation();
         //index_state_distribution();

         for (i = 0; i < nb_state; i++)
         {
            state_no_occurrence_probability(i);
            if (tree.state_characteristics != NULL)
               nb_value = tree.state_characteristics->_max_value;
            else
               nb_value = 1;

            state_first_occurrence_root_distribution(i, nb_value);
            state_first_occurrence_leaves_distribution(i, nb_value);

            // state_first_occurrence_root_distribution(i);
            // state_first_occurrence_leaves_distribution(i);

            state_leave_probability(memory, i);

            if (state_type[i] != 'a')
            {
               if (tree.state_characteristics != NULL)
                  nb_value = tree.state_characteristics->_max_value;
               else
                  nb_value = 1;

               state_sojourn_size_distribution(memory, i, nb_value);

               // state_sojourn_size_distribution(memory, i);
            }
            else
            {
               nonparametric_process[0]->absorption[i]= 1.;
               delete nonparametric_process[0]->sojourn_size[i];
               nonparametric_process[0]->sojourn_size[i]= NULL;
            }

            if (counting_flag)
            {
               state_nb_pattern_mixture(i, 'r');
               state_nb_pattern_mixture(i, 'o');
            }
         }
      }

      // computation of the characteristic distributions at observation level

      for(i = 1; i <= nb_output_process; i++)
      {
         if ((nonparametric_process[i] != NULL) && (((variable == I_DEFAULT) || (i == variable) &&
             ((!depth_flag) || ((depth_flag) && ((nonparametric_process[i]->depth == NULL) ||
             ((nonparametric_process[i]->depth != NULL) && (ddepth != *(nonparametric_process[i]->depth)))))))))
         {
            nonparametric_process[i]->create_characteristic(ddepth, true, counting_flag);

            if (memory == NULL)
               memory = memory_computation();
            // index_output_distribution(i);

            for (j = 0; j < nonparametric_process[i]->nb_value; j++)
            {
               output_no_occurrence_probability(i, j);
               if (tree.characteristics[i-1] != NULL)
                  nb_value = tree.characteristics[i-1]->_max_value;
               else
                  nb_value = 1;

               output_first_occurrence_root_distribution(i, j, nb_value);
               output_first_occurrence_leaves_distribution(i, j, nb_value);

               output_leave_probability(memory, i, j);

               for (k= 0; k < nb_state; k++)
                  if ((nonparametric_process[i]->observation[k]->mass[j] > 0.) &&
                     ((state_type[k] != 'a') || (nonparametric_process[i]->observation[k]->mass[j] < 1.)))
                     break;

               if (k < nb_state)
               {
                  if (tree.characteristics[i-1] != NULL)
                     nb_value = tree.characteristics[i-1]->_max_value;
                  else
                     nb_value = 1;

                  output_sojourn_size_distribution(memory, i, j, nb_value);
               }
               else
               {
                  nonparametric_process[i]->absorption[j]= 1.;
                  delete nonparametric_process[i]->sojourn_size[j];
                  nonparametric_process[i]->sojourn_size[j]= NULL;
               }

               if (counting_flag)
               {
                  output_nb_zones_mixture(i, j);
                  output_nb_occurrences_mixture(i, j);
               }
            }
         }
      }
      delete [] memory;
   }
}
/*****************************************************************
 *
 *  Return the number of ordered children involved
 *  in the branching processes as factors (due to dependencies
 *  regarding these states, excluding parent state)
 *  in MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_children_branching() const
{
   if (generation_distributions == NULL)
      return 0;
   else
      return generation_distributions->nb_children_branching;
}

/*****************************************************************
 *
 *  Return the number of other factors involved
 *  in the branching process in MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_factors() const
{
   if (generation_distributions == NULL)
      return 0;
   else
      return generation_distributions->nb_factors;
}

/*****************************************************************
 *
 *  Return the number of possible values for each factor
 *  in MarkovOutTree
 *
 **/

GenerationProcess::index_array MarkovOutTree::get_factor_values() const
{
   GenerationProcess::index_array res(0);

   if ((generation_distributions == NULL)
       && (generation_distributions->factor_values != NULL))
      res = *generation_distributions->factor_values;

   return res;
}

/*****************************************************************
 *
 *  Return the number of generation processes in MarkovOutTree
 *  (combining all possible values of every factor)
 *
 **/

unsigned int MarkovOutTree::get_nb_generation_values() const
{
   if (generation_distributions == NULL)
      return 0;
   else
      return generation_distributions->nb_generation_values;
}

/*****************************************************************
 *
 *  Return the number of factors for generation processes
 *  (including parent state) in MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_generation() const
{
   if (generation_distributions == NULL)
      return 0;
   else
      return generation_distributions->nb_generation;
}

/*****************************************************************
 *
 *  Return the number of variable order Markov chains
 *  in MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_vomc() const
{ return nb_vomc; }

/*****************************************************************
 *
 *  Return the number of ordered children in MarkovOutTree
 *
 **/

unsigned int MarkovOutTree::get_nb_ordered_children() const
{ return unordered_rank - 1; }

/*****************************************************************
 *
 *  Return the minimal offsets in generation processes
 *  for MarkovOutTrees
 *
 **/

std::vector<unsigned int> MarkovOutTree::get_generation_min_offsets() const
{
   const unsigned int nb_generation_values = get_nb_generation_values();
   unsigned int i, s; // configuration of factors
   std::vector<unsigned int> res(nb_generation_values, UINT_MAX),
                             &offset = generation_distributions->generation[0]->offset;

   for(s = 0; s < nb_generation_values; s++)
   {
      offset = generation_distributions->generation[s]->offset;
      for(i = 0; i < nb_state; i++)
         res[s] = min(res[s], offset[i]);
   }
   return res;
}

/*****************************************************************
 *
 *  Return the marginal distributions associated with a generation process
 *  for MarkovOutTree class (a new instance is allocated)
 *
 **/

std::vector<Distribution*>
MarkovOutTree::get_distribution(StatError &error,
                                const GenerationProcess::index_array &fact) const
{
   register int i;
   int icode;
   Distribution *marginal = NULL;
   std::vector<Distribution*> res(0);

   error.init();

   if (generation_distributions == NULL)
      error.update(STAT_TREES_error[TREESTATR_NON_EXISTING_GENERATION_PROCESS]);
   else
   {
      icode = generation_distributions->code(fact);
      res.resize(nb_state);
      for(i = 0; i < nb_state; i++)
      {
         marginal = generation_distributions->generation[icode]->get_marginal(i);
         if (marginal != NULL)
            res[i] = marginal;
         else
            res[i] = NULL;
      }
   }

   return res;
}

/*****************************************************************
 *
 *  Return joint distribution of the descendants' types
 *  for a given factor in MarkovOutTrees
 *  (return a pointer; object should not be deallocated)
 *
 **/

const DiscreteMultivariateDistribution*
MarkovOutTree::get_generation_ptr(StatError &error,
                                  const GenerationProcess::index_array &fact) const
{
   register int i;
   int icode;
   DiscreteMultivariateDistribution *dist = NULL;

   error.init();

   if (generation_distributions == NULL)
      error.update(STAT_TREES_error[TREESTATR_NON_EXISTING_GENERATION_PROCESS]);
   else
   {
      icode = generation_distributions->code(fact);
      dist = generation_distributions->generation[icode];
   }

   return dist;
}

/*****************************************************************
 *
 *  Copy operator for MarkovOutTree class
 *
 **/

void MarkovOutTree::copy(const MarkovOutTree& markov, bool data_flag)
{
   register int i;

   // Chain::copy(markov) must be used before (or Chain::Chain)
   // as well as remove()

   if ((data_flag) && (markov.markov_data != NULL))
      markov_data = new MarkovOutTreeData(*(markov.markov_data), false, true);
   else
      markov_data = NULL;

   // copy observation processes
   nb_output_process = markov.nb_output_process;

   nonparametric_process = new CategoricalTreeProcess*[nb_output_process + 1];
   discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process + 1];
   continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process + 1];
   if ((markov.nonparametric_process != NULL) && (markov.nonparametric_process[0] != NULL))
      nonparametric_process[0] = new CategoricalTreeProcess(*markov.nonparametric_process[0]);

   for(i = 0; i < nb_output_process; i++)
   {
      if ((markov.nonparametric_process != NULL) && (markov.nonparametric_process[i+1] != NULL))
         nonparametric_process[i+1] = new CategoricalTreeProcess(*markov.nonparametric_process[i+1]);
      else
         nonparametric_process[i+1] = NULL;
      if ((markov.discrete_parametric_process != NULL) && (markov.discrete_parametric_process[i+1] != NULL))
         discrete_parametric_process[i+1] = new DiscreteParametricProcess(*markov.discrete_parametric_process[i+1]);
      else
         discrete_parametric_process[i+1] = NULL;
      if ((markov.continuous_parametric_process != NULL) && (markov.continuous_parametric_process[i+1] != NULL))
         continuous_parametric_process[i+1]  = new ContinuousParametricProcess(*markov.continuous_parametric_process[i+1]);
      else
         continuous_parametric_process[i+1] = NULL;
   }

   // copy models for ordered children and branching processes
   state_copy(markov);
}

/*****************************************************************
 *
 *  Destructor for MarkovOutTree class
 *
 **/

void MarkovOutTree::remove()
{
   register int i;

   // do not call Chain::remove, which then must be used before

   if (markov_data != NULL)
   {
      delete markov_data;
      markov_data = NULL;
   }

   if (nonparametric_process[0] != NULL)
   {
      delete nonparametric_process[0];
      nonparametric_process[0] = NULL;
   }

   for(i = 0; i < nb_output_process; i++)
   {
      if (nonparametric_process[i+1] != NULL)
      {
         delete nonparametric_process[i+1];
         nonparametric_process[i+1] = NULL;
      }
      if (discrete_parametric_process[i+1] != NULL)
      {
         delete discrete_parametric_process[i+1];
         discrete_parametric_process[i+1] = NULL;
      }
      if (continuous_parametric_process[i+1] != NULL)
      {
         delete continuous_parametric_process[i+1];
         continuous_parametric_process[i+1] = NULL;
      }
   }
   if (nonparametric_process != NULL)
   {
      delete [] nonparametric_process;
      nonparametric_process = NULL;
   }
   if (discrete_parametric_process != NULL)
   {
      delete [] discrete_parametric_process;
      discrete_parametric_process = NULL;
   }
   if (continuous_parametric_process != NULL)
   {
      delete [] continuous_parametric_process;
      continuous_parametric_process = NULL;
   }

   if (vomc_models != NULL)
   {
      for(i = 0; i < nb_vomc; i++)
         if (vomc_models[i] != NULL)
         {
            delete vomc_models[i];
            vomc_models[i] = NULL;
         }
      delete [] vomc_models;
      vomc_models = NULL;
   }
   if (generation_distributions != NULL)
   {
      delete generation_distributions;
      generation_distributions = NULL;
   }
}

/*****************************************************************
 *
 *  Copy operator for state process in MarkovOutTree class
 *
 **/

void MarkovOutTree::state_copy(const MarkovOutTree& markov)
{
   register int i;

   // copy models for ordered children
   // which must be uninitialized
   unordered_rank = markov.unordered_rank;
   nb_vomc = markov.nb_vomc;

   if ((nb_vomc > 0) && (markov.vomc_models != NULL))
   {
      vomc_models = new VariableOrderMarkov*[nb_vomc];
      for(i = 0; i < nb_vomc; i++)
         if (markov.vomc_models[i] != NULL)
            vomc_models[i] = new VariableOrderMarkov(*markov.vomc_models[i]);
   }
   else
      vomc_models = NULL;

   // copy branching processes
   // which must be uninitialized
   if (markov.generation_distributions != NULL)
      generation_distributions = new GenerationProcess(*markov.generation_distributions);
   else
      generation_distributions = NULL;
}

/*****************************************************************
 *
 *  Initialize output processes in MarkovOutTree class
 *  using the number of values for each process
 *  (I_DEFAULT for continuous processes)
 *
 **/

void MarkovOutTree::output_process_init(const int* nb_value)
{
   register int var;

   for(var = 1; var <= nb_output_process; var++)
   {
     if (nb_value[var-1] == I_DEFAULT)
     {
        nonparametric_process[var] = NULL;
        discrete_parametric_process[var] = NULL;
        continuous_parametric_process[var] = new ContinuousParametricProcess(nb_state);
     }
     else
        if (nb_value[var-1] <= NB_OUTPUT)
        {
           nonparametric_process[var] = new CategoricalTreeProcess(nb_state , nb_value[var-1], true);
           discrete_parametric_process[var] = NULL;
           continuous_parametric_process[var] = NULL;
        }
        else
        {
           nonparametric_process[var] = NULL;
           discrete_parametric_process[var] = new DiscreteParametricProcess(nb_state , (int)(nb_value[var-1] * SAMPLE_NB_VALUE_COEFF));
           continuous_parametric_process[var] = NULL;
        }

   }
}

/*****************************************************************
 *
 *  Copy operator for output process in MarkovOutTree class
 *
 **/

void MarkovOutTree::observation_copy(int inb_output_process,
                                     CategoricalProcess * const * npprocess,
                                     DiscreteParametricProcess * const * dpprocess,
                                     ContinuousParametricProcess * const * cpprocess)
{
   register int i;

   // copy ouput models
   // which must be uninitialized
   nb_output_process = inb_output_process;

   nonparametric_process = new CategoricalTreeProcess*[nb_output_process + 1];
   discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process + 1];
   continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process + 1];

   nonparametric_process[0] = new CategoricalTreeProcess(nb_state, nb_state);
   discrete_parametric_process[0] = NULL;
   continuous_parametric_process[0] = NULL;

   for(i = 1; i <= nb_output_process; i++)
   {
      if ((npprocess != NULL) && (npprocess[i - 1] != NULL))
      {
         nonparametric_process[i] = new CategoricalTreeProcess(*npprocess[i - 1]);
         discrete_parametric_process[i] = NULL;
         continuous_parametric_process[i] = NULL;
      }
      else
      {
         if ((dpprocess != NULL) && (dpprocess[i - 1] != NULL))
         {
            nonparametric_process[i] = NULL;
            discrete_parametric_process[i] = new DiscreteParametricProcess(*dpprocess[i - 1]);
            continuous_parametric_process[i] = NULL;
         }
         else
            if ((cpprocess != NULL) && (cpprocess[i - 1] != NULL))
            {
               nonparametric_process[i] = NULL;
               discrete_parametric_process[i] = NULL;
               continuous_parametric_process[i] = new ContinuousParametricProcess(*cpprocess[i - 1]);
            }
            else
            {
               nonparametric_process[i] = NULL;
               discrete_parametric_process[i] = NULL;
               continuous_parametric_process[i] = NULL;
            }
      }
   }
}

/*****************************************************************
 *
 *  Print MarkovIndOutTree and associated data structure
 *  using an output stream, a MarkovOutTreeData object,
 *  a flag on the level of detail, a flag on the file use
 *  and a flag on the hidden nature of model
 *
 **/

ostream& MarkovOutTree::ascii_write(ostream& os, const MarkovOutTreeData* otrees,
                                    bool exhaustive, bool file_flag,
                                    bool hidden) const
{
   register int i, j;
   int buff, width, variable, cumul_size;
   FrequencyDistribution **observation_dist = NULL,
                         *marginal_distribution = NULL;
   Histogram *marginal_histogram = NULL, **observation_histo = NULL;
   GenerationProcess::mvHistogram **empirical_observation = NULL;
   TreeCharacteristics *characteristics = NULL;

   os << STAT_TREES_word[TREESTATW_MARKOV_OUT_TREE] << endl;

   // Print number of states and initial probabilities
   os << "\n" << nb_state << " " << STAT_word[STATW_STATES] << endl;

   // compute column width

   width = column_width(nb_state, initial);
   for(i = 0;i < nb_row;i++)
   {
      buff = column_width(nb_state , transition[i]);
      if (buff > width)
         width = buff;
   }
   width += ASCII_SPACE;

   os << "\n";

   // print initial probabilities
   os << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;

   os << initial[0];
   for(i = 1; i < nb_state; i++)
       os << setw(width) << initial[i];

   os << endl << endl;

   os << unordered_rank - 1 << " " ;
   if (unordered_rank  > 2)
      os << STAT_TREES_word[TREESTATW_ORDERED_CHILDREN];
   else
      os << STAT_TREES_word[TREESTATW_ORDERED_CHILD];

   os << endl << endl;

   os << nb_vomc << " " ;
   if (nb_vomc > 1)
      os << STAT_TREES_word[TREESTATW_VARIABLE_ORDER_MARKOV_CHAINS];
   else
      os << STAT_TREES_word[TREESTATW_VARIABLE_ORDER_MARKOV_CHAIN];

   os << endl << endl;

   // print variable order markov models
   if (nb_vomc == 1)
   {
      // initial probabilities for *vomc_models[0] and *this
      // must match
      for(i = 0; i < nb_state; i++)
         vomc_models[0]->initial[i] = initial[i];
      if ((otrees != NULL) && (otrees->vomc_data != NULL))
         vomc_models[0]->ascii_write(os, otrees->vomc_data[0], exhaustive,
                                     file_flag, false);
      else
         vomc_models[0]->ascii_write(os, NULL, exhaustive,
                                     file_flag, false);
      os << endl;
   }
   if (nb_vomc > 1)
      for(i = 0; i < nb_vomc; i++)
      {
         os << STAT_TREES_word[TREESTATW_PARENT_STATE] << " " << i << " "
            << STAT_TREES_word[TREESTATW_VARIABLE_ORDER_MARKOV_CHAIN];
         if ((otrees != NULL) && (otrees->vomc_data != NULL))
            vomc_models[i]->ascii_write(os, otrees->vomc_data[i], exhaustive,
                                        file_flag, false);
         else
            vomc_models[i]->ascii_write(os, NULL, exhaustive,
                                        file_flag, false);
         os << endl;
      }


   if ((otrees != NULL) && (otrees->_type[0] == STATE)
       && (nonparametric_process[0]->size != NULL))
   {
      variable = 0;
      characteristics = otrees->characteristics[variable];
   }

   nonparametric_process[0]->ascii_print(os, 0, NULL, characteristics, exhaustive, file_flag);

   if (generation_distributions != NULL)
   {
      if ((otrees != NULL) && (otrees->mot_reestimation != NULL))
         empirical_observation = otrees->mot_reestimation->get_generation_reestim_ptr();
      generation_distributions->ascii_print(os, empirical_observation, exhaustive, file_flag);
   }

   // print the (characteristic ?) distributions of each
   // observed process
   // if (hidden)
   if (true)
   {
      for(i = 1; i <= nb_output_process; i++)
      {
         if (discrete_parametric_process[i] != NULL)
         {
            if (discrete_parametric_process[i]->weight != NULL)
            {
               os << "\n";
               if (file_flag)
                  os << "# ";
               os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << endl;
               if (file_flag)
                  os << "# ";
               for (j = 0; j < nb_state; j++)
                  os << discrete_parametric_process[i]->weight->mass[j] << "  ";
               os << endl;
            }

            if (discrete_parametric_process[i]->restoration_weight != NULL)
            {
               os << "\n";
               if (file_flag)
                 os << "# ";
               os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << endl;
               if (file_flag)
                 os << "# ";
               for (j = 0;j < nb_state;j++)
                  os << discrete_parametric_process[i]->restoration_weight->mass[j] << "  ";
               os << endl;
            }
            break;
         }
         else if (continuous_parametric_process[i] != NULL)
         {
            if (continuous_parametric_process[i]->weight)
            {
               os << "\n";
               if (file_flag)
                 os << "# ";
               os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << endl;
               if (file_flag)
                 os << "# ";
               for (j = 0;j < nb_state;j++)
                  os << continuous_parametric_process[i]->weight->mass[j] << "  ";
               os << endl;
            }

            if (continuous_parametric_process[i]->restoration_weight != NULL)
            {
               os << "\n";
               if (file_flag)
                 os << "# ";
               os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << endl;
               if (file_flag)
                 os << "# ";
               for (j = 0;j < nb_state;j++)
                  os << continuous_parametric_process[i]->restoration_weight->mass[j] << "  ";
               os << endl;
            }
            break;
         }
      }
      os << "\n" << nb_output_process << " "
         << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
   }
   else
   {
      os << "\n" << nb_output_process << " ";
      if (nb_output_process > 1)
         os << STAT_word[STATW_OUTPUT_PROCESSES] << endl;
      else
         os << STAT_word[STATW_OUTPUT_PROCESS] << endl;
   }

   // print distributions associated with each observation process

   for(i = 1; i <= nb_output_process; i++)
   {
      os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

      if (true)
      {
         os << " " << i;

         if (nonparametric_process[i] != NULL)
            os << " : " << STAT_word[STATW_NONPARAMETRIC];
         else
         {
            if (discrete_parametric_process[i] != NULL)
               os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC];
            else
              os << " : " << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
         }
      }
      os << endl;

      if (otrees != NULL)
      {
         if (otrees->_type[0] == STATE)
         {
            if ((nb_output_process > 1) && (variable < nb_output_process - 1))
               variable = i;
            else
               variable = i - 1;
         }
         else
            variable = i - 1;

         if (otrees->characteristics[variable] != NULL)
         {
            characteristics = otrees->characteristics[variable];
            marginal_distribution = otrees->characteristics[variable]->marginal_distribution;
            marginal_histogram = otrees->characteristics[variable]->marginal_histogram;
         }
         else
         {
            characteristics = NULL;
            marginal_histogram = NULL;
         }
         if (otrees->observation_distribution != NULL)
            observation_dist = otrees->observation_distribution[variable];
         else
            observation_dist = NULL;
         if (otrees->observation_histogram != NULL)
            observation_histo = otrees->observation_histogram[variable];
         else
            observation_histo = NULL;
      }

      if (nonparametric_process[i] != NULL)
         nonparametric_process[i]->ascii_print(os, i, observation_dist, characteristics,
                                               exhaustive, file_flag);
      else
      {
         if (discrete_parametric_process[i])
            discrete_parametric_process[i]->ascii_print(os, observation_dist,
                                                        marginal_distribution, exhaustive,
                                                        file_flag);
         else
           continuous_parametric_process[i]->ascii_print(os, observation_histo, observation_dist,
                                                         marginal_histogram,
                                                         marginal_distribution,
                                                         exhaustive, file_flag);
      }
   }

   if (otrees != NULL)
   {
      int nb_parameter;
      double information, likelihood, hidden_likelihood;

      if (hidden)
         nb_parameter = nb_parameter_computation(MIN_PROBABILITY);
      else
         nb_parameter = nb_parameter_computation(.0);

      hidden_likelihood = otrees->hidden_likelihood;
      likelihood = otrees->likelihood;

      // print the quantities for which the characteristic distributions
      // are invariant - if any

      os << "\n";
      if (file_flag)
         os << "# ";
      os << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      otrees->hsize->ascii_characteristic_print(os, false, file_flag);

      if (exhaustive)
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << "   | " << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
         otrees->hsize->ascii_print(os, file_flag);
      }

      os << "\n";
      if (file_flag)
        os << "# ";
      os << STAT_TREES_label[TREESTATL_CUMULATIVE_SIZE] << ": " << otrees->cumul_size_computation() << endl;

      // print the tree information quantity in the iid case

      information= otrees->iid_information_computation();

      os << "\n";
      if (file_flag)
         os << "# ";
      cumul_size= otrees->cumul_size_computation();
      os << STAT_TREES_label[TREESTATL_TREES_IID_INFORMATION] << ": " << information << " ("
         << information / cumul_size << ")" << endl;

      // print the likelihood

      if ((hidden) && (hidden_likelihood != D_INF))
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << STAT_TREES_label[TREESTATL_STATE_TREES_LIKELIHOOD] << ": " << hidden_likelihood << "   ("
            << STAT_label[STATL_NORMALIZED] << ": " << hidden_likelihood / cumul_size << ")" << endl;
      }

      if (likelihood != D_INF)
      {
         os << "\n";
         if (file_flag)
           os << "# ";
         os << STAT_TREES_label[TREESTATL_OBSERVED_TREES_LIKELIHOOD] << ": " << likelihood << "   ("
            << STAT_label[STATL_NORMALIZED] << ": " << likelihood / cumul_size << ")" << endl;
      }

      // print AIC, AICc and BIC
      if (likelihood != D_INF)
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
            << 2 * (likelihood - nb_parameter) << endl;

         if (nb_parameter < cumul_size-1)
         {
            os << "\n";
            if (file_flag)
               os << "# ";
            os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
               << 2*(likelihood-(double)(nb_parameter*cumul_size) /
                 (double)(cumul_size-nb_parameter-1)) << endl;
         }

         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
            << 2 * likelihood - nb_parameter * log((double)cumul_size) << endl;

         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << "): "
            << 2 * likelihood-penalty_computation(MIN_PROBABILITY) << endl;
      }

      // print ICL
      if (otrees->hidden_likelihood != D_INF)
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
            << 2*hidden_likelihood-nb_parameter*log((double)cumul_size) << endl;
         // otrees->likelihood == completed likelihood

         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << "): "
            << 2*hidden_likelihood-penalty_computation(MIN_PROBABILITY) << endl;
      }


      /* if ((likelihood != D_INF) && (nb_component == 1))
      {
         if (nb_parameter < cumul_size-1)
         {
            os << "\n";
            if (file_flag)
              os << "# ";
            os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
               << 2 * (likelihood - (double)(nb_parameter * cumul_size) /
                 (double)(cumul_size-nb_parameter-1)) << endl;
         }

         os << "\n";
         if (file_flag)
           os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
            << 2 * likelihood - nb_parameter * log((double)cumul_size) << endl;
      } */
   }

   return os;
}

ostream& MarkovOutTree::spreadsheet_write(ostream& os, const MarkovOutTreeData * tree,
                                          const Test * test) const
{ return os; }

bool MarkovOutTree::plot_write(const char * prefix, const char * title,
                               const MarkovOutTreeData * otrees) const
{
   bool status = true;

   return status;
}

/*****************************************************************
 *
 *  Matplotlib output of the MarkovOutTree model using data
 *
 **/

MultiPlotSet* MarkovOutTree::get_plotable(const MarkovOutTreeData * otrees) const
{
   register int var;
   int index = 0, nb_plot_set = 0;
   TreeCharacteristics *characteristics = NULL;
   Statiskit::Marginal::Multivariate::CountTable **histo = NULL;
   // DiscreteMultivariateReestimation<int> **histo = NULL;
   FrequencyDistribution **observation_dist = NULL; // discrete observation distribution
   Histogram **observation_histo = NULL; // continuous observation distribution
   MultiPlotSet *plot_set = NULL;

   // compute the number of views

   if (generation_distributions != NULL)
   {
      if ((otrees != NULL) && (otrees->mot_reestimation != NULL))
         histo = otrees->mot_reestimation->get_generation_reestim_ptr();
      nb_plot_set += generation_distributions->nb_plot_set_computation(histo);
   }
   for(var = 1; var <= nb_output_process; var++)
   {
      if (otrees != NULL)
      {
         if (otrees->observation_distribution != NULL)
            observation_dist = otrees->observation_distribution[var-1];
         else
            observation_dist = NULL;

         if (otrees->observation_histogram != NULL)
            observation_histo = otrees->observation_histogram[var-1];
         else
            observation_histo = NULL;

         if (otrees->characteristics[var-1] != NULL)
            characteristics = otrees->characteristics[var-1];
         else
            characteristics = NULL;
      }

      if (nonparametric_process[var] != NULL)
         nb_plot_set += nonparametric_process[var]->nb_plot_set_computation(var, observation_dist,
                                                                            characteristics, (otrees != NULL) ? otrees->hdepth : NULL);
      else
      {
         if ((otrees != NULL)
             && ((otrees->observation_distribution != NULL) || (otrees->observation_histogram != NULL)))
             nb_plot_set += nb_state;
         else
             nb_plot_set++;
         if ((discrete_parametric_process[var] != NULL) && (otrees != NULL)
              && (otrees->characteristics != NULL) && (otrees->characteristics[var-1]->marginal_distribution != NULL))
         {
            if ((discrete_parametric_process[var]->weight != NULL) &&
                (discrete_parametric_process[var]->mixture != NULL))
               nb_plot_set += 2;
            if ((discrete_parametric_process[var]->restoration_weight != NULL)
                && (discrete_parametric_process[var]->restoration_mixture != NULL))
               nb_plot_set += 2;
         }
         else
            if ((continuous_parametric_process[var] != NULL) && (otrees != NULL)
                 && (otrees->characteristics != NULL) && ((otrees->characteristics[var-1]->marginal_histogram != NULL)
                                                          || (otrees->characteristics[var-1]->marginal_distribution != NULL)))
            {
               if (continuous_parametric_process[var]->weight != NULL)
                  nb_plot_set += 2;
               if (continuous_parametric_process[var]->restoration_weight != NULL)
                  nb_plot_set += 2;
            }
      }
   }

   index = 0;
   plot_set = new MultiPlotSet(nb_plot_set, nb_output_process + 1);
   plot_set->variable_nb_viewpoint[0] = 0;
   plot_set->border = "15 lw 0";

   if (generation_distributions != NULL)
      generation_distributions->plotable_write(*plot_set, index, histo);
    // print output processes
   for(var = 1; var <= nb_output_process; var++)
   {
      if (otrees != NULL)
      {
         if (otrees->observation_distribution != NULL)
            observation_dist = otrees->observation_distribution[var-1];
         else
            observation_dist = NULL;

         if (otrees->observation_histogram != NULL)
            observation_histo = otrees->observation_histogram[var-1];
         else
            observation_histo = NULL;

         if (otrees->characteristics[var-1] != NULL)
            characteristics = otrees->characteristics[var-1];
         else
            characteristics = NULL;
      }

      if (nonparametric_process[var] != NULL)
      {
         plot_set->variable_nb_viewpoint[var] = 0;
         nonparametric_process[var]->plotable_write(*plot_set, index, var, observation_dist ,
                                                    characteristics, (otrees != NULL) ? otrees->hdepth : NULL);
      }
      else
      {
         if (discrete_parametric_process[var] != NULL)
            discrete_parametric_process[var]->plotable_write(*plot_set, index, var, observation_dist,
                                                             ((otrees != NULL) ? otrees->characteristics[var-1]->marginal_distribution : NULL));
         else
            continuous_parametric_process[var]->plotable_write(*plot_set, index, var,
                                                               observation_histo, observation_dist,
                                                               ((otrees != NULL) ? otrees->characteristics[var-1]->marginal_histogram : NULL),
                                                               ((otrees != NULL) ? otrees->characteristics[var-1]->marginal_distribution : NULL));
      }
   }

   return plot_set;
}

int MarkovOutTree::nb_parameter_computation(double min_probability) const
{
   return I_DEFAULT;
}

/*****************************************************************
 *
 *  Compute the distribution of each memory for a MarkovOutTree
 *  given the tree size distribution
 *
 **/

double* MarkovOutTree::memory_computation() const
{
   register int i, j, k;
   int power[ORDER], state_index[ORDER];
   double *memory, *state_tree, *pstate_tree, *states, *pstates, **ptransition;


   memory = new double[nb_row];
   for(i = 0; i < nb_row; i++)
      memory[i] = 0.;

   i = 1;
   // this line to be corrected and generalized
   for(j = 0; j < unordered_rank-1; j++)
   {
      power[j] = i;
      i *= nb_state;
   }

   // initialization of the tree state probabilities

   state_tree = new double[nb_row];
   pstate_tree = new double[nb_row];

   pstates = pstate_tree;
   i = 0;

   for(j = 0; j < nb_row;j++)
   {
      // if (j == self_row[i])
      //   *pstates++ = initial[i++];
      // else
        *pstates++ = 0.;
   }

   // computation of the probability of each memory, depending on the index
   for(i = 1; i < nonparametric_process[0]->depth->nb_value-2; i++)
   {
     // computation of the state tree probabilities having a depth
     // == unordered_children_rank-1

     for(j = 0; j < unordered_rank-1; j++)
        state_index[j]= 0;

     states = state_tree;

     for(j = 0;j < nb_row;j++)
     {
        ptransition = transition;
        pstates = pstate_tree;
        for (k = 0; k <unordered_rank-2; k++)
        {
           ptransition += state_index[k] * power[k+1];
           pstates += state_index[k] * power[k+1];
        }

        *states = 0.;
        for(k = 0; k < nb_state; k++)
        {
           *states += *(*ptransition + state_index[unordered_rank-2]) * *pstates++;
           ptransition++;
        }
        states++;

        // update of the state indices

        for (k= 0; k < unordered_rank-1; k++)
        {
           if (state_index[k] < nb_state - 1)
           {
              state_index[k]++;
              break;
           }
           else
             state_index[k]= 0;
        }
     }

     // update of the tree state probabilities and
     // of the cumulative memory probabilities

     states = state_tree;
     pstates = pstate_tree;

     for(j= 0; j < nb_row; j++)
     {
        memory[j] += *states * (1. - nonparametric_process[0]->depth->cumul[i]);
        *pstates++ = *states++;
     }
   }

   delete [] state_tree;
   delete [] pstate_tree;

   return memory;
}


void MarkovOutTree::state_no_occurrence_probability(int state, double increment)
{}

void MarkovOutTree::state_first_occurrence_root_distribution(int state,
                                                             int min_nb_value,
                                                             double cumul_threshold)
{}

void MarkovOutTree::state_first_occurrence_leaves_distribution(int state, int min_nb_value,
                                                               double cumul_threshold)
{}

void MarkovOutTree::state_leave_probability(const double * memory, int state,
                                            double increment)
{}

void MarkovOutTree::state_sojourn_size_distribution(const double * memory, int state,
                                                    int min_nb_value,
                                                    double cumul_threshold)
{}

void MarkovOutTree::state_nb_pattern_mixture(int state, char pattern)
{}

void MarkovOutTree::output_no_occurrence_probability(int variable, int output,
                                                     double increment)
{}

void MarkovOutTree::output_first_occurrence_root_distribution(int variable, int output,
                                                              int min_nb_value,
                                                              double cumul_threshold)
{}

void MarkovOutTree::output_first_occurrence_leaves_distribution(int variable, int output,
                                                                int min_nb_value,
                                                                double cumul_threshold)
{}

void MarkovOutTree::output_leave_probability(const double * memory,
                                             int variable, int output,
                                             double increment)
{}

void MarkovOutTree::output_sojourn_size_distribution(const double * memory, int variable,
                                                     int output, int min_nb_value,
                                                     double cumul_threshold)
{}

void MarkovOutTree::output_nb_zones_mixture(int variable, int output)
{}

void MarkovOutTree::output_nb_occurrences_mixture(int variable, int output)
{}

void MarkovOutTree::state_marginal_distribution(const MarkovOutTreeData& trees,
                                                double_array_3d& state_marginal,
                                                int index) const
{}

double** MarkovOutTree::state_marginal_distribution(const Trees& trees,
                                                    int index) const
{
   double **marginal = NULL;

   return marginal;

}

double MarkovOutTree::likelihood_correction(const MarkovOutTreeData& trees) const
{
   return D_INF;
}

double MarkovOutTree::penalty_computation(double min_probability) const
{  return D_INF; }

/*****************************************************************
 *
 *  Parse a MarkovOutTree from an input stream using a StatError object,
 *  the stream, and the current line number
 *
 **/

MarkovOutTree* Stat_trees::markov_out_tree_parsing(StatError& error,
                                                   std::ifstream &in_file,
                                                   int &line)
{
   RWLocaleSnapshot locale("en");
   RWCString buffer, token;
   size_t position;
   bool status = true, lstatus; //, **logic_transition;
   register int i, j;
   char type = 'v';
   int nb_state = 0, unordered_rank_read = -1, nb_vomc_read = -1,
       line_save;
   unsigned int unordered_rank = 0, nb_vomc = 0;
   long value;
   double proba, cumul;
   Chain *chain = NULL;
   VariableOrderMarkov **vomc = NULL;
   VariableOrderMarkovChain *vomch = NULL;
   GenerationProcess *pgeneration = NULL;
   MarkovOutTree *markov = NULL;

   while (buffer.readLine(in_file, false))
   {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS)
         buffer.remove(position);

      // read number of states and initial probabilities
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull()))
      {
         switch (i)
         {
            // test number of states
            case 0 :
            {
               lstatus = locale.stringToNum(token, &value);
               if (lstatus)
               {
                  if ((value < 2) || (value > NB_STATE))
                     lstatus = false;
                  else
                     nb_state = value;
               }
               if (!lstatus)
               {
                  status = false;
                  error.update(STAT_parsing[STATP_NB_STATE], line, i + 1);
               }
               break;
            }

            // test for STATES keyword
            case 1 :
            {
               if (token != STAT_word[STATW_STATES])
               {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_STATES] , line , i + 1);
               }
               break;
            }
         }
         i++;
      }

      if (i > 0)
      {
         if (i != 2)
         {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
         }
         break;
      }
   }

   if (nb_state == 0)
   {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
   }

   if (status)
   {
      chain = new Chain('o', nb_state, nb_state , false);

      // read initial probabilities
      while (buffer.readLine(in_file , false))
      {
         line++;

#        ifdef DEBUG
         cout << line << "  " << buffer << endl;
#        endif

         position = buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);
         i = 0;

         RWCTokenizer next(buffer);

         while (!((token = next()).isNull()))
         {
            // test for INITIAL_PROBABILITIES keyword
            if (i == 0)
            {
               if (token != STAT_word[STATW_INITIAL_PROBABILITIES])
               {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_word[STATW_INITIAL_PROBABILITIES], line);
               }
               i++;
            }
            if (i > 0)
            {
               if (i != 1)
               {
                  status = false;
                  error.update(STAT_parsing[STATP_FORMAT] , line);
               }
               break;
            }
         }
         if (i == 1)
            break;
      }

      while (buffer.readLine(in_file , false))
      {
         line++;

#        ifdef DEBUG
         cout << line << "  " << buffer << endl;
#        endif

         position = buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);
         i = 0;

         RWCTokenizer next(buffer);

         cumul = 0.;

         while (!((token = next()).isNull())) // parse 1 single line
         {
            if (i < nb_state)
            {
               lstatus = locale.stringToNum(token, &proba);
               if (lstatus)
               {
                  if ((proba < 0.) || (proba > 1. - cumul + DOUBLE_ERROR))
                     lstatus = false;
                  else
                  {
                     cumul += proba;
                     chain->initial[i] = proba;
                  }
               }

               if (!lstatus)
               {
                  status = false;
                  error.update(STAT_parsing[STATP_INITIAL_PROBABILITY], line, i + 1);
               }
            }
            i++;
         }

         if (i > 0)
         {
            if (i != nb_state)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT], line);
            }

            if (cumul < 1. - DOUBLE_ERROR)
            {
               status = false;
               error.update(STAT_parsing[STATP_PROBABILITY_SUM] , line);
            }
            break;
         }
      } // end while (buffer.readLine(in_file , false))
   }


   if (status)
   {
      // read unordered_rank;
      while (buffer.readLine(in_file, false))
      {
         line++;

   #     ifdef DEBUG
         cout << line << "  " << buffer << endl;
   #     endif

         position = buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);

         i = 0;

         RWCTokenizer next(buffer);

         while (!((token = next()).isNull()))
         {
            switch (i)
            {
               // test unordered_rank
               case 0 :
               {
                  lstatus = locale.stringToNum(token, &value);
                  if (lstatus)
                  {
                     if (value < 0)
                        lstatus = false;
                     else
                     {
                         unordered_rank_read = value;
                         unordered_rank = (unsigned int)(unordered_rank_read + 1);
                     }
                  }
                  if (!lstatus)
                  {
                     status = false;
                     error.update(STAT_TREES_parsing[TREESTATP_ORDERED_CHILDREN], line, i + 1);
                  }
                  break;
               }

               // test for keyword STATW_ORDERED_CHILD || STATW_ORDERED_CHILDREN

               case 1 :
               {
                  if (((unordered_rank_read > 1) && (token != STAT_TREES_word[TREESTATW_ORDERED_CHILDREN]))
                      || ((unordered_rank_read <= 1) && (token != STAT_TREES_word[TREESTATW_ORDERED_CHILD])))
                  {
                     status = false;
                     if (unordered_rank_read > 1)
                        error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_ORDERED_CHILDREN], line, i + 1);
                     else
                        error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_ORDERED_CHILD], line, i + 1);
                  }
                  break;
               }
            }
            i++;
         }

         if (i > 0)
         {
            if (i != 2)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT] , line);
            }
            break;
         }
      }

      if (unordered_rank_read == -1)
      {
         status = false;
         error.update(STAT_parsing[STATP_FORMAT], line);
      }
   } // end read unordered_rank;

   if (unordered_rank == 1)
      nb_vomc = 0;

   if ((status) && (unordered_rank > 1))
   {
      // read nb_vomc;
      while (buffer.readLine(in_file, false))
      {
         line++;

   #     ifdef DEBUG
         cout << line << "  " << buffer << endl;
   #     endif

         position = buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);

         i = 0;

         RWCTokenizer next(buffer);

         while (!((token = next()).isNull()))
         {
            switch (i)
            {
               // test nb_vomc
               case 0 :
               {
                  lstatus = locale.stringToNum(token, &value);
                  if (lstatus)
                  {
                     if ((value < 0) || (value > unordered_rank))
                        lstatus = false;
                     else
                     {
                         nb_vomc_read = value;
                         nb_vomc = (unsigned int)(nb_vomc_read);
                         line_save = line;
                     }
                  }
                  if (!lstatus)
                  {
                     status = false;
                     error.update(STAT_TREES_parsing[TREESTATP_NB_VOMC], line, i + 1);
                  }
                  break;
               }

               // test for keyword STATW_VARIABLE_ORDER_MARKOV_CHAIN || STATW_VARIABLE_ORDER_MARKOV_CHAINS

               case 1 :
               {
                  if (((nb_vomc > 1) && (token != STAT_TREES_word[TREESTATW_VARIABLE_ORDER_MARKOV_CHAINS]))
                      || ((nb_vomc <= 1) && (token != STAT_TREES_word[TREESTATW_VARIABLE_ORDER_MARKOV_CHAIN])))
                  {
                     status = false;
                     if (nb_vomc > 1)
                        error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_VARIABLE_ORDER_MARKOV_CHAINS], line, i + 1);
                     else
                        error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_TREES_word[TREESTATW_VARIABLE_ORDER_MARKOV_CHAIN], line, i + 1);
                  }
                  break;
               }
            }
            i++;
         }

         if (i > 0)
         {
            if (i != 2)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT] , line);
            }
            break;
         }
      }

   } // end read nb_vomc;

   if (nb_vomc == 0)
   {
      if (unordered_rank != 1)
      {
         status = false;
         error.update(STAT_parsing[TREESTATP_NB_VOMC], line_save);
      }
   }
   else
   {
      if ((nb_vomc == 1) || (nb_vomc == nb_state))
      {
         if (unordered_rank < 2)
         {
            status = false;
            error.update(STAT_parsing[TREESTATP_NB_VOMC], line_save);
         }
      }
      else
      {
         status = false;
         error.update(STAT_parsing[TREESTATP_NB_VOMC], line_save);
      }
   }

   if (nb_vomc > 0)
   {
      vomc = new VariableOrderMarkov*[nb_vomc];
      for(j = 0; j < nb_vomc; j++)
      {
         // read MARKOV_CHAIN keyword
         while (buffer.readLine(in_file , false))
         {
            line++;
#           ifdef DEBUG
            cout << line << "  " << buffer << endl;
#           endif

            position = buffer.first('#');
            if (position != RW_NPOS)
               buffer.remove(position);
            i = 0;

            RWCTokenizer next(buffer);

            while (!((token = next()).isNull()))
            {

               // test for (EQUILIBRIUM) MARKOV_CHAIN keyword

               if (i == 0)
               {
                  if (token == SEQ_word[SEQW_MARKOV_CHAIN])
                     type = 'o';
                  else
                  {
                     if (token == SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN])
                        type = 'e';
                     else
                     {
                        status = false;
                        ostringstream correction_message;
                        correction_message << SEQ_word[SEQW_MARKOV_CHAIN] << " or "
                                           << SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN];
                        error.correction_update(STAT_parsing[STATP_KEY_WORD], (correction_message.str()).c_str() , line);
                     }
                  }
               }
               i++;
            }

            if (i > 0)
            {
               if (i != 1)
               {
                  status = false;
                  error.update(STAT_parsing[STATP_FORMAT] , line);
               }
               break;
            }
         } // end while

         line_save = line;
         if (type != 'v')
         {
            vomch = variable_order_markov_parsing(error, in_file, line, type);
            if (vomch != NULL)
            {
               vomc[j] = new VariableOrderMarkov(vomch, 0, 0);
               delete vomch;
               vomch = NULL;
            }
         }
         else
            vomc[j] = NULL;

         type = 'v';
         if (vomc[j] != NULL)
         {
            vomc[j]->build_memory_transition();
            vomc[j]->build_previous_memory();
            if (vomc[j]->nb_state != nb_state)
               error.update(STAT_parsing[STATP_NB_STATE], line_save);
            if (vomc[j]->get_nb_output_process() != 0)
               error.update(STAT_parsing[STATP_NB_STATE], line_save, i + 1);
         }
      } // end for
   }

   // parse generation processes
   line_save = line;
   pgeneration = generation_process_parsing(error, in_file, line, nb_state);


   if (pgeneration != NULL)
   {
      if (pgeneration->get_nb_children_branching() > (unordered_rank - 1))
      {
         status = false;
         error.update(STAT_parsing[TREESTATP_NB_CHILDREN_BRANCHING], line_save);
      }
      markov = new MarkovOutTree(*chain, 0, NULL, NULL, NULL,
                                 unordered_rank, nb_vomc,
                                 vomc, pgeneration);

   }
   if (pgeneration != NULL)
   {
      delete pgeneration;
      pgeneration = NULL;
   }
   if (chain != NULL)
   {
      delete chain;
      chain = NULL;
   }

   if (vomc != NULL)
   {
      for(i = 0; i < nb_vomc; i++)
         if (vomc[i] != NULL)
         {
            delete vomc[i];
            vomc[i] = NULL;
         }
      delete [] vomc;
      vomc = NULL;
   }
   return markov;
}

/*****************************************************************
 *
 *  Create a MarkovOutTree from a file
 *  using a StatError object, the path,
 *  and the tree size (used to compute the characteristic distributions)
 *
 **/

MarkovOutTree* Stat_trees::markov_out_tree_ascii_read(StatError& error,
                                                      const char * path,
                                                      int size,
                                                      double cumul_threshold)
{
   RWLocaleSnapshot locale("en");
   RWCString buffer, token;
   size_t position;
   bool status = false, lstatus = true;
   register int i;
   long value;
   int line, nb_output_process = I_DEFAULT, output_process_type, index, nb_state;
   VariableOrderMarkov **vomc = NULL;
   CategoricalProcess **nonparametric_observation = NULL;
   DiscreteParametricProcess **discrete_parametric_observation = NULL;
   ContinuousParametricProcess **continuous_parametric_observation = NULL;
   MarkovOutTree *markov = NULL,
                 *markov_cp = NULL; // result without output processes
   ifstream in_file(path);

   markov = NULL;
   error.init();

   if (!in_file)
      error.update(STAT_error[STATR_FILE_NAME]);
   else
   {
      status = true;
      line = 0;

      if (size < 2)
      {
         status= false;
         error.update(STAT_TREES_error[TREESTATR_SMALL_TREE_SIZE]);
      }

      if (size > MAX_SIZE)
      {
         status= false;
         error.update(STAT_TREES_error[TREESTATR_BIG_TREE_SIZE]);
      }

      while (buffer.readLine(in_file , false))
      {
         line++;

 #       ifdef DEBUG
         cout << line << "  " << buffer << endl;
 #       endif

         position = buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);
         i = 0;

         RWCTokenizer next(buffer);

         while (!((token = next()).isNull()))
         {
            // test for MARKOV_OUT_TREE keyword
            if (i == 0)
            {
               if (token != STAT_TREES_word[TREESTATW_MARKOV_OUT_TREE])
               {
                  status = false;
                  ostringstream correction_message;
                  correction_message << STAT_TREES_word[TREESTATW_MARKOV_OUT_TREE];
                  error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                          (correction_message.str()).c_str() , line);
               }
            }
            i++;
         }

         if (i > 0)
         {
            if (i != 1)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT] , line);
            }
            break;
         }
      }

      if (status)
      {
         // read MarkovOutTree (without output processes)

         markov_cp = markov_out_tree_parsing(error, in_file, line);

         if (markov_cp != NULL)
         {
            nb_state = markov_cp->nb_state;

            // read observation processes
            nb_output_process = I_DEFAULT;

            nonparametric_observation = NULL;
            discrete_parametric_observation = NULL;
            continuous_parametric_observation = NULL;

            while (buffer.readLine(in_file , false))
            {
               line++;

#              ifdef DEBUG
               cout << line << "  " << buffer << endl;
#              endif

               position = buffer.first('#');
               if (position != RW_NPOS)
                  buffer.remove(position);

               i = 0;

               RWCTokenizer next(buffer);

               while (!((token = next()).isNull()))
               {
                  switch (i)
                  {
                     // test number of output processes
                     case 0 :
                     {
                        lstatus = locale.stringToNum(token, &value);
                        if (lstatus)
                        {
                           if ((value < 0) || (value > NB_OUTPUT_PROCESS))
                              lstatus = false;
                           else
                               nb_output_process = value;
                        }

                        if (!lstatus)
                        {
                           status = false;
                           error.update(STAT_parsing[STATP_NB_OUTPUT_PROCESS], line, i + 1);
                        }
                        break;
                     }

                     // test for OUTPUT_PROCESS(ES) keyword

                     case 1 :
                     {
                        if (token != STAT_word[nb_output_process <= 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES])
                        {
                           status = false;
                           error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                                   STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES], line, i + 1);
                        }
                        break;
                     }
                  }
                  i++;
               }

               if (i > 0)
               {
                  if (i != 2)
                  {
                     status = false;
                     error.update(STAT_parsing[STATP_FORMAT], line);
                  }
                  break;
               }
            } // end while readlines

            if (nb_output_process == I_DEFAULT)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT], line);
            }

            else
            {
               nonparametric_observation = new CategoricalProcess*[nb_output_process];
               discrete_parametric_observation = new DiscreteParametricProcess*[nb_output_process];
               continuous_parametric_observation = new ContinuousParametricProcess*[nb_output_process];

               for (i = 0; i < nb_output_process; i++)
               {
                  nonparametric_observation[i] = NULL;
                  discrete_parametric_observation[i] = NULL;
                  continuous_parametric_observation[i] = NULL;
               }

               index = 0;

               while (buffer.readLine(in_file , false))
               {
                  line++;

#                 ifdef DEBUG
                  cout << line << "  " << buffer << endl;
#                 endif

                  position = buffer.first('#');
                  if (position != RW_NPOS)
                     buffer.remove(position);
                  i = 0;

                  RWCTokenizer next(buffer);

                  while (!((token = next()).isNull()))
                  {
                     switch (i)
                     {

                        // test mot cle OUTPUT_PROCESS
                        case 0 :
                        {
                           output_process_type = CATEGORICAL;

                           if (token == STAT_word[STATW_OUTPUT_PROCESS])
                             index++;
                           else
                           {
                              status = false;
                              error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OUTPUT_PROCESS] , line , i + 1);
                           }
                           break;
                        }

                        // test index of output process

                        case 1 :
                        {
                           lstatus = locale.stringToNum(token , &value);
                           if ((lstatus) && ((value != index) || (value > nb_output_process)))
                             lstatus = false;

                           if (!lstatus)
                           {
                              status = false;
                              error.update(STAT_parsing[STATP_OUTPUT_PROCESS_INDEX] , line , i + 1);
                           }
                           break;
                        }

                        // test separator

                        case 2 :
                        {
                           if (token != ":")
                           {
                              status = false;
                              error.update(STAT_parsing[STATP_SEPARATOR], line, i + 1);
                           }
                           break;
                        }

                        // test for NONPARAMETRIC / DISCRETE_PARAMETRIC / CONTINUOUS_PARAMETRIC keywords

                        case 3 :
                        {
                           if (token == STAT_word[STATW_NONPARAMETRIC])
                             output_process_type = CATEGORICAL;
                           else
                           {
                              if ((token == STAT_word[STATW_DISCRETE_PARAMETRIC]) ||
                                    (token == STAT_word[STATW_PARAMETRIC]))
                                 output_process_type = DISCRETE_PARAMETRIC;
                              else
                              {
                                 if (token == STAT_word[STATW_CONTINUOUS_PARAMETRIC])
                                    output_process_type = CONTINUOUS_PARAMETRIC;
                                 else
                                 {
                                    status = false;
                                    ostringstream correction_message;
                                    correction_message << STAT_word[STATW_NONPARAMETRIC] << " or "
                                                       << STAT_word[STATW_DISCRETE_PARAMETRIC] << " or "
                                                       << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
                                    error.correction_update(STAT_parsing[STATP_KEY_WORD],
                                                            (correction_message.str()).c_str(),
                                                            line, i + 1);
                                 }
                              }
                           }
                           break;
                        }
                     } // end switch
                     i++;
                  } // end while parse line

                  if (i > 0)
                  {
                     if (i != 4)
                     {
                        status = false;
                        error.update(STAT_parsing[STATP_FORMAT] , line);
                     }

                     switch (output_process_type)
                     {
                        case CATEGORICAL :
                        {
                           nonparametric_observation[index - 1] = categorical_observation_parsing(error, in_file, line,
                                                                                                  nb_state, HIDDEN_MARKOV,
                                                                                                  true);
                           if (nonparametric_observation[index - 1] == NULL)
                             status = false;
                           break;
                        }

                        case DISCRETE_PARAMETRIC :
                        {
                           discrete_parametric_observation[index - 1] = discrete_observation_parsing(error, in_file, line,
                                                                                                     nb_state, HIDDEN_MARKOV,
                                                                                                     cumul_threshold);
                           if (discrete_parametric_observation[index - 1] == NULL)
                              status = false;
                           break;
                        }

                        case CONTINUOUS_PARAMETRIC :
                        {
                           continuous_parametric_observation[index - 1] = continuous_observation_parsing(error, in_file, line,
                                                                                                         nb_state, VON_MISES);
                           // VON_MISES may be erroneous

                           if (!continuous_parametric_observation[index - 1])
                              status = false;
                           break;
                        }
                     } // end switch
                  } // end if (i > 0)
               } // end while readline

               // check whether there is unrelevant information
               // in the sequel
               while (buffer.readLine(in_file , false))
               {
                  line++;

            #     ifdef DEBUG
                  cout << line << " " << buffer << endl;
            #     endif

                  position = buffer.first('#');
                  if (position != RW_NPOS)
                     buffer.remove(position);
                  if (!(buffer.isNull()))
                  {
                     status = false;
                     error.update(STAT_parsing[STATP_FORMAT] , line);
                  }
               }

               if (index != nb_output_process)
               {
                  status = false;
                  error.update(STAT_parsing[STATP_FORMAT] , line);
               }

               if (status)
               {
                  markov = new MarkovOutTree(markov_cp, nb_output_process,
                                             nonparametric_observation,
                                             discrete_parametric_observation,
                                             continuous_parametric_observation,
                                             size);

#                 ifdef DEBUG
                  markov->ascii_write(cout);
#                 endif

               }

               delete markov_cp;
               markov_cp = NULL;

               for (i = 0;i < nb_output_process;i++)
               {
                  delete nonparametric_observation[i];
                  delete discrete_parametric_observation[i];
                  delete continuous_parametric_observation[i];
               }
               delete [] nonparametric_observation;
               delete [] discrete_parametric_observation;
               delete [] continuous_parametric_observation;
            }
         }  // end (markov_cp != NULL)
      } // end if status
   } // end if (!in_file)
   return markov;
}

/*****************************************************************
 *
 * Left (bit) shift operator of a MarkovOutTree
 *
 **/

std::ostream& Stat_trees::operator<<(std::ostream& os, const MarkovOutTree& markov)
{ return markov.ascii_write(os, markov.markov_data); }

/*****************************************************************
 *
 *  Constructor of MarkovOutTreeData class
 *  using the number of integral and float variables,
 *  and the number of trees (composed of a single vertex)
 *
 **/

MarkovOutTreeData::MarkovOutTreeData(unsigned int inb_integral,
                                     unsigned int inb_float,
                                     unsigned int inb_trees)
 : Trees(inb_integral, inb_float, inb_trees)
 , markov(NULL)
 , chain_data(NULL)
 , vomc_data(NULL)
// , empirical_observation(NULL)
 , hdepth(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(0)
 , state_trees(NULL)
 , virt_vertices(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
 , mot_reestimation(NULL)
{
   Typed_edge_one_int_tree::value default_state_value;
   value default_value;
   key v;
   unsigned int t;

   if (inb_trees > 0)
   {

      default_state_value.Int() = I_DEFAULT;
      default_value.reset(inb_integral, inb_float);

      state_trees = new Typed_edge_one_int_tree*[inb_trees];
      virt_vertices = new virtual_vdic*[inb_trees];

      for(t = 0; t < inb_trees; t++)
      {
         trees[t] = new tree_type(0, 0, default_value);
         v = trees[t]->add_vertex(default_value);
         state_trees[t] = new Typed_edge_one_int_tree(0, 0, default_state_value);
         assert(v == state_trees[t]->add_vertex(default_state_value));
         virt_vertices[t] = new virtual_vdic();
      }
   }
}

/*****************************************************************
 *
 *  Copy constructor of MarkovOutTreeData class using a flag
 *  on copying markov_data and a flag on copying the characteristics
 *
 **/

MarkovOutTreeData::MarkovOutTreeData(const MarkovOutTreeData& trees,
                                     bool model_flag,
                                     bool characteristic_flag)
 : Trees(trees)
 , markov(NULL)
 , chain_data(NULL)
 , vomc_data(NULL)
// , empirical_observation(NULL)
 , hdepth(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(0)
 , state_trees(NULL)
 , virt_vertices(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
 , mot_reestimation(NULL)
{ copy(trees, model_flag, characteristic_flag); }

/*****************************************************************
 *
 *  Constructor of MarkovOutTreeData class
 *  using the number of trees, the type of each variable
 *  and the multidimensional observed trees
 *
 **/

MarkovOutTreeData::MarkovOutTreeData(int inb_trees,
                                     const int * itype,
                                     Default_tree **otrees)
 : Trees(inb_trees, itype, otrees)
 , markov(NULL)
 , chain_data(NULL)
 , vomc_data(NULL)
 , hdepth(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(0)
 , state_trees(NULL)
 , virt_vertices(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
 , mot_reestimation(NULL)
{
   unsigned int t;

   virt_vertices = new virtual_vdic*[inb_trees];

   for(t = 0; t < inb_trees; t++)
      virt_vertices[t] = new virtual_vdic();
}

/*****************************************************************
 *
 *  Constructor of MarkovOutTreeData class
 *  using a Trees object and the potential state variable
 *
 **/

MarkovOutTreeData::MarkovOutTreeData(const Trees& otrees, int istate_variable)
 : Trees()
 , markov(NULL)
 , chain_data(NULL)
 , vomc_data(NULL)
 , hdepth(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(0)
 , state_trees(NULL)
 , virt_vertices(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
 , mot_reestimation(NULL)
{
   unsigned int t;

   _nb_integral = otrees.get_nb_int();
   _nb_float = otrees.get_nb_float();
   _nb_trees = otrees.get_nb_trees();

   virt_vertices = new virtual_vdic*[_nb_trees];

   for(t = 0; t < _nb_trees; t++)
      virt_vertices[t] = new virtual_vdic();

   if (_nb_integral + _nb_float <= 0)
   {
      trees = new tree_type*[_nb_trees];
      for(t = 0; t < _nb_trees; t++)
         trees[t] = new tree_type();

   }
   else
   {
      if (istate_variable == I_DEFAULT)
         Trees::copy(otrees, true);
      else
      {
         // delete istate_variable from output variables
         // and add it as state variable
         bool status = true;
         unsigned int s = 0, var;
         StatError error;
         Trees *selected_variables = NULL;
         int_array ivariable = NULL;
         Unlabelled_typed_edge_tree *tmp_utree = NULL;
         Int_fl_container i;
         One_int_container si;
         Trees::vertex_iterator it, end;

         if (_nb_integral + _nb_float > 1)
         {
            // copy non-state variables
            ivariable = new int[_nb_integral + _nb_float - 1];
            for(var = 0; var < _nb_integral + _nb_float; var++)
               if (var != istate_variable)
                  ivariable[s++] = var + 1;

            selected_variables = otrees.select_variable(error, _nb_integral + _nb_float - 1,
                                                        ivariable, true);
            delete [] ivariable;
            ivariable = NULL;
            if (error.get_nb_error() > 0)
               status = false;
            else
            {
               Trees::remove();
               Trees::copy(*selected_variables, true);
            }
            delete selected_variables;
            selected_variables = NULL;
         }
         else
            _nb_integral--; // no output variable
         if (status)
         {
            // copy state_variable
            ivariable = new int[1];
            ivariable[0] = istate_variable+1;

            selected_variables = otrees.select_variable(error, 1, ivariable, true);
            delete [] ivariable;
            ivariable = NULL;
            if (error.get_nb_error() > 0)
            {
               status = false;
               delete selected_variables;
               selected_variables = NULL;
            }
         }
         if (status)
         {
            // _nb_integral--; should have been handled by "copy"
            state_trees = new Typed_edge_one_int_tree*[_nb_trees];

            for(t = 0; t < _nb_trees; t++)
            {
               state_trees[t] = new Typed_edge_one_int_tree(0, 0);

               tmp_utree = otrees.get_tree_ptr(t)->get_structure();
               state_trees[t]->set_structure(*tmp_utree);
               Tree_tie::tie(it, end) = state_trees[t]->vertices();
               while (it < end)
               {
                  i = otrees.get_tree_ptr(t)->get(*it);
                  // add state variable
                  si.Int() = i.Int(istate_variable);
                  _nb_states = max(_nb_states, si.Int());
                  state_trees[t]->put(*it++, si);
               }
               delete tmp_utree;
               tmp_utree = NULL;
            }
            _nb_states++;

            build_characteristics();
            build_size_frequency_distribution();
            build_nb_children_frequency_distribution();
            build_observation_frequency_distribution();
         }
      } // end (istate_variable != I_DEFAULT)
   }
}

/*****************************************************************
 *
 *  Destructor for MarkovOutTreeData class
 *
 **/

MarkovOutTreeData::~MarkovOutTreeData()
{
   remove();
   Trees::remove();
}

/*****************************************************************
 *
 *  Return the size of a given tree of MarkovOutTreeData
 *
 **/

unsigned int MarkovOutTreeData::get_size(int index) const
{
   assert((index >= 0) && (index < _nb_trees));
   if (state_trees != NULL)
      return state_trees[index]->get_size();
   else
      return Trees::get_size(index);
}


/*****************************************************************
 *
 *  Allocate the FrequencyDistributions associated with
 *  the output distributions in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::create_observation_frequency_distribution(int nb_state)
{
   register int var, j;

   if (_nb_integral > 0)
   {
      if (observation_distribution == NULL)
      {
         observation_distribution = new FrequencyDistribution**[_nb_integral];
         for(var= 0; var < _nb_integral; var++)
         {
            observation_distribution[var] = new FrequencyDistribution*[nb_state];
            for(j = 0; j < nb_state; j++)
               observation_distribution[var][j] = new FrequencyDistribution(get_max_int_value(var)+1);
         }
      }
      else // the number of values may have changed
         for(var = 0; var < _nb_integral; var++)
         {
            if (observation_distribution[var] != NULL)
               for(j = 0; j < nb_state; j++)
               {
                  if (observation_distribution[var][j] != NULL)
                  {
                     delete observation_distribution[var][j];
                     observation_distribution[var][j] = NULL;
                  }
                  observation_distribution[var][j] = new FrequencyDistribution(get_max_int_value(var)+1);
               }
            else
            {
               observation_distribution[var] = new FrequencyDistribution*[nb_state];
               for(j= 0; j < nb_state; j++)
                  observation_distribution[var][j] = new FrequencyDistribution(get_max_int_value(var)+1);
            }
         }
   }
   else // delete the potential frequency distributions
   {
      if (observation_distribution != NULL)
         for(var = 0; var < _nb_integral; var++)
         {
            if (observation_distribution[var] != NULL)
            {
               for(j = 0; j < nb_state; j++)
               {
                  if (observation_distribution[var][j] != NULL)
                  {
                     delete observation_distribution[var][j];
                     observation_distribution[var][j] = NULL;
                  }
               }
               delete [] observation_distribution[var];
               observation_distribution[var] = NULL;
            }
         }
      observation_distribution = NULL;
   }
}

/*****************************************************************
 *
 *  Compute the FrequencyDistributions associated with
 *  the output distributions in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::observation_frequency_distribution_computation()
{
   register int var;

   for(var = 0; var < _nb_integral; var++)
      observation_frequency_distribution_computation(var);
}

/*****************************************************************
 *
 *  Allocate and compute the FrequencyDistributions associated with
 *  the output distributions in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::build_observation_frequency_distribution()
{
   build_state_characteristics();
   create_observation_frequency_distribution(state_characteristics->marginal_distribution->nb_value);
   observation_frequency_distribution_computation();
}

/*****************************************************************
 *
 *  Allocate and compute the FrequencyDistributions associated with
 *  state process in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::build_state_frequency_distribution(StatError &error,
                                                           unsigned int inb_children_branching,
                                                           unsigned int inb_factors,
                                                           const GenerationProcess::index_array *ifactor_values,
                                                           unsigned int inb_generation,
                                                           unsigned int iunordered_rank,
                                                           unsigned int inb_vomc,
                                                           const std::vector<unsigned int>* factor_indices)
{
   unsigned int i;
   bool status = true;
   ostringstream error_message;

   update_frequency_distributions();

   error.init();

   if (((inb_vomc > 0) && (iunordered_rank < 2))
       || ((inb_vomc < 1) && (iunordered_rank > 1)))
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_NB_VOMC] << " (: "
                    << inb_vomc << ") or " << STAT_TREES_parsing[TREESTATP_ORDERED_CHILDREN]  << " (: "
                    << iunordered_rank << ")";
      error.update((error_message.str()).c_str());
   }
   if ((inb_vomc > 0) && (inb_vomc != 1) && (inb_vomc != _nb_states))
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_NB_VOMC] << "s : "
                    << inb_vomc;
      error.update((error_message.str()).c_str());
   }
   if ((inb_generation - inb_children_branching - inb_factors < 0)
       || (inb_generation - inb_children_branching - inb_factors > 1))
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_NB_GENERATION_PROCESS] << "(: "
                    << inb_generation << ") or " << STAT_TREES_parsing[TREESTATP_NB_CHILDREN_BRANCHING] << "(: "
                    << inb_children_branching << ") or " << STAT_TREES_parsing[TREESTATP_NB_FACTORS] << ")";
      error.update((error_message.str()).c_str());
   }
   if (inb_factors > 0)
   {
      if ((ifactor_values == NULL) || (factor_indices == NULL))
      {
         status = false;
         error_message << STAT_TREES_parsing[TREESTATP_NB_FACTORS]
                       << " or missing values of factors";
         error.update((error_message.str()).c_str());
      }
      else
         for(i=0; i < inb_factors; i++)
         {
            if ((*factor_indices)[i] >= _nb_integral)
            {
               status = false;
               if (error.get_nb_error() > 0)
                  error_message << "; ";
               error_message << STAT_error[STATR_VARIABLE_INDEX]
                             << ": " << (*factor_indices)[i];
               error.update((error_message.str()).c_str());
            }
            if ((*ifactor_values)[i] < _max_value.Int(i))
            {
               status = false;
               if (error.get_nb_error() > 0)
                  error_message << "; ";
               error_message << STAT_TREES_parsing[TREESTATP_FACTOR_VALUE]
                             << ": " << (*ifactor_values)[i];
               error.update((error_message.str()).c_str());
            }
         }
   }
   if (inb_children_branching > iunordered_rank - 1)
   {
      status = false;
      error_message << STAT_TREES_parsing[TREESTATP_ORDERED_CHILDREN]  << " (: "
                    << iunordered_rank - 1 << ") or " << STAT_TREES_parsing[TREESTATP_NB_CHILDREN_BRANCHING] << "(: "
                    << inb_children_branching << ")";
      error.update((error_message.str()).c_str());
   }

   if (status)
   {
      build_chain_reestimation('o', _nb_states);
      build_markov_reestimation(inb_children_branching, inb_factors, ifactor_values,
                                inb_generation, iunordered_rank, inb_vomc, factor_indices);
   }
}

/*****************************************************************
 *
 *  Return the number of states of MarkovOutTreeData
 *
 **/

int MarkovOutTreeData::get_nb_states() const
{ return _nb_states; }


/*****************************************************************
 *
 *  Return model associated with the data in MarkovOutTreeData
 *  (a new instance is allocated)
 *
 **/

MarkovOutTree* MarkovOutTreeData::get_markov() const
{
   MarkovOutTree *res = NULL;

   if (markov != NULL)
      res = new MarkovOutTree(*markov);

   return markov;
}


/*****************************************************************
 *
 *  Return model associated with the data in MarkovOutTreeData
 *  (return a pointer; object should not be deallocated)
 *
 **/

MarkovOutTree* MarkovOutTreeData::get_markov_ptr() const
{ return markov; }

/*****************************************************************
 *
 *  Return the set of state trees in MarkovOutTreeData
 *  (a new instance is allocated)
 *
 **/

MarkovOutTreeData::ptOne_int_tree_set MarkovOutTreeData::get_state_trees() const
{
   int t;
   ptOne_int_tree_set res_trees = NULL;

   if (this->state_trees != NULL)
   {
      res_trees = new Typed_edge_one_int_tree*[_nb_trees];
      for(t = 0; t < _nb_trees; t++)
      {
         if (this->state_trees[t] != NULL)
            res_trees[t] = new Typed_edge_one_int_tree(*(this->state_trees[t]));
         else
            res_trees[t] = NULL;
      }
   }
   return res_trees;
}

/*****************************************************************
 *
 *  Return the set of state trees in MarkovOutTreeData
 *  (return a pointer; object should not be deallocated)
 *
 **/

MarkovOutTreeData::ptOne_int_tree_set MarkovOutTreeData::get_state_trees_ptr() const
{ return state_trees; }

/*****************************************************************
 *
 *  Return the marginal distributions associated with a generation process
 *  for MarkovOutTreeData class (a new instance is allocated)
 *
 **/

std::vector<DiscreteDistributionData*>
MarkovOutTreeData::get_distribution_data(StatError &error,
                                         const MarkovOutTreeData::index_array &fact) const
{
   bool status = true, markov_dist = false;
   register int i;
   StatError markov_error;
   // DiscreteMultivariateReestimation<int> *reestim = NULL;
   mvHistogram *reestim = NULL;
   std::vector<DiscreteDistributionData*> res(0);
   FrequencyDistribution *cfd = NULL;
   uHistogram_ptr_type chisto;
   std::vector<Distribution*> dist(0);
   Distribution* cdist = NULL;

   error.init();

   if (mot_reestimation == NULL)
   {
      status = false;
      error.update(STAT_TREES_error[TREESTATR_CHARACTERISTICS_NOT_COMPUTED]);
   }
   else
   {
      reestim = mot_reestimation->get_generation_reestim_ptr(error, fact);
      if (error.get_nb_error() > 0)
         status = false;
#     ifdef DEBUG
      print_stat_histogram(*reestim, cout);
#     endif
   }
   if (status)
   {
      res.resize(_nb_states);
      if (markov != NULL)
      {
         dist = markov->get_distribution(markov_error, fact);

         if (markov_error.get_nb_error() == 0)
            markov_dist = true;
      }

      for(i = 0; i < _nb_states; i++)
      {
         chisto = Statiskit::marginalizing(*reestim, i);
         cfd = Scalar_stat_histogram_to_frequency_distribution(*chisto);
         // chisto = reestim->get_marginal_frequency(i);
         if (markov_dist)
            cdist = dist[i];
         else
            cdist = NULL;
         res[i] = new DiscreteDistributionData(*cfd, cdist);
         delete cfd;
         cfd = NULL;
         if (markov_dist)
         {
            delete dist[i];
            dist[i] = NULL;
         }
      }
   }

   return res;
}

/*****************************************************************
 *
 *  Return observed vectors associated with a generation process
 *  together with their weights for MarkovOutTree class
 *
 **/

std::vector<std::pair<std::vector<unsigned int>, unsigned int> >
MarkovOutTreeData::get_joint_distribution(StatError &error,
                                          const GenerationProcess::index_array &fact) const
{
   bool status = true;
   register int i;
   int icode;
   Statiskit::Marginal::Univariate::CountTable::const_iterator it;
   // DiscreteMultivariateReestimation<int> *reestim = NULL;
   mvHistogram *reestim = NULL;
   std::vector<std::vector<unsigned int> > elements;
   std::pair<std::vector<unsigned int>, unsigned int> p; // pair vector, weight
   std::vector<std::pair<std::vector<unsigned int>, unsigned int> > res(0);

   error.init();

   if (mot_reestimation == NULL)
   {
      status = false;
      error.update(STAT_TREES_error[TREESTATR_CHARACTERISTICS_NOT_COMPUTED]);
   }
   else
   {
      reestim = mot_reestimation->get_generation_reestim_ptr(error, fact);
      if (error.get_nb_error() > 0)
         status = false;
   }
   if (status)
   {
      elements = Stat_histogram_get_elements(*reestim);
      res.resize(elements.size());
      for(i = 0; i < elements.size(); i++)
      {
         p.first = elements[i];
         p.second = Stat_trees::Stat_histogram_get_frequency(*reestim,
                                                             elements[i]);
         // p.second = (unsigned int)(*it).second;
         res[i] = p;
      }
   }
   return res;
}


/*****************************************************************
 *
 *  Print a summary of the data in MarkovOutTreeData class
 *
 **/

ostream& MarkovOutTreeData::ascii_write(ostream& os, bool exhaustive) const
{ return ascii_write(os, exhaustive, false); }

/*****************************************************************
 *
 *  Matplotlib output for MarkovOutTreeData class
 *
 **/

MultiPlotSet* MarkovOutTreeData::get_plotable() const
{
   register int var;
   int index = 0, nb_plot_set = 0;
   DiscreteMultivariateReestimation<int> **histo = NULL;
   FrequencyDistribution **observation_dist = NULL; // discrete observation distribution
   Histogram **observation_histo = NULL; // continuous observation distribution
   MultiPlotSet *plot_set = NULL;

   if (markov != NULL)
      return markov->get_plotable(this);
   else
   {
      // compute the number of views
      nb_plot_set += generation_nb_plot_set_computation();

      for(var = 0; var < _nb_integral; var++)
      {
         if (observation_distribution[var] != NULL)
            nb_plot_set += _nb_states;
         else
            nb_plot_set++;
      }

      index = 0;
      plot_set = new MultiPlotSet(nb_plot_set, _nb_integral + 1);
      plot_set->variable_nb_viewpoint[0] = 0;
      plot_set->border = "15 lw 0";

      generation_plotable_write(*plot_set, index);
       // print output processes
      for(var = 0; var < _nb_integral; var++)
         discrete_output_plotable_write(*plot_set, index, var);
      return plot_set;
      // also consider Trees::get_plotable();
   }
}

/*****************************************************************
 *
 *  Update the FrequencyDistributions of the joint distribution
 *  of children states in a MarkovOutTreeData using the number
 *  of states, the number of children states, of factors
 *  and parent children determining the number inb_generation
 *  of branching processes. Compute the sequences of ordered children.
 *  build_nb_children_frequency_distribution and min_max_value_computation
 *  must be called before (or min_offset and hnb_children must be up-to-date)
 *  If the arguments are incompatible with previously defined model,
 *  an error is returned
 *
 **/

void MarkovOutTreeData::update_markov_reestimation(StatError& error,
                                                   unsigned int inb_vomc,
                                                   unsigned int inb_ordered_children,
                                                   unsigned int inb_children_branching,
                                                   const std::vector<unsigned int>& factor_variables,
                                                   bool parent_dependent)
{
   bool status = true, first_error = true;
   const unsigned int iunordered_rank = inb_ordered_children+1,
                      inb_factors = factor_variables.size();
   unsigned int inb_generation = 0, inb_generation_values = 1, i;
   ostringstream error_message;
   GenerationProcess::index_array *ifactor_max_values = NULL;
   std::vector<int> generation_types_estimation(0), generation_types(0);
   std::vector<bool> skipped_variables(_nb_integral+_nb_float, false);

   error.init();

   if (markov != NULL)
   {
      if (inb_children_branching != markov->get_nb_children_branching())
      {
         status = false;
         if (first_error)
         {
            error_message << "Parameters set for generation process are not "
                          << "compatible with current model: " << endl;
            first_error = false;
         }
         error_message << "Number of ordered children on which "
                       << "depend the generation processes should be: "
                       << markov->get_nb_children_branching() << endl;
      }
      if (inb_factors != markov->get_nb_factors())
      {
         status = false;
         if (first_error)
         {
            error_message << "Parameters set for generation process are not "
                          << "compatible with current model: " << endl;
            first_error = false;
         }
         error_message << "Number discrete variables on which "
                       << "depend the generation processes should be: "
                       << markov->get_nb_factors() << endl;

      }
      if (inb_factors > 0)
      {
         ifactor_max_values = new GenerationProcess::index_array(inb_factors, 0);
         *ifactor_max_values = markov->get_factor_values();
         for(i=0; i < inb_factors; i++)
            if ((*ifactor_max_values)[i] != _max_value.Int(factor_variables[i])+1)
            {
               status = false;
               if (first_error)
               {
                  error_message << "Parameters set for generation process are not "
                                << "compatible with current model: " << endl;
                  first_error = false;
               }
               error_message << "Bad maximal value for factor " << _max_value.Int(factor_variables[i])+1
                             << "; should be " << (*ifactor_max_values)[i] << endl;

            }
      }

      if (factor_variables.size() > 0)
      {
         if (inb_factors != factor_variables.size())
         {
            status = false;
            if (first_error)
            {
               error_message << "Parameters set for generation process are not "
                             << "compatible with current model: " << endl;
               first_error = false;
            }
            error_message << STAT_TREES_parsing[TREESTATP_FACTOR_VALUE]
                          << ": " << factor_variables.size()
                          << "; should be " << inb_factors << endl;
         }
         for(i=0; i < inb_factors; i++)
            if (factor_variables[i] >= _nb_integral)
            {
               status = false;
               if (first_error)
               {
                  error_message << "Parameters set for generation process are not "
                                << "compatible with current model: " << endl;
                  first_error = false;
               }
               error_message << "Bad value for factor variable: "
                             << factor_variables[i]
                             << "; should be < " << _nb_integral << endl;
            }
         if (*ifactor_max_values != markov->get_factor_values())
         {
            status = false;
            if (first_error)
            {
               error_message << "Parameters set for generation process are not "
                             << "compatible with current model: " << endl;
               first_error = false;
            }
            error_message << "Values of discrete variables on which "
                          << "depend the generation processes should be: "
                          << markov->get_factor_values() << endl;
         }
      }

      if (inb_generation != markov->get_nb_generation())
      {
         status = false;
         if (first_error)
         {
            error_message << "Parameters set for generation process are not "
                          << "compatible with current model: " << endl;
            first_error = false;
         }
         error_message << "Total number of variables on which "
                       << "depend the generation processes should be: "
                       << markov->get_nb_generation() << endl;
      }

      if (iunordered_rank != markov->get_nb_ordered_children()+1)
      {
         status = false;
         if (first_error)
         {
            error_message << "Parameters set for generation process are not "
                          << "compatible with current model: " << endl;
            first_error = false;
         }
         error_message << "Number of ordered children "
                       << "in state process should be: "
                       << markov->get_nb_ordered_children() << endl;
      }

      if (inb_vomc != markov->get_nb_vomc())
      {
         status = false;
         if (first_error)
         {
            error_message << "Parameters set for generation process are not "
                          << "compatible with current model: " << endl;
            first_error = false;
         }
         error_message << "number of variable order Markov chains "
                       << "for ordered children should be: "
                       << markov->get_nb_vomc() << endl;
      }
   }
   else
      // markov == NULL
      if (inb_factors > 0)
      {
         ifactor_max_values = new GenerationProcess::index_array(inb_factors, 0);
         for(i=0; i < inb_factors; i++)
            (*ifactor_max_values)[i] = _max_value.Int(factor_variables[i])+1;
      }


   if (status)
   {
      if (mot_reestimation != NULL)
      {
         delete mot_reestimation;
         mot_reestimation = NULL;
      }
      status = markov_out_tree_check_estimation(error, inb_vomc, inb_ordered_children,
                                                inb_children_branching, factor_variables,
                                                parent_dependent, generation_types,
                                                inb_generation, inb_generation_values,
                                                generation_types_estimation,
                                                ifactor_max_values, skipped_variables);

      if (status)
         build_markov_reestimation(inb_children_branching, inb_factors, ifactor_max_values,
                                   inb_generation, iunordered_rank, inb_vomc, &factor_variables);
   }
   else
      error.update((error_message.str()).c_str());

   delete ifactor_max_values;
   ifactor_max_values = NULL;
}

/*****************************************************************
 *
 *  Return a given state tree in MarkovOutTreeData
 *  (a new instance is allocated)
 *
 **/

Typed_edge_one_int_tree* MarkovOutTreeData::get_state_tree(int itree) const
{
   Typed_edge_one_int_tree *res_tree;

   assert(itree < _nb_trees);

   if (this->state_trees != NULL)
   {
      if (this->state_trees[itree] != NULL)
         res_tree = new Typed_edge_one_int_tree(*(this->state_trees[itree]));
      else
         res_tree = NULL;
   }
   else
      res_tree = NULL;

   return res_tree;
}

/*****************************************************************
 *
 *  Return a MarkovOutTreeData containing the states
 *  as a variable. If no state tree is present, return a copy
 *  of this.
 *
 **/

MarkovOutTreeData* MarkovOutTreeData::get_state_markov_out_tree_data() const
{
   register int offset;
   const short int added_int_variables = 1;
   MarkovOutTreeData *res = NULL;
   Unlabelled_typed_edge_tree *tmp_utree = NULL;
   bool return_this = false;
   int inb_variables, inb_trees, var, t;
   int *itype = NULL;
   Int_fl_container i(_nb_integral+added_int_variables, 0), s;
   state_tree_type::vertex_iterator it, end;
   Default_tree **otrees = NULL;

   if ((this->state_trees == NULL) || (markov == NULL))
   {
      // state tree is not present
      return_this = true;
      res = new MarkovOutTreeData(*this, true, true);
   }
   else
   {
      inb_trees = _nb_trees;
      inb_variables = _nb_integral + 0 + added_int_variables;
      itype = new int[inb_variables];

      for(var = 0; var < added_int_variables; var++)
         itype[var] = STATE;

      offset = added_int_variables;
      for(var = 0; var < _nb_integral; var++)
         itype[var+offset] = _type[var];

      offset += _nb_integral;
      for(var = 0; var < _nb_float; var++)
         itype[var+offset] = _type[var+_nb_integral];

      otrees = new Default_tree*[inb_trees];
      for(t = 0; t < inb_trees; t++)
      {
         tmp_utree = state_trees[t]->get_structure();
         otrees[t] = new Default_tree(_nb_integral+added_int_variables,
                                      _nb_float+0,
                                      trees[t]->root(), 1);
         otrees[t]->set_structure(*tmp_utree, i);
         Tree_tie::tie(it, end) = state_trees[t]->vertices();
         while (it < end)
         {
            // output variables
            if (_nb_integral +_nb_float > 0)
               s = trees[t]->get(*it);
            // add state variable
            i.Int(0) = (this->state_trees[t]->get(*it)).Int();
            // copy existing integer variables
            for(var = 0; var < _nb_integral; var++)
               i.Int(var+added_int_variables) = s.Int(var);
            // copy existing floating variables
            for(var = 0; var < _nb_float; var++)
               i.Double(var+0) = s.Double(var);
            otrees[t]->put(*it++, i);
         }
         delete tmp_utree;
         tmp_utree = NULL;
      }

      res = new MarkovOutTreeData(inb_trees, itype, otrees);
      res->markov = new MarkovOutTree(*markov, false);

      res->state_trees = new Typed_edge_one_int_tree*[inb_trees];
      for(t = 0; t < inb_trees; t++)
         res->state_trees[t] = new Typed_edge_one_int_tree(*(this->state_trees[t]));
      res->likelihood = likelihood;
      res->hidden_likelihood = hidden_likelihood;
      res->_nb_states = _nb_states;

      // res->chain_data= new ChainData(*res, 0, 1, markov);
      res->chain_data = new ChainData(markov->type, markov->nb_state,
                                      markov->nb_state);
      res->build_characteristics();
      res->build_size_frequency_distribution();
      res->build_nb_children_frequency_distribution();
      res->build_observation_frequency_distribution();

      if (vomc_data != NULL)
      {
         res->vomc_data = new VariableOrderMarkovData*[1];
         res->vomc_data[0] = new VariableOrderMarkovData(*vomc_data[0]);
      }

      if (hdepth != NULL)
         res->hdepth = new FrequencyDistribution(*hdepth);

      if (virt_vertices != NULL)
      {
         res->virt_vertices = new virtual_vdic*[_nb_trees];
         for(t = 0; t < _nb_trees; t++)
         {
             if (virt_vertices[t] != NULL)
                res->virt_vertices[t] =
                   new virtual_vdic(*virt_vertices[t]);
             else
                res->virt_vertices[t] = NULL;
         }
      }

      res->markov->characteristic_computation(*res, true);

      if (mot_reestimation != NULL)
         res->mot_reestimation = new MarkovOutTreeReestimation<int>(*mot_reestimation);


      delete [] itype;
      itype = NULL;
      for(t = 0; t < inb_trees; t++)
      {
         delete otrees[t];
         otrees[t] = NULL;
      }
      delete [] otrees;
      otrees = NULL;
   }
   return res;

}


/*****************************************************************
 *
 *  Test for virtuality of a vertex in MarkovTreeOutData
 *  (true: either virtual or with censored number of children)
 *
 **/

bool MarkovOutTreeData::is_virtual(int itree, int ivertex) const
{
   virtual_vdic *vdic = NULL;

   assert((itree >= 0) && (itree < _nb_trees)  && (state_trees != NULL));

   if (state_trees[itree]->get_nb_children(ivertex) > 0)
      return false;
   else
      if (virt_vertices[itree]->count(ivertex) == 0)
      // by default, vertices not in dictionary are not virtual
         return false;
      else
      // key ivertex is in dictionary
      {
         vdic = virt_vertices[itree];
         return (*vdic)[ivertex];
      }
}

/*****************************************************************
 *
 *  Return dictionary of virtual vertices for a given tree
 *  in MarkovTreeOutData (a new instance is allocated)
 *
 **/

std::map<int, bool>* MarkovOutTreeData::get_virtual_vertices(int itree) const
{
   virtual_vdic *vdic = NULL;

   assert((itree >= 0) && (itree < _nb_trees));

   vdic = new virtual_vdic(*virt_vertices[itree]);

   return vdic;
}

/*****************************************************************
 *
 *  Return dictionary of virtual vertices for a given tree
 *  in MarkovTreeOutData (return a pointer; object should not be deallocated)
 *
 **/

const std::map<int, bool> * MarkovOutTreeData::get_virtual_vertices_ptr(int itree) const
{
   assert((itree >= 0) && (itree < _nb_trees));

   const virtual_vdic *vdic = virt_vertices[itree];

   return vdic;
}

/*****************************************************************
 *
 *  Set virtual status of a vertex in MarkovTreeOutData
 *
 **/

void MarkovOutTreeData::set_virtual_vertex(int itree, int ivertex, bool bvirtual) const
{
   assert((itree >= 0) && (itree < _nb_trees));
   assert((ivertex >= 0) &&
           (((state_trees[itree] != NULL) && (ivertex < state_trees[itree]->get_size()))
            || (trees[itree] != NULL) && (ivertex < trees[itree]->get_size())));

   virt_vertices[itree]->insert(pair<int, bool>(ivertex, bvirtual));

}

/*****************************************************************
 *
 *  Update histogram for the number of children
 *  and minimal / maximal values of variables,
 *  taking into account virtual vertices
 *
 **/

void MarkovOutTreeData::update_frequency_distributions()
{
   this->build_nb_children_frequency_distribution();
   this->min_max_value_computation();
}

/*****************************************************************
 *
 *  Get every possible combination of factors for generation processes
 * (including parent state) in MarkovTreeOutData.
 * Parent state comes first (?) if included
 **/

std::vector<GenerationProcess::index_array> MarkovOutTreeData::get_factor_combinations() const
{
   std::vector<GenerationProcess::index_array> res(0);

   if (mot_reestimation != NULL)
      res = mot_reestimation->get_factor_combinations();

   return res;
}

/*****************************************************************
 *
 *  Return the set of state trees if MarkovTreeData
 *  (return a pointer; object should not be deallocated)
 *
 **/

Typed_edge_one_int_tree* MarkovOutTreeData::get_state_tree_ptr(int itree) const
{
   assert((itree < _nb_trees) && (this->state_trees != NULL));
   return state_trees[itree];
}

/*****************************************************************
 *
 *  Select subset of trees of MarkovOutTreeData by identifiers
 *  using a StatError object, the number of trees, tree identifiers
 *  and a flag on keeping or rejecting the selected trees
 *
 **/

MarkovOutTreeData* MarkovOutTreeData::select_individual(StatError& error,
                                                        int inb_trees,
                                                        int_array iidentifier,
                                                        bool keep) const

{
   bool status = true;
   unsigned int t;
   int_array index = NULL;
   Trees *res_trees = NULL;
   MarkovOutTreeData *res = NULL;
   ptOne_int_tree_set res_state_trees = NULL;

   res_trees = Trees::select_individual(error, inb_trees, iidentifier, keep);
   if (error.get_nb_error() > 0)
   {
      status = false;
      return res;
   }
   else
   {
      res = new MarkovOutTreeData(*res_trees, true, true);
      delete res_trees;
      res_trees = NULL;
   }

   if (status)
      if (state_trees != NULL)
      {
         res_state_trees = state_trees_select_individual(error, inb_trees, iidentifier,
                                                         index, keep);
         if (error.get_nb_error() > 0)
         {
            if (res != NULL)
            {
               delete res;
               res = NULL;
            }
            if (res_state_trees != NULL)
            {
               delete res_state_trees;
               res_state_trees = NULL;
            }
            status = false;
            return res;
         }
         else
         {
            res->state_trees = res_state_trees;
            // note that if there is no ouput variable, _max_size, _max_depth
            // hsize, hnb_children, trees will be incorrect
            if (_nb_integral + _nb_float == 0)
            {
               res->min_max_value_computation();
               res->max_size_computation();
               res->max_depth_computation();
               res->max_nb_children_computation();
               res->cumul_size_computation();
               res->build_size_frequency_distribution();
               res->cumul_nb_children_computation();
               res->build_nb_children_frequency_distribution();
            }

            if (virt_vertices != NULL)
            {
               res->virt_vertices = new virtual_vdic*[inb_trees];
               for(t = 0; t < inb_trees; t++)
                  if (virt_vertices[index[t]] != NULL)
                     res->virt_vertices[t] = new virtual_vdic(*virt_vertices[index[t]]);
            }
            res->_nb_states = _nb_states;
            delete [] index;
            index = NULL;
         }
      }

   return res;
}


/*****************************************************************
 *
 *  Copy operator for MarkovOutTreeData class
 *
 **/

void MarkovOutTreeData::copy(const MarkovOutTreeData& otrees, bool model_flag,
                             bool characteristic_flag)
{
   register int i, j, var, t;
   // Trees::copy(otrees) must be used before (or Trees::Trees)
   // as well as remove()

   hidden_likelihood = otrees.hidden_likelihood;
   likelihood = otrees.likelihood;
   _nb_states = otrees._nb_states;

   if ((model_flag) && (otrees.markov != NULL))
      markov = new MarkovOutTree(*otrees.markov, false);
   else
      markov = NULL;

   if (otrees.chain_data != NULL)
      chain_data = new ChainData(*otrees.chain_data);
   else
      chain_data = NULL;

   if (otrees.vomc_data != NULL)
   {
      vomc_data = new VariableOrderMarkovData*[1];
      vomc_data[0] = new VariableOrderMarkovData(*(otrees.vomc_data[0]));
   }
   else
      vomc_data = NULL;

/*   if (otrees.empirical_observation != NULL)
   {
      empirical_observation = new DiscreteMultivariateReestimation<int>*[_nb_states];
      for(i = 0; i < _nb_states; i++)
      {
         if (otrees.empirical_observation[i] != NULL)
            empirical_observation[i] = new DiscreteMultivariateReestimation<int>(*otrees.empirical_observation[i]);
      }
   }
   else
      empirical_observation = NULL; */

   if (otrees.hdepth != NULL)
      hdepth = new FrequencyDistribution(*otrees.hdepth);
   else
      hdepth = NULL;


   if (otrees.state_trees != NULL)
   {
      state_trees = new Typed_edge_one_int_tree*[_nb_trees];
      for(t = 0; t < _nb_trees; t++)
      {
          if (otrees.state_trees[t] != NULL)
             state_trees[t] =
                new Typed_edge_one_int_tree(*otrees.state_trees[t]);
          else
             state_trees[t] = NULL;
      }
   }
   else
      state_trees = NULL;

   if (otrees.virt_vertices != NULL)
   {
      virt_vertices = new virtual_vdic*[_nb_trees];
      for(t = 0; t < _nb_trees; t++)
      {
          if (otrees.virt_vertices[t] != NULL)
             virt_vertices[t] =
                new virtual_vdic(*otrees.virt_vertices[t]);
          else
             virt_vertices[t] = NULL;
      }
   }
   else
      virt_vertices = NULL;

   if (characteristic_flag)
   {
      if (otrees.state_characteristics != NULL)
         state_characteristics = new TreeCharacteristics(*otrees.state_characteristics);
      else
         state_characteristics = NULL;

      if (otrees.observation_distribution != NULL)
      {
         observation_distribution = new FrequencyDistribution**[_nb_integral];
         for(var = 0; var < _nb_integral; var++)
            if (otrees.observation_distribution[var] != NULL)
            {
               observation_distribution[var] = new FrequencyDistribution*[state_characteristics->marginal_distribution->nb_value];
               for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
               {
                  if (otrees.observation_distribution[var][j] != NULL)
                     observation_distribution[var][j] = new FrequencyDistribution(*otrees.observation_distribution[var][j]);
                  else
                     observation_distribution[var][j] = NULL;
               }
            }
            else
               observation_distribution[var] = NULL;
      }
      else
         observation_distribution = NULL;

      if (otrees.observation_histogram != NULL)
      {
         observation_histogram = new Histogram**[_nb_integral];
         for(var = 0; var < _nb_integral; var++)
            if (otrees.observation_histogram[var] != NULL)
            {
               observation_histogram[var] = new Histogram*[state_characteristics->marginal_distribution->nb_value];
               for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
               {
                  if (otrees.observation_histogram[var][j] != NULL)
                     observation_histogram[var][j] = new Histogram(*otrees.observation_histogram[var][j]);
                  else
                     observation_histogram[var][j] = NULL;
               }
            }
            else
               observation_histogram[var] = NULL;
      }
      // ptHistogram_array_2d observation_histogram;

      if (otrees.mot_reestimation != NULL)
         mot_reestimation = new MarkovOutTreeReestimation<int>(*otrees.mot_reestimation);
      else
         mot_reestimation = NULL;

   } // end if (characteristic_flag)
   else
   {
      if (state_characteristics != NULL)
         delete state_characteristics;
      if (observation_distribution != NULL)
         delete observation_distribution;
      if (observation_histogram != NULL)
         delete observation_histogram;
      if (mot_reestimation != NULL)
         delete mot_reestimation;
      state_characteristics = NULL;
      observation_distribution = NULL;
      observation_histogram = NULL;
      mot_reestimation = NULL;
   }
}

/*****************************************************************
 *
 *  Destructor for MarkovOutTreeData class
 *
 **/

void MarkovOutTreeData::remove()
{
   register int i, j, var, t;
   // Trees::remove() must be used before


   if (markov != NULL)
   {
      delete markov;
      markov = NULL;
   }

   if (chain_data != NULL)
   {
      delete chain_data;
      chain_data = NULL;
   }

   if (vomc_data != NULL)
   {
      if (vomc_data[0] != NULL)
      {
         delete vomc_data[1]; // to be changed
         vomc_data[0] = NULL;
      }
      delete [] vomc_data;
      vomc_data = NULL;
   }

/*   if (empirical_observation != NULL)
   {
      for(i = 0; i < _nb_states; i++)
         if (empirical_observation[i] != NULL)
         {
            delete empirical_observation[i];
            empirical_observation[i] = NULL;
         }
      delete [] empirical_observation;
      empirical_observation = NULL;
   }*/

   if (hdepth != NULL)
   {
      delete hdepth;
      hdepth = NULL;
   }

   if (state_trees != NULL)
   {
      for(t = 0; t < _nb_trees; t++)
          if (state_trees[t] != NULL)
          {
             delete state_trees[t];
             state_trees[t] = NULL;
          }
      delete [] state_trees;
      state_trees = NULL;
   }

   if (virt_vertices != NULL)
   {
      for(t = 0; t < _nb_trees; t++)
          if (virt_vertices[t] != NULL)
          {
             delete virt_vertices[t];
             virt_vertices[t] = NULL;
          }
      delete [] virt_vertices;
      virt_vertices = NULL;
   }

   if (observation_distribution != NULL)
   {
      assert(state_characteristics != NULL);
      for(var = 0; var < _nb_integral; var++)
         if (observation_distribution[var] != NULL)
         {
            for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
               if (observation_distribution[var][j] != NULL)
               {
                  delete observation_distribution[var][j];
                  observation_distribution[var][j] = NULL;
               }
            delete [] observation_distribution[var];
            observation_distribution[var] = NULL;
         }
      delete [] observation_distribution;
      observation_distribution = NULL;
   }

   if (observation_histogram != NULL)
   {
      assert(state_characteristics != NULL);
      for(var = 0; var < _nb_integral; var++)
         if (observation_histogram[var] != NULL)
         {
            for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
               if (observation_histogram[var][j] != NULL)
               {
                  delete observation_histogram[var][j];
                  observation_histogram[var][j] = NULL;
               }
            delete [] observation_histogram[var];
            observation_histogram[var] = NULL;
         }
      delete [] observation_histogram;
      observation_histogram = NULL;
   }

   if (state_characteristics != NULL)
   {
      delete state_characteristics;
      state_characteristics = NULL;
   }

   if (mot_reestimation != NULL)
   {
      delete mot_reestimation;
      mot_reestimation = NULL;
   }
}

/*****************************************************************
 *
 *  Compute the number of states of a MarkovOutTreeData
 *  if unknown
 *
 **/

void MarkovOutTreeData::nb_state_computation()
{
   int t, res = 0;
   Typed_edge_one_int_tree::vertex_iterator it, end;
   value v;

   for(t = 0; t < _nb_trees; t++)
   {
        Tree_tie::tie(it, end) = this->state_trees[t]->vertices();
        while (it < end)
           res = max(res, (this->state_trees[t]->get(*it++)).Int());
   }
   _nb_states = res+1;
}

/*****************************************************************
 *
 *  Compute min and max values variables in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::min_max_value_computation()
{
   if (trees == NULL)
   {
      _min_value.reset(0,0);
      _max_value.reset(0,0);
   }
   else
      Trees::min_max_value_computation();
}

/*****************************************************************
 *
 *  Compute the max size of trees in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::max_size_computation()
{
   _max_size = 0;
   int t;

   if (state_trees != NULL)
   {
      for(t = 0; t < _nb_trees; t++)
         if (state_trees[t] != NULL)
            _max_size = max(_max_size, state_trees[t]->get_size());
   }
   else
      Trees::max_size_computation();
}

/*****************************************************************
 *
 *  Compute the max depth of trees in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::max_depth_computation()
{
   _max_depth = 0;
   int t;

   if (state_trees != NULL)
   {
      for(t = 0; t < _nb_trees; t++)
         if (state_trees[t] != NULL)
            _max_depth = max(_max_depth, state_trees[t]->get_depth());
   }
   else
      Trees::max_depth_computation();
}

/*****************************************************************
 *
 *  Compute the maximal number of children in MarkovOutTreeData
 *
 **/

int MarkovOutTreeData::max_nb_children_computation()
{
   typedef state_tree_type::vertex_iterator vertex_iterator;

   int t;
   unsigned int max_nb_children = 0;
   vertex_iterator it, end;

   if (state_trees != NULL)
   {
      for(t = 0; t < _nb_trees; t++)
         if (state_trees[t] != NULL)
         {
            Tree_tie::tie(it, end) = state_trees[t]->vertices();
            while(it < end)
               max_nb_children= max(max_nb_children,
                                    state_trees[t]->get_nb_children(*it++));
         }
   }
   else
      max_nb_children = Trees::max_nb_children_computation();
   return max_nb_children;
}

/*****************************************************************
 *
 *  Compute the total tree size in MarkovOutTreeData
 *
 **/

int MarkovOutTreeData::cumul_size_computation() const
{  // computation of the cumulated tree sizes
   // i.e. the total number of vertices
   int res = 0;
   register int t;

   if (state_trees != NULL)
      for(t = 0; t < _nb_trees; t++)
         res += state_trees[t]->get_size();
   else
      res = Trees::cumul_size_computation();

   return res;
}

/*****************************************************************
 *
 *  Compute the FrequencyDistribution for total tree size
 *  in MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::build_size_frequency_distribution()
{  // computation of the size frequency distribution
   register int t;

   if (state_trees != NULL)
   {
      if (hsize != NULL)
         delete hsize;

      hsize = new FrequencyDistribution(_max_size+1);
      hsize->nb_element = _nb_trees;
      for(t = 0; t < _nb_trees; t++)
         (hsize->frequency[state_trees[t]->get_size()])++;

      hsize->nb_value_computation();
      hsize->offset_computation();
      hsize->max_computation();
      hsize->mean_computation();
      hsize->variance_computation();
   }
  else
     Trees::build_size_frequency_distribution();
}

/*****************************************************************
 *
 *  Compute the total number of children number in MarkovOutTreeData
 *
 **/

int MarkovOutTreeData::cumul_nb_children_computation() const
{
   // computation of the cumulated number of children of each node
   int res = 0;
   register int t;
   vertex_iterator it, end;

   if (state_trees != NULL)
   {
      for(t = 0; t < _nb_trees; t++)
         if (state_trees[t] != NULL)
         {
            Tree_tie::tie(it, end) = state_trees[t]->vertices();
            while (it < end)
               res += state_trees[t]->get_nb_children(*it++);
         }
   }
   else
      res = Trees::cumul_nb_children_computation();
   return res;
}

/*****************************************************************
 *
 *  Compute the FrequencyDistributions of the
 *  number of children number in MarkovOutTreeData,
 *  taking into account virtual vertices
 *
 **/

void MarkovOutTreeData::build_nb_children_frequency_distribution()
{  // computation of frequency distribution for the number of children
   register int t;
   int nb_vertices = cumul_size_computation();
   vertex_iterator it, end;

   if (state_trees != NULL)
   {
      if (hnb_children != NULL)
         delete hnb_children;

      hnb_children = new FrequencyDistribution(max_nb_children_computation()+1);
      hnb_children->nb_element = nb_vertices;
      for(t = 0; t < _nb_trees; t++)
         if (state_trees[t] != NULL)
         {
            Tree_tie::tie(it, end) = state_trees[t]->vertices();
            while (it < end)
            {
               if (!(this->is_virtual(t, *it)))
                  (hnb_children->frequency[state_trees[t]->get_nb_children(*it)])++;
               it++;
            }
         }

      hnb_children->nb_value_computation();
      hnb_children->offset_computation();
      hnb_children->max_computation();
      hnb_children->mean_computation();
      hnb_children->variance_computation();
   }
   else
      Trees::build_nb_children_frequency_distribution();
}

/*****************************************************************
 *
 *  Print MarkovOutTreeData using an output stream,
 *  a flag on the detail level and one on the comment level
 *
 **/

ostream& MarkovOutTreeData::ascii_write(ostream& os, bool exhaustive,
                                        bool comment_flag) const
{
   if (state_trees != NULL)
   {
      os << 1 << " " << STAT_word[STATW_STATE] << " PROCESS";

      os << "   ";
      if (comment_flag)
         os << "# ";

      os << "(" << _nb_states << " " << STAT_word[STATW_STATES]
         << ")" << endl;

      os << "\n";
      if (comment_flag)
         os << "# ";

      os << STAT_label[STATL_STATE] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";

      if ((state_characteristics != NULL)
           && (state_characteristics->marginal_distribution != NULL))
      {
         state_characteristics->marginal_distribution->ascii_characteristic_print(os, exhaustive, comment_flag);

         if ((state_characteristics->_max_value <= ASCII_NB_VALUE) || (exhaustive))
         {
            os << "\n";
            if (comment_flag)
               os << "# ";

            os << "   | " << STAT_label[STATL_STATE] << " "
               << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
            state_characteristics->marginal_distribution->ascii_print(os, comment_flag);
            os << "\n";
         }
      }
   }

   this->Trees::ascii_write(os, exhaustive, comment_flag);
}


/*****************************************************************
 *
 *  Return the number of views (i.e. the size) in Matplotlib output
 *  due to generation processes in MarkovOutTreeData
 *
 **/

unsigned int MarkovOutTreeData::generation_nb_plot_set_computation() const
{
   if (mot_reestimation == NULL)
      return 0;
   else
      return mot_reestimation->get_nb_generation_values() * _nb_states;
}

/*****************************************************************
 *
 *  Write views due to generation processes into a MultiPlotSet object
 *  for MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::generation_plotable_write(MultiPlotSet &plot_set,
                                                  int &index) const
{
   // plot_set[state][view]
   register int j, k, p;
   unsigned int i, var;
   const int nb_plot_set = generation_nb_plot_set_computation();
   unsigned int nb_generation, nb_factors, nb_children_branching;
   ostringstream title, legend;
   index_array *decodev = NULL; // configuration of factors
   Statiskit::Marginal::Multivariate::CountTable **histo = NULL;
   DiscreteMultivariateReestimation<int> *chisto = NULL;
   MultiPlotSet *generation_plot = NULL;

   title.str("");
   title << STAT_TREES_label[TREESTATL_GENERATION_PROCESS] << ": ";

   if (nb_plot_set > 0)
      histo = mot_reestimation->get_generation_reestim_ptr();
   if (histo != NULL)
   {
      nb_generation = mot_reestimation->get_nb_generation();
      nb_factors = mot_reestimation->get_nb_factors();
      nb_children_branching = mot_reestimation->get_nb_children_branching();

      title << " " << STAT_label[STATL_FIT];

      plot_set.title = title.str();

      decodev = new index_array(mot_reestimation->get_nb_generation(), 0);

      i = 0;
      while (i < nb_plot_set)
      {
         plot_set.variable[index] = 0; // I_DEFAULT ?
         plot_set.viewpoint[index] = 0; // I_DEFAULT ? STATE_PROCESS ?

         title.str("");
         k = i / _nb_states; // generation_value (i.e. configuration of factors)
         // (*plot_set)[i].resize(nb_state);
         j = i % _nb_states; // unordered children state

         // should be called every nb_state only
         if (j == 0)
            if (histo[k] != NULL)
            {
               chisto = stat_histogram_to_DiscreteMultivariateReestimation(*histo[k]);
               generation_plot = chisto->get_plotable();
               delete chisto;
               chisto = NULL;
            }

         plot_set[index].resize((*generation_plot)[j].size());
         for(p = 0; p < (*generation_plot)[j].size(); p++)
            plot_set[index][p] = (*generation_plot)[j][p];
         legend.str("");
         legend << STAT_label[STATL_STATE] << " " << j << " "
                << STAT_TREES_label[TREESTATL_GENERATION]
                << " " << STAT_label[STATL_DISTRIBUTION] << ": ";
         title << legend.str();
         // (*plot_set)[i].legend = legend.str();

         // recode k using factor_base
         *decodev = mot_reestimation->decode(k);
         j = 0;
         if (nb_generation == nb_children_branching + nb_factors + 1)
         {
            // parent state is a factor
            title << STAT_TREES_label[TREESTATL_PARENT] << " " << STAT_label[STATL_STATE] << " "
                  << (*decodev)[j] << " ";
             j++;
          }
          k = 0;
          while (k < nb_children_branching)
          {
             title << STAT_TREES_word[TREESTATW_CHILD] << " " << k + 1 << " "
                   << STAT_label[STATL_STATE] << " " << (*decodev)[j] << " ";
             k++;
             j++;
          }
          k = 0;
          while (k < nb_factors)
          {
             title << STAT_TREES_word[TREESTATW_FACTOR] << " " << k << " "
                   << STAT_TREES_word[TREESTATW_VALUE] << " " << (*decodev)[j] << " ";
             k++;
             j++;
          }
          // title << STAT_label[STATL_DISTRIBUTION];
          plot_set[index].title = title.str();
          if (j == 0)
             if ((histo != NULL) && (histo[k] != NULL))
             {
                delete generation_plot;
                generation_plot = NULL;
             }
          i++;
          index++;
       } // end while i


      delete decodev;
      decodev = NULL;
   }
}

/*****************************************************************
 *
 *  Write views due to discrete output processes into a MultiPlotSet object
 *  for MarkovOutTreeData, for a given variable
 *
 **/

void MarkovOutTreeData::discrete_output_plotable_write(MultiPlotSet &plot_set,
                                                       int &index,
                                                       int var) const

{
   // plot_set[state][view]
   unsigned int val;
   double scale;
   ostringstream title, legend;

   plot_set.variable_nb_viewpoint[var] = 1;

   if (observation_distribution[var] != NULL)
   {
      for(val = 0; val < _nb_states; val++)
      {
         plot_set.variable[index] = var + 1;
         plot_set.viewpoint[index] = OBSERVATION;

         title.str("");
         title << STAT_label[STATL_OUTPUT_PROCESS] << " " << var;
         plot_set[index].title = title.str();

         if (observation_distribution[var] != NULL)
         {
            plot_set[index].xrange = Range(0 , observation_distribution[var][val]->nb_value - 1);
            if (observation_distribution[var][val]->nb_value - 1 < TIC_THRESHOLD)
               plot_set[index].xtics = 1;

            plot_set[index].resize(1);

            if (observation_distribution[var][val]->nb_element > 0)
            {
               plot_set[index].yrange = Range(0, observation_distribution[var][val]->max * YSCALE);
               legend.str("");
               legend << STAT_label[STATL_STATE] << " " << val << " "
                      << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
               plot_set[index][0].legend = legend.str();
               plot_set[index][0].style = "impulses";

               observation_distribution[var][val]->plotable_frequency_write(plot_set[index][0]);
            }
            else
            {
               scale = 1.;
               plot_set[index].yrange = Range(0, 1.);
            }
         }
         index++;
      }
   }

}


/*****************************************************************
 *
 *  Compute the FrequencyDistributions associated with
 *  the output distributions of a MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::observation_frequency_distribution_computation(int ivariable)
{
   register int t, s;
   int val;
   Typed_edge_one_int_tree::vertex_iterator it, end;

   assert(ivariable < _nb_integral);

   for(t = 0; t < _nb_trees; t++)
      if (state_trees[t] != NULL)
      {
         Tree_tie::tie(it, end) = state_trees[t]->vertices();
         while (it < end)
         {
            s = (trees[t]->get(*it)).Int(ivariable);
            s = (this->state_trees[t]->get(*it)).Int();
            val = (trees[t]->get(*it)).Int(ivariable);
            observation_distribution[ivariable][s]->frequency[val]++;
            it++;
         }
      }

   for(s = 0; s < state_characteristics->marginal_distribution->nb_value; s++)
   {
      observation_distribution[ivariable][s]->nb_value_computation();
      observation_distribution[ivariable][s]->offset_computation();
      observation_distribution[ivariable][s]->nb_element_computation();
      observation_distribution[ivariable][s]->max_computation();
      observation_distribution[ivariable][s]->mean_computation();
      observation_distribution[ivariable][s]->variance_computation();
   }
}

/*****************************************************************
 *
 * Compute the characteristic quantity FrequencyDistributions
 * for the state variable of a MarkovOutTreeData
 *
 **/

void MarkovOutTreeData::build_state_characteristics()
{

   Typed_edge_one_int_tree **otrees1 = new Typed_edge_one_int_tree*[_nb_trees];
   int t; // i

   if (_nb_states == I_DEFAULT)
      nb_state_computation();

   if (state_characteristics != NULL)
   {
      delete state_characteristics;
      state_characteristics = NULL;
   }
   for(t = 0; t < _nb_trees; t++)
      otrees1[t] = this->state_trees[t];

   state_characteristics = new TreeCharacteristics(0, _nb_states-1, _max_size,
                                                   _max_depth, _nb_trees, otrees1,
                                                   0, true);

   delete [] otrees1;
   otrees1 = NULL;
}

/*****************************************************************
 *
 *  Compute the FrequencyDistributions of initial states
 *
 **/

void MarkovOutTreeData::build_chain_reestimation(char itype, int inb_state)
{
   unsigned int s, t;

   if (chain_data == NULL)
      chain_data = new ChainData(itype, inb_state, inb_state);

   for(s = 0; s < inb_state; s++)
      chain_data->initial[s] = 0;

   for(t = 0; t < _nb_trees; t++)
      if (state_trees[t] != NULL)
      {
         s = (state_trees[t]->get(state_trees[t]->root())).Int();
         chain_data->initial[s]++;
      }
}

/*****************************************************************
 *
 *  Compute the FrequencyDistributions of the joint distribution
 *  of children states in a MarkovOutTreeData using the number
 *  of states, the number of children states, of factors
 *  and parent children determining the number inb_generation
 *  of branching processes. Compute the sequences of ordered children.
 *  build_nb_children_frequency_distribution and min_max_value_computation
 *  must be called before (or min_offset and hnb_children must be up-to-date)
 *
 **/

void MarkovOutTreeData::build_markov_reestimation(unsigned int inb_children_branching,
                                                  unsigned int inb_factors,
                                                  const GenerationProcess::index_array *ifactor_values,
                                                  unsigned int inb_generation,
                                                  unsigned int iunordered_rank,
                                                  unsigned int inb_vomc,
                                                  const std::vector<unsigned int>* factor_indices)
{
   bool init_from_markov = true; // initialize bounds from Markov generation distributions
   unsigned int s;
   std::vector<unsigned int> min_offset, max_nb_values;

   if (mot_reestimation == NULL)
      mot_reestimation = new MarkovOutTreeReestimation<int>(_nb_states,
                                                            inb_children_branching,
                                                            inb_factors,
                                                            ifactor_values,
                                                            inb_generation,
                                                            iunordered_rank,
                                                            inb_vomc);
   if (markov != NULL)
      min_offset = markov->get_generation_min_offsets();
   else
      init_from_markov = false;
   if (!init_from_markov)
   {
      min_offset.resize(1);
      min_offset[0] = hnb_children->offset;
   }
   max_nb_values.resize(min_offset.size());
   for(s = 0; s < max_nb_values.size(); s++)
      max_nb_values[s] = hnb_children->nb_value + 1;

   mot_reestimation->init(min_offset, max_nb_values);
   mot_reestimation->build_characteristics(*this, I_DEFAULT, factor_indices);

}

/*****************************************************************
 *
 *  Permute the states of a MarkovOutTreeData
 *  based on a given permutation perm.
 *  Validity of permutation must be ensured before call
 *
 **/

void MarkovOutTreeData::state_permutation(int* perm)
{
   register int i, var, t;
   TreeCharacteristics* pstate_char = new TreeCharacteristics[_nb_states];
   ptFrequencyDistribution_array_2d pobservation_dist = new FrequencyDistribution**[_nb_integral];
   ptHistogram_array_2d pobservation_histo = new Histogram**[_nb_integral];
   Typed_edge_one_int_tree::vertex_iterator it, end;
   One_int_container c;

   for(var = 0; var < _nb_integral; var++)
   {
      pobservation_dist[var] = new FrequencyDistribution*[_nb_states];
      pobservation_histo[var] = new Histogram*[_nb_states];
      for(i = 0; i < _nb_states; i++)
      {
         pobservation_dist[var][perm[i]] = observation_distribution[var][i];
         pobservation_histo[var][perm[i]] = observation_histogram[var][i];
      }
      for(i = 0; i < _nb_states; i++)
      {
         observation_distribution[var][i]= pobservation_dist[var][i];
         observation_histogram[var][i]= pobservation_histo[var][i];
      }
      delete [] pobservation_dist[var];
      pobservation_dist[var] = NULL;
      delete [] pobservation_histo[var];
      pobservation_histo[var] = NULL;
   }

   delete [] pstate_char;
   pstate_char = NULL;
   delete [] pobservation_dist;
   pobservation_dist = NULL;
   delete [] pobservation_histo;
   pobservation_histo = NULL;

   if (state_trees != NULL)
   {
      for(t= 0; t < _nb_trees; t++)
      {
         Tree_tie::tie(it, end)= state_trees[t]->vertices();
         while (it < end)
         {
            c = state_trees[t]->get(*it);
            c.Int() = perm[c.Int()];
            state_trees[t]->put(*it++, c);
         }
      }
   }
   build_state_characteristics();
   // vomc_data and mot_reestimation should be permuted as well
}

/*****************************************************************
 *
 *  Select subset of state trees from their indices in MarkovOutTreeData
 *  using a StatError object, the number of trees, a tree identifier
 *  and a flag on keeping or rejecting the selected trees.
 *  Mapping between new and old identifiers is stored in index.
 *
 **/


MarkovOutTreeData::ptOne_int_tree_set
MarkovOutTreeData::state_trees_select_individual(StatError& error,
                                                 int inb_trees,
                                                 int_array iidentifier,
                                                 int_array& index,
                                                 bool keep) const

{
   bool status = true, *selected_trees;
   register int t, j;
   int max_identifier;
   int_array identifier;
   ptOne_int_tree_set res = NULL;

   error.init();

   if ((inb_trees < 1) || (inb_trees > (keep ? _nb_trees : _nb_trees-1))
       || (state_trees == NULL))
   {
      status = false;
      error.update(STAT_TREES_error[TREESTATR_NB_TREES]);
   }
   else
   {
      max_identifier = 1;
      for(t = 0; t < inb_trees; t++)
         if (iidentifier[t] > max_identifier)
            max_identifier = iidentifier[t];

      selected_trees = new bool[max_identifier+1];
      identifier = new int[max_identifier+1];
      for(t = 0; t <= max_identifier; t++)
      {
         selected_trees[t] = false;
         identifier[t] = t;
      }
      for(t = 0; t < inb_trees; t++)
      {
         for(j = 0; j < _nb_trees; j++)
             if (iidentifier[t] == identifier[j])
                break;

         if (j == _nb_trees)
         {
            status = false;
            ostringstream error_message;
            error_message << iidentifier[t] << ": " << STAT_TREES_error[TREESTATR_TREE_IDENTIFIER];
            error.update((error_message.str()).c_str());
         }
         else
            if (selected_trees[iidentifier[t]])
            {
               status = false;
               ostringstream error_message;
               error_message << STAT_TREES_label[TREESTATL_TREE] << " " << iidentifier[t] << " "
                             << STAT_error[STATR_ALREADY_SELECTED];
               error.update((error_message.str()).c_str());
            }
            else
               selected_trees[iidentifier[t]]= true;
      }
      delete [] selected_trees;
   }

   if (status)
   {
      index = identifier_select(_nb_trees, identifier, inb_trees, iidentifier, keep);
      res = new Typed_edge_one_int_tree*[inb_trees];
      for(t = 0; t < inb_trees; t++)
         res[t] = new Typed_edge_one_int_tree(*(state_trees[index[t]]));
   }
   return res;
}

/*****************************************************************
 *
 *  Left (bit) shift operator of a MarkovOutTreeData
 *
 **/

std::ostream& Stat_trees::operator<<(std::ostream& os, const MarkovOutTreeData& data)
{ return data.ascii_write(os, false); }
