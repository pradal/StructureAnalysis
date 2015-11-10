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
 *       $Id: multivariate_distributions.cpp 9554 2010-09-20 17:10:30Z jbdurand $
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

#include <numeric>
#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"   // definition of DiscreteParametricModel class
#include "stat_tool/markovian.h"   // definition of DiscreteParametricProcess class
#include "stat_tool/stat_label.h"   // definition of DiscreteParametricProcess class
                          // and D_DEFAULT_PARAMETER replaced by D_DEFAULT
#include "stat_tool/vectors.h"

#include "tree_labels.h"
#include "statiskit/core/data/design.h"
#include "statiskit/core/data/marginal/multivariate.h"
#include "mv_histogram_tools.h"

#include "multivariate_distribution.h"   // definition of DiscreteParametricModel class

#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include <assert.h>

extern char* label(const char*);

using namespace stat_tool;
using namespace Stat_trees;


/*****************************************************************
 *
 *  Default constructor of DiscreteMultivariateDistribution class
 *
 **/

DiscreteMultivariateDistribution::DiscreteMultivariateDistribution()
 : nb_value()
 , offset()
 , mean()
 , variance()
 , mass(NULL)
 , marginals()
 , nb_variable(0)
 , total_nb_value(0)
 , complement(0.)
 , nb_parameter(I_DEFAULT)
 , val_strides(NULL)
 , counters(NULL)
{
   counters = new index_array;
   val_strides = new index_array;
   mass = new valarray<double>(0);
}

/*****************************************************************
 *
 *  Copy constructor of DiscreteMultivariateDistribution class
 *
 **/

DiscreteMultivariateDistribution::DiscreteMultivariateDistribution(const DiscreteMultivariateDistribution& dist)
 : nb_value()
 , offset()
 , mean()
 , variance()
 , mass(NULL)
 , marginals()
 , nb_variable(0)
 , total_nb_value(0)
 , complement(0.)
 , nb_parameter(I_DEFAULT)
 , val_strides(NULL)
 , counters(NULL)
{ copy(dist); }

/*****************************************************************
 *
 *  Copy a DiscreteMultivariateDistribution and return a pointer
 *
 **/

DiscreteMultivariateDistribution* DiscreteMultivariateDistribution::copy_ptr() const
{
   DiscreteMultivariateDistribution *res = new DiscreteMultivariateDistribution(*this);

   return res;
}

/*****************************************************************
 *
 *  Destructor for DiscreteMultivariateDistribution class
 *
 **/

DiscreteMultivariateDistribution::~DiscreteMultivariateDistribution()
{ remove(); }

/*****************************************************************
 *
 *  Return the probability associated with a set of indices
 *  (or D_DEFAULT if the value is not stored in *mass)
 *
 **/

double DiscreteMultivariateDistribution::get_mass(const index_array& indices) const
{
   bool return_default = false;
   unsigned int i, mono_index = 0; // index in mass

   assert((indices.size() == nb_variable) && (total_nb_value > 0));
   assert(mass != NULL);

   if (val_strides->empty())
      val_strides_computation();

   for(i = 0; i < nb_variable; i++)
   {
      if (indices[i] >= nb_value[i])
         return_default = true;
      mono_index += (*val_strides)[i] * indices[i];
   }

   if(mono_index >= total_nb_value)
      return_default = false;

   if ((mono_index < mass->size()) && (!return_default))
      return (*mass)[mono_index];
   else
      return D_DEFAULT;
}

/*****************************************************************
 *
 *  Get either the probability value for a set of indices,
 *  or the log probability (D_DEFAULT if the value is not stored in *mass)
 *
 **/

double DiscreteMultivariateDistribution::get_mass(const index_array& indices, bool log_computation) const
{
   double res = get_mass(indices);

   if ((res > 0) && (log_computation))
      return log(res);
   else
      return res;
}

/*****************************************************************
 *
 *  Set the probability associated with a set of indices
 *  and return true iif the value could actually be set
 *
 **/

bool DiscreteMultivariateDistribution::set_mass(const index_array& indices,
                                                double value) const
{
   bool set = true;
   unsigned int i, mono_index = 0; // index in mass
   const unsigned int eff_nb_value = min(total_nb_value, MULTIVARIATE_DISTRIBUTION_MAX_VALUE+1);

   assert((indices.size() == nb_variable) && (total_nb_value > 0));
   assert(mass != NULL);

   if (mass->size() != eff_nb_value)
      mass->resize(eff_nb_value, D_DEFAULT);

   if (val_strides->empty())
      val_strides_computation();

   for(i = 0; i < nb_variable; i++)
   {
      if (indices[i] >= nb_value[i])
         set = false;
      mono_index += (*val_strides)[i] * indices[i];
   }

   if ((mono_index < eff_nb_value) && (set))
      (*mass)[mono_index] = value;

   return set;

}

/*****************************************************************
 *
 *  Return marginal distribution of a given variable
 *  in DiscreteMultivariateDistribution
 *  (a new instance is allocated)
 *
 **/

Distribution* DiscreteMultivariateDistribution::get_marginal(unsigned int ivariable) const
{
   assert((ivariable < nb_variable) && (marginals.size() == nb_variable));

   if (marginals[ivariable] == NULL)
      return NULL;
   else
      return new Distribution(*marginals[ivariable]);
}

/*****************************************************************
 *
 *  Init a loop on indices
 *
 **/

DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_init() const
{
   unsigned int i;

   assert(counters != NULL);

   if (nb_variable > 0)
      counters->resize(nb_variable);

   for(i = 0; i < nb_variable; i++)
      (*counters)[i] = 0;

   return *counters;
}

/*****************************************************************
 *
 *  Init a loop on indices, starting at a given index
 *
 **/

DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_init(const index_array& begin_index) const
{
   unsigned int i;

   assert(counters != NULL);
   assert(begin_index.size() == nb_variable);

   if (nb_variable > 0)
      counters->resize(nb_variable);

   for(i = 0; i < nb_variable; i++)
   {
      assert(begin_index[i] < nb_value[i]);
      (*counters)[i] = begin_index[i];
   }

   return *counters;
}

/*****************************************************************
 *
 *  Loop on indices (one step)
 *
 **/

DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_next() const
{
   bool update = false; // should some indices be reset ?
   unsigned int i = nb_variable - 1;

   assert(counters != NULL);

   if ((counters->size() == nb_variable) && (nb_value.size() == nb_variable))
   {
      // update last index
      (*counters)[i]++;
      if ((*counters)[i] == nb_value[i])
      {
         // update other indices if the maximum is reached
         update = true;
         (*counters)[i] = 0;
      }
      while (update)
      {
         // update every index that reached its maximum
         i--;
         (*counters)[i]++;
         if ((*counters)[i] == nb_value[i])
            (*counters)[i] = 0;
         else
            update = false;
         if (i == 0)
            update = false;
      }
   }
   else
      counters->resize(0);

   return *counters;
}

/*****************************************************************
 *
 *  Init a loop on indices (lowest indices first)
 *
 **/

/*DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_init_low() const
{
   /// first element of counters represents the current maximal index valu_e
   unsigned int i;

   assert(counters != NULL);

   if (nb_variable > 0)
      counters->resize(nb_variable+1);

   for(i = 0; i < nb_variable+1; i++)
      (*counters)[i] = 0;

   return *counters;
}*/

/*****************************************************************
 *
 *  Init a loop on indices, starting at a given index
 *  (lowest indices first)
 *
 **/

/*DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_init_low(const index_array& begin_index) const
{
   unsigned int i, max_index = 0;
   index_array res(nb_variable);

   assert(counters != NULL);
   assert(begin_index.size() == nb_variable);

   if (nb_variable > 0)
      counters->resize(nb_variable+1);

   for(i = 0; i < nb_variable; i++)
   {
      assert(begin_index[i] < nb_value[i]);
      (*counters)[i+1] = begin_index[i];
      res[i] = begin_index[i];
      max_index = max(max_index, begin_index[i]);
   }
   (*counters)[0] = max_index;

   return res;
}*/

/*****************************************************************
 *
 *  Loop on indices (one step, lowest indices first)
 *
 **/
DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_next_low() const
{
   unsigned int i;
   index_array min_array(nb_variable);

   for(i = 0; i < nb_variable; i++)
      min_array[i] = 0;

   return loop_next_low(min_array);
}

/*
DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_next_low() const
{
   bool update = false, back_update = false, // should some indices be reset ?
        first_hit = true, // first time a new value is hit
        last_value = true; // current vector is nb_value - 1
   unsigned int i, j,
                vmax = 0, ifmax, ilmax, // maximum index value and its first and last positions
                fn0; // first non-null index

   assert(counters != NULL);

   if ((counters->size() == nb_variable) && (nb_value.size() == nb_variable))
   {
      // find maximal indices, or those at their highest possible value
      ifmax = 0;
      for(i = 0; i < nb_variable; i++)
      {
         if (((*counters)[i] == vmax) || ((*counters)[i] == nb_value[i]-1))
            ilmax = i;

         if ((*counters)[i] > vmax)
         {
            ifmax = i;
            ilmax = i;
            vmax = (*counters)[i];
         }
      }

      for(i = 0; i < nb_variable; i++)
         first_hit = first_hit && ((*counters)[i] == vmax);

      for(i = 0; i < nb_variable; i++)
         last_value = last_value && ((*counters)[i] == nb_value[i]-1);

      if (last_value)
      {
         for(i = 0; i < nb_variable; i++)
            (*counters)[i] = 0;
         return *counters;
      }

      // increase counters for every index after ifmax
      if ((ifmax == nb_variable - 1) && (nb_variable > 2) && (ilmax == nb_variable - 1))
         // vmax is at nb_variable - 1 only: increase at nb_variable - 2
         i = nb_variable - 2;
      else
         i = nb_variable - 1;
      if (first_hit)
      {
         (*counters)[i]++;
         for(j = 0; j < i; j++)
            (*counters)[j] = 0;
      }
      // if (((*counters)[i] >= vmax) && !(first_hit)) // not the first time vmax is reached
      if (!first_hit) // not the first time vmax is reached
      {
         if (i == nb_variable - 2)
         {
            (*counters)[i]++;
            // find first non-null index
            fn0 = 0;
            first_hit = ((*counters)[fn0] == 0);
            while (first_hit)
            {
               fn0++;
               if ((*counters)[fn0] > 0)
                  first_hit = false;
            }
            // if (((*counters)[i] == vmax) && (i == fn0))
            if (((*counters)[i] == vmax))
               (*counters)[nb_variable-1] = 0;
            // should be min(vmax, nb_values[nb_variable-2]) and
            // then search for the lowest index that can be at vmax
         }
         else
         {
            (*counters)[i]++;
            back_update = false;
            if ((*counters)[i] > min(vmax, nb_value[i]-1))
            {
               update = true;
               back_update = true;
               (*counters)[i] = 0;
            }
            while ((update) || (back_update))
            {
               // update every index that reached its maximum
               i--;

               if (update)
               {
                  (*counters)[i]++;

                  if ((*counters)[i] <= min(vmax, nb_value[i]-1))
                    update = false;
                  else
                    (*counters)[i] = 0;
               }

               if ((*counters)[i] == vmax)
                  back_update = false;

               if (i == 0)
               {
                  update = false;
                  back_update = false;
               }
            }

            if (i == 0)
               back_update = ((*counters)[i] != vmax);

            if (back_update)
               (*counters)[nb_variable-1] = vmax;
               // should be min(vmax, nb_values[nb_variable-1]) and
               // then search for the highest index that can be at vmax

         }
      }
   }
   else
      counters->resize(0);

   return *counters;
}
*/

/*****************************************************************
 *
 *  Loop on indices (one step, lowest indices first), respecting
 *  some given min indices
 *
 **/

DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::loop_next_low(const index_array& min_index) const
{
   bool update = false, back_update = false, // should some indices be reset ?
        first_hit = true, // first time a new value is hit
        last_value = true; // current vector is nb_value - 1
   unsigned int i, j,
                vmax = 0, ifmax, ilmax, // maximum index value and its first and last positions
                fn0; // first non-null index

   assert(counters != NULL);

   if ((counters->size() == nb_variable) && (nb_value.size() == nb_variable))
   {
      for(i = 0; i < nb_variable; i++)
         assert((*counters)[i] >= min_index[i]);

      // find maximal indices, or those at their highest possible value
      ifmax = 0;
      for(i = 0; i < nb_variable; i++)
      {
         if (((*counters)[i] == vmax) || ((*counters)[i] == nb_value[i]-1))
            ilmax = i;

         if ((*counters)[i] > vmax)
         {
            ifmax = i;
            ilmax = i;
            vmax = (*counters)[i];
         }
      }

      for(i = 0; i < nb_variable; i++)
         first_hit = first_hit && ((*counters)[i] == vmax);

      for(i = 0; i < nb_variable; i++)
         last_value = last_value && ((*counters)[i] == nb_value[i]-1);

      if (last_value)
      {
         for(i = 0; i < nb_variable; i++)
            (*counters)[i] = 0;
         return *counters;
      }

      // increase counters for every index after ifmax
      if ((ifmax == nb_variable - 1) && (nb_variable > 2) && (ilmax == nb_variable - 1))
         // vmax is at nb_variable - 1 only: increase at nb_variable - 2
         i = nb_variable - 2;
      else
         i = nb_variable - 1;
      if (first_hit)
      {
         (*counters)[i]++;
         for(j = 0; j < i; j++)
            (*counters)[j] = min_index[j];
      }

      if (!first_hit) // not the first time vmax is reached
      {
         if (i == nb_variable - 2)
         {
            (*counters)[i]++;
            // find first non-minimal index
            fn0 = 0;
            first_hit = ((*counters)[fn0] == min_index[fn0]);
            while (first_hit)
            {
               fn0++;
               if ((*counters)[fn0] > min_index[fn0])
                  first_hit = false;
            }
            if (((*counters)[i] == vmax))
               (*counters)[nb_variable-1] = min_index[nb_variable-1];
            // should be min(vmax, nb_values[nb_variable-2], min_index[nb_variable-1]) and
            // then search for the lowest index that can be at vmax
         }
         else
         {
            (*counters)[i]++;
            back_update = false;
            if ((*counters)[i] > min(vmax, nb_value[i]-1))
            {
               update = true;
               back_update = true;
               (*counters)[i] = min_index[i];
            }
            while ((update) || (back_update))
            {
               // update every index that reached its maximum
               i--;

               if (update)
               {
                  (*counters)[i]++;

                  if ((*counters)[i] <= min(vmax, nb_value[i]-1))
                    update = false;
                  else
                    (*counters)[i] = min_index[i];
               }

               if ((*counters)[i] == vmax)
                  back_update = false;

               if (i == 0)
               {
                  update = false;
                  back_update = false;
               }
            }

            if (i == 0)
               back_update = ((*counters)[i] != vmax);

            if (back_update)
               (*counters)[nb_variable-1] = vmax;
               // should be min(vmax, nb_values[nb_variable-1]) and
               // then search for the highest index that can be at vmax

         }
      }
   }
   else
      counters->resize(0);

   return *counters;
}

/*****************************************************************
 *
 *  Destructor for DiscreteMultivariateDistribution class
 *
 **/

void DiscreteMultivariateDistribution::remove()
{
   unsigned int i;

   for(i = 0; i < marginals.size(); i++)
      if (marginals[i] != NULL)
      {
         delete marginals[i];
         marginals[i] = NULL;
      }
   if (counters != NULL)
   {
      delete counters;
      counters = NULL;
   }
   if (mass != NULL)
   {
      delete mass;
      mass = NULL;
   }

   if (val_strides != NULL)
   {
      delete val_strides;
      val_strides = NULL;
   }
}

/*****************************************************************
 *
 *  Copy operator for DiscreteMultivariateDistribution class
 *
 **/

void DiscreteMultivariateDistribution::copy(const DiscreteMultivariateDistribution& dist)
{
   unsigned int i;

   nb_variable = dist.nb_variable;
   nb_value = dist.nb_value;
   total_nb_value = dist.total_nb_value;
   offset = dist.offset;
   complement = dist.complement;
   mean = dist.mean;
   variance = dist.variance;

   marginals.resize(dist.marginals.size());
   for(i = 0; i < marginals.size(); i++)
      marginals[i] = new Distribution(*dist.marginals[i]);

   nb_parameter = dist.nb_parameter;

   if ((mass != NULL) && (dist.mass != NULL))
       *mass = *dist.mass;
   if ((mass == NULL) && (dist.mass != NULL))
       mass = new valarray<double>(*dist.mass);

   if ((val_strides != NULL) && (val_strides != NULL))
       *val_strides = *dist.val_strides;
   if ((val_strides == NULL) && (dist.val_strides != NULL))
       val_strides = new index_array(*dist.val_strides);

   if ((counters != NULL) && (dist.counters != NULL))
       *counters = *dist.counters;
   if ((counters == NULL) && (dist.counters != NULL))
       counters = new index_array(*dist.counters);
}

/*****************************************************************
 *
 *  Compute the val_strides (products for each dimension > i
 *  of the number of values for this dimension).
 *  Useful to manipulate val_arrays
 *
 **/

void DiscreteMultivariateDistribution::val_strides_computation() const
{
   unsigned int i;

   assert(val_strides != NULL);

   val_strides->resize(nb_variable);

   (*val_strides)[nb_variable-1] = 1;

   for(i = nb_variable-1; i > 0; i--)
      (*val_strides)[i-1] = (*val_strides)[i] * nb_value[i-1];
}

/*****************************************************************
 *
 *  Print a DiscreteMultivariateDistribution into an output stream, using
 *  a flag on the level of detail
 *
 **/

ostream& DiscreteMultivariateDistribution::ascii_write(ostream &os, bool exhaustive) const
{ return ascii_print(os, false, exhaustive, false, false, NULL); }

/*****************************************************************
 *
 *  Print the nature of distribution and number of variables
 *
 **/

ostream& DiscreteMultivariateDistribution::header_ascii_print(ostream &os) const
{
   os << STAT_TREES_multivariate_distribution_word[DISCRETE_MULTIVARIATE] << endl;
   os << nb_variable << " " << STAT_word[STATW_VARIABLES] << endl;

   return os;
}

/*****************************************************************
 *
 *  Print the header and the whole set of probabilities
 *  for DiscreteMultivariateDistribution
 *
 **/

ostream& DiscreteMultivariateDistribution::ascii_print(ostream &os, bool comment_flag,
                                                       bool exhaustive, bool file_flag, bool cumul_flag,
                                                       const Statiskit::Marginal::Multivariate::CountTable *histo) const
{
   unsigned int i;
   FrequencyDistribution h, *m = NULL;
   Statiskit::Marginal::Multivariate::CountTable::marginal_type::ptr_type chisto;

   header_ascii_print(os);
   os << endl;
   probability_ascii_print(os, comment_flag, cumul_flag);

   if ((histo != NULL) && (exhaustive) && (marginals.size() == nb_variable))
      for(i = 0; i < nb_variable; i++)
         if (marginals[i] != NULL)
         {
            os << "\n";
            if (comment_flag)
               os << "# ";
            os << "   | " << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_TREES_label[TREESTATL_GENERATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << " | " << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_TREES_label[TREESTATL_GENERATION] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
            chisto = Statiskit::marginalizing(*histo, i);
            m = Scalar_stat_histogram_to_frequency_distribution(*chisto);
            h.copy(*m);
            delete m;
            m = NULL;
            marginals[i]->ascii_print(os, comment_flag, cumul_flag, true,
                                      &h);

         }
}

/*****************************************************************
 *
 *  Print the whole set of probabilities
 *  for DiscreteMultivariateDistribution
 *
 **/

ostream& DiscreteMultivariateDistribution::probability_ascii_print(ostream &os, bool comment_flag,
                                                                   bool cumul_flag) const
{
   unsigned int i, // index in mass
                var;
   index_array indices;

   if (comment_flag)
      os << "# " ;

   for(var = 0; var < nb_variable; var++)
      os << STAT_word[STATW_VARIABLE] << var << " | ";
   os << " " << STAT_word[STATW_PROBABILITY] << endl;

   if (total_nb_value < MULTIVARIATE_DISTRIBUTION_MAX_VALUE)
   {
      indices = loop_init();
      for(i = 0; i < total_nb_value; i++)
      {
         if (comment_flag)
            os << "# " ;
         os << "\t";
         for(var = 0; var < nb_variable; var++)
            os << indices[var] << " | ";
         os << " " << get_mass(indices) << endl;
         //indices = loop_next_low();
         indices = loop_next();
      }
   }
   else
   {
      indices = loop_init(offset);
      for(i = 0; i < MULTIVARIATE_DISTRIBUTION_MAX_VALUE; i++)
      {
         if (comment_flag)
            os << "# " ;
         os << "\t";
         for(var = 0; var < nb_variable; var++)
            os << indices[var] << " | ";
         os << " " << get_mass(indices) << endl;
         indices = loop_next_low();
      }
   }
   return os;
}

/*****************************************************************
 *
 *  Print the whole set of probabilities for DiscreteMultivariateDistribution
 *  in a spreadsheet fashion
 *
 **/

ostream& DiscreteMultivariateDistribution::spreadsheet_print(ostream &os) const
{ return os; }

/*****************************************************************
 *
 *  Read a DiscreteMultivariateDistribution from a file
 *
 **/

DiscreteMultivariateDistribution* Stat_trees::discrete_multivariate_ascii_read(StatError &error, const char *path,
                                                                               double cumul_threshold)
{
   RWCString buffer;
   size_t position;
   bool status;
   int line;
   DiscreteMultivariateDistribution *pdist = NULL;
   ifstream in_file(path);

   error.init();

   if (in_file == NULL)
      error.update(STAT_error[STATR_FILE_NAME]);
   else
   {
      status = true;
      line = 0;

      pdist = discrete_multivariate_parsing(error, in_file, line,
                                            MCOMPOUND_MULTINOMIAL,
                                            cumul_threshold);

      if (pdist == NULL)
         status = false;

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
   }
   return pdist;
}

/*****************************************************************
 *
 *  Matplotlib output for DiscreteMultivariateDistribution class
 *  using potential data
 *
 **/

MultiPlotSet* DiscreteMultivariateDistribution::get_plotable(const Statiskit::Marginal::Multivariate::CountTable* histo) const
{
   MultiPlotSet *plot_set = NULL;

   return plot_set;
}

/*****************************************************************
 *
 *  Variable permutation for DiscreteMultivariateDistribution class
 *  (checking or not validity of permutation)
 *
 **/

void DiscreteMultivariateDistribution::variable_permutation(StatError& error,
                                                            int* perm, bool check_flag)
{
   unsigned int i, j;
   bool status = true;
   // each element of check_perm must be used exactly once

   if (check_flag)
   {
      // check permutation
      bool * const check_perm = new bool[nb_variable];
      error.init();

      for(i = 0; i < nb_variable; i++)
         check_perm[i] = false; // indicates if ith element was used in perm

      for(i = 0; i < nb_variable; i++)
      {
         if (check_perm[perm[i]])
            status = false;
         else
            check_perm[perm[i]] = true;
      }
      if (!status)
         error.update(STAT_TREES_error[TREESTATR_NO_PERMUTATION]);
   }
   if (status)
   {
      assert(nb_variable > 0);

      if (counters != NULL)
      {
          counters->resize(nb_variable);

          for(i = 0; i < nb_variable; i++)
             (*counters)[i] = 0;
      }

      std::vector<unsigned int> pnb_value(nb_variable);
      std::vector<unsigned int> poffset(nb_variable);
      std::vector<double> pmean(nb_variable);
      std::vector<double> pvariance(nb_variable);
      std::vector<Distribution*> pmarginals(nb_variable);

      for(i = 0; i < nb_variable; i++)
      {
         pnb_value[perm[i]] = nb_value[i];
         nb_value[i] = pnb_value[perm[i]];
         poffset[perm[i]] = offset[i];
         offset[i] = poffset[i];
      }
      if (mean.size() == nb_variable)
         for(i = 0; i < nb_variable; i++)
         {
            pmean[perm[i]] = mean[i];
            mean[i] = pmean[i];
         }
      if (variance.size() == nb_variable)
         for(i = 0; i < nb_variable; i++)
         {
            pvariance[perm[i]] = variance[i];
            variance[i] = pvariance[i];
         }

      if (mass != NULL)
      {
         std::valarray<double> pmass(*mass);
         index_array pval_strides = *val_strides;
         index_array indices(nb_variable), pindices(nb_variable);
         // total_nb_value = prod(nb_value[i]) remains unchanged
         val_strides_computation();

         for(i = 0; i < mass->size(); i++)
         {
            // compute old and new indices associated with element i
            indices[0] = i / (*val_strides)[0];
            pindices[perm[0]] = indices[0];
            for(j = 1; j < nb_variable; j++)
            {
               indices[j] = indices[j-1] / (*val_strides)[j];
               pindices[perm[j]] = indices[j];
            }
            // mass associated with old indices:
            // pmass[i]
            set_mass(pindices, pmass[i]);
         }
      }
   }
}

/*****************************************************************
 *
 *  Parse a DiscreteMultivariateDistribution from a file
 *
 **/

DiscreteMultivariateDistribution* Stat_trees::discrete_multivariate_parsing(StatError &error, std::ifstream &in_file,
                                                                            int &line, int last_ident,
                                                                            double cumul_threshold,
                                                                            int min_inf_bound)
{
   unsigned int nb_variable;
   RWLocaleSnapshot locale("en");
   RWCString buffer, token;
   size_t position;
   bool status = true, lstatus;
   register int i, j;
   int ident = I_DEFAULT;
   long inf_bound, sup_bound = I_DEFAULT,
        value;
   // double parameter = D_DEFAULT , probability = D_DEFAULT;
   DiscreteMultivariateDistribution *dist;


   dist = NULL;

   while (buffer.readLine(in_file, false))
   {
      line++;

  #   ifdef DEBUG
      cout << line << " " << buffer << endl;
  #   endif

      position = buffer.first('#');
      if (position != RW_NPOS)
         buffer.remove(position);
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull()))
      {
         // search for DISCRETE_MULTIVARIATE_DISTRIBUTION

         // test nom de la loi

         if (i == 0)
         {
            if (token != STAT_TREES_multivariate_distribution_word[DISCRETE_MULTIVARIATE])
            {
               status = false;
               error.update(STAT_parsing[STATP_DISTRIBUTION_NAME] , line , i + 1);
            }
         }
         i++;
      }
      if (i > 0)
      {
         if (i != 1)
         {
            status= false;
            error.update(STAT_parsing[STATP_FORMAT], line);
         }
         break;
       }
   } // end while

   if (status)
   {
      // read number of variables
      i = 0;
      while (buffer.readLine(in_file, false))
      {
         line++;

     #   ifdef DEBUG
         cout << line << " " << buffer << endl;
     #   endif

         position = buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);

         RWCTokenizer next(buffer);

         while (!((token = next()).isNull()))
         {
            switch (i)
            {
               // test number of variables
               case 0 :
               {
                  lstatus = locale.stringToNum(token, &value);
                  if (lstatus)
                     if (value < 2)
                        lstatus = false;
                     else
                        nb_variable = value;

                  if (!lstatus)
                  {
                     status = false;
                     error.update(STAT_parsing[STATP_NB_VARIABLE] , line , i + 1);
                  }
                  break;
               }

               // test for VARIABLES keyword

               case 1 :
               {
                  if (token != STAT_word[STATW_VARIABLES])
                  {
                     status = false;
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_word[STATW_VARIABLES],
                                             line, i + 1);
                  }
                  break;
               }
            } // end switch
            i++;
         } // end while
         if (i > 0)
         {
            if (i != 2)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT] , line);
            }
            break;
         }
      } // end while
   } // end if



   if (status)
   {
      // read nature of distribution
      i = 0;
      while (buffer.readLine(in_file, false))
      {
         // read new line
         line++;

     #   ifdef DEBUG
         cout << line << " " << buffer << endl;
     #   endif

         position = buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);

         RWCTokenizer next(buffer);

         // parse a single line
         while (!((token = next()).isNull()))
         {
            if (i == 0)
            {
               for (j = DISCRETE_MULTIVARIATE; j <= last_ident; j++)
                  if (token == STAT_TREES_multivariate_distribution_word[j])
                  {
                     ident = j;
                     break;
                  }
            }
            if (i != 0)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            // test nature of distribution
            if (j == MIID)
               dist = iid_discrete_multivariate_parsing(error,
                                                        in_file,
                                                        line,
                                                        nb_variable,
                                                        cumul_threshold,
                                                        min_inf_bound);

            if ((j >= MPOISSON) && (j <= MMULTINOMIAL))
               // the remainder of the line should be sent together
               // with other arguments
               dist = discrete_multivariate_parametric_parsing(error,
                                                               in_file, line,
                                                               next,
                                                               j, nb_variable,
                                                               cumul_threshold,
                                                               min_inf_bound);

            if (j == MCOMPOUND_MULTINOMIAL)
               dist = multinomial_compound_parsing(error, in_file, line,
                                                   next, nb_variable,
                                                   cumul_threshold,
                                                   min_inf_bound);

            if (j == last_ident + 1)
               dist = discrete_multivariate_probability_parsing(error, in_file,
                                                                line, nb_variable,
                                                                cumul_threshold,
                                                                min_inf_bound);
            i++;
         } // end while
         if (i != 1)
         {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
         }
         else
            break;
      } // end while
   } // end if

   return dist;
}

/*****************************************************************
 *
 *  Parse probabilities in DiscreteMultivariateDistribution from a file
 *
 **/

DiscreteMultivariateDistribution* Stat_trees::discrete_multivariate_probability_parsing(StatError &error, std::ifstream &in_file,
                                                                                        int &line,
                                                                                        unsigned int inb_variable,
                                                                                        double cumul_threshold,
                                                                                        int min_inf_bound)
{
   DiscreteMultivariateDistribution *pdist = NULL;

   return pdist;
}

/*****************************************************************
 *
 *  Default constructor of DiscreteMultivariateParametric class
 *
 **/

DiscreteMultivariateParametric::DiscreteMultivariateParametric()
 : DiscreteMultivariateDistribution()
 , rparameter(NULL)
 , sup_bound()
 , parameter()
 , probability()
 , param_marginals()
 , ident(I_DEFAULT)
 , inf_bound(I_DEFAULT)
{}


/*****************************************************************
 *
 *  Constructor of DiscreteMultivariateParametric class using
 *  the number of variables, the family of distribution,
 *  the (common) inf. bound, the sup. bounds, the set of parameters
 *  and a threshold on cumulated probabilities
 *
 **/

DiscreteMultivariateParametric::DiscreteMultivariateParametric(int inb_variable,
                                                               int iident,
                                                               unsigned int iinf_bound,
                                                               const std::vector<unsigned int>& isup_bound,
                                                               const std::vector<double>& vparameter,
                                                               const std::vector<double>& vprobability,
                                                               double cumul_threshold)
 : DiscreteMultivariateDistribution()
 , rparameter(NULL)
 , param_marginals()
{ init(inb_variable, iident, iinf_bound, isup_bound, vparameter, vprobability, cumul_threshold); }

/*****************************************************************
 *
 *  Copy constructor of DiscreteMultivariateParametric class
 *
 **/

DiscreteMultivariateParametric::DiscreteMultivariateParametric(const DiscreteMultivariateParametric& dist)
 : DiscreteMultivariateDistribution(dist)
 , rparameter(NULL)
 , sup_bound()
 , parameter()
 , probability()
 , param_marginals()
 , ident(I_DEFAULT)
 , inf_bound(I_DEFAULT)
{ copy(dist); }

/*****************************************************************
 *
 *  Copy a DiscreteMultivariateParametric and return a pointer
 *
 **/

DiscreteMultivariateParametric* DiscreteMultivariateParametric::copy_ptr() const
{
   DiscreteMultivariateParametric *res = new DiscreteMultivariateParametric(*this);

   return res;
}

/*****************************************************************
 *
 *  Destructor for DiscreteMultivariateParametric class
 *
 **/

DiscreteMultivariateParametric::~DiscreteMultivariateParametric()
{ remove(); }

/*****************************************************************
 *
 *  Set distribution type and parameters
 *  for DiscreteMultivariateParametric class
 *
 **/

void DiscreteMultivariateParametric::init(int inb_variable, int iident, unsigned int iinf_bound,
                                          const std::vector<unsigned int>& isup_bound,
                                          const std::vector<double>& vparameter,
                                          const std::vector<double>& vprobability,
                                          double cumul_threshold)
{
     unsigned int i,j,k;
     bool all_null;
     double s;
     std::vector<double> iprobability;

     remove();

     sup_bound = isup_bound;
     parameter = vparameter;
     probability = vprobability;
     ident = iident;
     inf_bound = iinf_bound;

     nb_variable = inb_variable;
     nb_value.resize(nb_variable);
     offset.resize(nb_variable);
     complement = 0;


     switch (ident) {
         case MNEGATIVE_MULTINOMIAL : {
             assert(parameter.size() == 1);
             assert(parameter[0] > 0.);
             assert(probability.size() == nb_variable);

             rparameter = new vector<double>(1);
             (*rparameter)[0] = 1.; // (1 - sum of (1-probabilities))
             for(i = 0; i < nb_variable; i++)
                 (*rparameter)[0] -= (1-probability[i]);
             assert(((*rparameter)[0] < 1.) && ((*rparameter)[0] > 0.));

             mean = parametric_mean_computation();
             variance = parametric_variance_computation();
             param_marginals.resize(nb_variable);
             sup_bound.resize(0);

             for(i = 0; i < nb_variable; i++){
                 offset[i] = inf_bound;
                 nb_value[i] = (int)round(inf_bound + (mean[i] - inf_bound + sqrt(variance[i])) * 20.);
                 if (nb_value[i] == inf_bound)
                     nb_value[i]++;

         if (probability[i] < 1.)
            param_marginals[i] = new DiscreteParametric(NEGATIVE_BINOMIAL,
                                                        inf_bound, I_DEFAULT,
                                                        parameter[0],
                                                        (*rparameter)[0] / ((*rparameter)[0] + 1. - probability[i]),
                                                        cumul_threshold);
            else
               param_marginals[i] = new DiscreteParametric(BINOMIAL,
                                                           inf_bound, inf_bound+1,
                                                           D_DEFAULT, 1.,
                                                           cumul_threshold);
             }
             break;
         }
         case MPOISSON : {
             assert(parameter.size() == nb_variable+1);
                         assert(parameter[0] >= 0);
             for(i = 0; i < nb_variable; i++){
                 assert(parameter[i+1] > 0.);
             }

             mean = parametric_mean_computation();
             variance = parametric_variance_computation();
             param_marginals.resize(nb_variable);
             sup_bound.resize(0);

             for(i = 0; i < nb_variable; i++){
                 offset[i] = inf_bound;
                 nb_value[i] = (int)round(inf_bound + (mean[i] - inf_bound + sqrt(variance[i])) * 20.);
                 if(nb_value[i] == inf_bound)
                     nb_value[i]++;
                 param_marginals[i] = new DiscreteParametric(POISSON,
                                                                                                         inf_bound, I_DEFAULT,
                                                                                                         parameter[0]+parameter[i+1],
                                                                                                         D_DEFAULT, cumul_threshold);
             }
             break;
         }
         case MMULTINOMIAL : {
             assert(parameter[0] >= inf_bound*nb_variable);
             s = 0;
             for(i = 0; i < nb_variable; i++){
                 assert(probability[i] > 0);
                 s += probability[i];
             }
             if(s != 1){
                 for(i = 0; i < nb_variable; i++){
                     probability[i] /= s;
                 }
                 s = 0;
                 for(i = 0; i < nb_variable-1; i++){
                     s += probability[i];
                 }
                 probability[nb_variable-1] = 1-s;
             }

             mean = parametric_mean_computation();
             variance = parametric_variance_computation();
             param_marginals.resize(nb_variable);
             sup_bound.resize(nb_variable);

             for(i = 0; i < nb_variable; i++){
                 offset[i] = inf_bound;
                 sup_bound[i] = parameter[0]-(nb_variable-1)*inf_bound;
                 nb_value[i] = sup_bound[i] + 1;
                 param_marginals[i] = new DiscreteParametric(BINOMIAL,
                                                                                                             inf_bound, sup_bound[i],
                                                                                                             D_DEFAULT,
                                                                                                             probability[i], cumul_threshold);
             }
             break;
         }
         default:
             assert(false);
             break;
     }
     total_nb_value = 1;
     all_null = (nb_value[0] == 0);
     for(i = 0; i < nb_variable; i++){
         total_nb_value *= nb_value[i];
         if (!(all_null) | (nb_value[0] > 0))
             all_null = false;
     }
     if (all_null)
         total_nb_value = 0;

     if (total_nb_value < MULTIVARIATE_DISTRIBUTION_MAX_VALUE)
         mass_computation();
         /*
          * MUNIFORM
          * */
}

/*****************************************************************
 *
 *  Compute the number of free parameters in DiscreteMultivariateParametric
 *
 **/

unsigned int DiscreteMultivariateParametric::nb_parameter_computation() const
{
     unsigned int bnb_parameter;
     switch (ident) {
         case MNEGATIVE_MULTINOMIAL :
             bnb_parameter = nb_variable + 1;
             break;
         case MPOISSON :
             bnb_parameter = nb_variable + 1;
             break;
         case MMULTINOMIAL :
             bnb_parameter = nb_variable;
         default :
             bnb_parameter = 0;
             break;
     }
     /*
      * MUNIFORM
      * */

     return bnb_parameter;
}

/*****************************************************************
 *
 *  Update the number of free parameters in DiscreteMultivariateParametric
 *
 **/

void DiscreteMultivariateParametric::nb_parameter_update()
{
   nb_parameter = nb_parameter_computation();
}

/*****************************************************************
 *
 * Get inferior bound of DiscreteMultivariateParametric distribution
 *
 **/

unsigned int DiscreteMultivariateParametric::get_inf_bound() const
{ return inf_bound; }

/*****************************************************************
 *
 * Get superior bounds of DiscreteMultivariateParametric distribution
 *
 **/

std::vector<unsigned int> DiscreteMultivariateParametric::get_sup_bound() const
{ return sup_bound; }

/*****************************************************************
 *
 * Get parameter of DiscreteMultivariateParametric distribution
 *
 **/

std::vector<double> DiscreteMultivariateParametric::get_parameter() const
{ return parameter; }


/*****************************************************************
 *
 * Get probabilities of DiscreteMultivariateParametric distribution
 *
 **/

std::vector<double> DiscreteMultivariateParametric::get_probability() const
{ return probability; }

/*****************************************************************
 *
 *  Get the probability value for a set of indices
 *
 **/

double DiscreteMultivariateParametric::get_mass(const index_array& indices) const
{ return get_mass(indices, false); }

/*****************************************************************
 *
 *  Get either the probability value for a set of indices,
 *  or the log probability
 *
 **/

double DiscreteMultivariateParametric::get_mass(const index_array& indices, bool log_computation) const
{
     double res = DiscreteMultivariateDistribution::get_mass(indices);

     if (res < 0){
         switch (ident) {
             case MNEGATIVE_MULTINOMIAL :
                 res = negative_multinomial_get_mass(indices, log_computation);
                 break;
             case MPOISSON :
                 res = mpoisson_get_mass(indices, log_computation);
                 break;
             case MMULTINOMIAL :
                 res = mmultinomial_get_mass(indices, log_computation);
                 break;
             default :
                 break;
         }

         if (log_computation)
             set_mass(indices, exp(res));
         else
             set_mass(indices, res);

             /*
              * MUNIFORM
              * */
     }
         else {
             if(log_computation)
                 res = log(res);
         }
   return res;
}

/*****************************************************************
 *
 *  Analytic computation of expectations in DiscreteMultivariateParametric
 *
 **/

std::vector<double> DiscreteMultivariateParametric::parametric_mean_computation() const
{
     std::vector<double> res;
     unsigned int i;

     res.resize(nb_variable);

     switch (ident) {
         case MNEGATIVE_MULTINOMIAL :
             for(i = 0; i < nb_variable; i++)
                 res[i] = inf_bound + parameter[0] * (1-probability[i]) / (*rparameter)[0];
             break;
         case MPOISSON :
             for(i = 0; i < nb_variable; i++)
                 res[i] = inf_bound + parameter[i+1]+parameter[0];
             break;
         case MMULTINOMIAL :
             for(i = 0; i < nb_variable; i++)
                 res[i] = (parameter[0]-(nb_variable-1)*inf_bound)*probability[i];
             break;
         default :
             res = mean;
             break;
     }
     /*
      * MUNIFORM
      * */
     return res;
}

/*****************************************************************
 *
 *  Analytic computation of variances in DiscreteMultivariateParametric
 *
 **/

std::vector<double> DiscreteMultivariateParametric::parametric_variance_computation() const
{
     std::vector<double> res;
     unsigned int i;

     res.resize(nb_variable);

     switch (ident) {
         case MNEGATIVE_MULTINOMIAL :
             for(i = 0; i < nb_variable; i++)
                 res[i] = parameter[0] * (1. - probability[i]) * (1. - probability[i] +(*rparameter)[0]) / ((*rparameter)[0] * (*rparameter)[0]);
             break;
         case MPOISSON :
             for(i = 0; i < nb_variable; i++)
                 res[i] = parameter[i+1] + parameter[0];
             break;
         case MMULTINOMIAL :
             for(i = 0; i < nb_variable; i++)
                 res[i] = (parameter[0]-nb_variable*inf_bound) * probability[i] * (1 - probability[i]);
             break;
         default :
             res = variance;
             break;
     }

     /*
      * MUNIFORM
      * */

     return res;
}

/*****************************************************************
 *
 *  Analytic computation of covariance between two given variables
 *  in DiscreteMultivariateParametric
 *
 **/

double DiscreteMultivariateParametric::parametric_covariance_computation(unsigned int ivariable1,
                                                                         unsigned int ivariable2) const
{
     double res = D_INF;
     std::vector<double> variance;
     assert((ivariable1 < nb_variable) && (ivariable2 < nb_variable));
     if (ivariable1 == ivariable2){
         variance = parametric_variance_computation();
         res = variance[ivariable1];
     } else {
         switch (ident) {
             case MNEGATIVE_MULTINOMIAL :
                 res = parameter[0] * (1. - probability[ivariable1]) * (1. - probability[ivariable2]) / ((*rparameter)[0] * (*rparameter)[0]);
                 break;
             case MPOISSON :
                 res = parameter[0];
                 break;
             case MMULTINOMIAL :
                 res = -(parameter[0]-nb_variable*inf_bound) * probability[ivariable1] * probability[ivariable2];
                 break;
             default :
                 break;
         }
         /*
          * MUNIFORM
          * */
     }

     return res;
}

/*****************************************************************
 *
 *  Conditional distribution of variable i given variables 0,...,i-1
 *
 **/

DiscreteParametric* DiscreteMultivariateParametric::conditional_distribution(unsigned int ivariablei,
                                                                             const index_array& values) const
{
     DiscreteParametric *res = NULL;

     switch (ident) {
         case MIID : {
             res = extract_marginal(ivariablei);
             break;
         }
         case MNEGATIVE_MULTINOMIAL : {
             res = negative_multinomial_conditional(ivariablei, values);
             break;
         }
         case MMULTINOMIAL : {
             res = mmultinomial_conditional(ivariablei, values);
             break;
         }
         default :
             break;
     }

     return res;
}

/*****************************************************************
 *
 *  Compute marginal distributions in DiscreteMultivariateParametric
 *
 **/

void DiscreteMultivariateParametric::compute_param_marginals(double cumul_threshold)
{
     unsigned int i,j,k;
     double p;
     std::vector<double> iprobability;

     param_marginals.resize(nb_variable);

     switch (ident) {
         case MNEGATIVE_MULTINOMIAL :
             for(i = 0; i < nb_variable; i++)
                 param_marginals[i] = new DiscreteParametric(NEGATIVE_BINOMIAL,
                                                                                                         inf_bound, I_DEFAULT,
                                                                                                         parameter[0],
                                                                                                         (*rparameter)[0] / ((*rparameter)[0] + 1. - probability[i]),
                                                                                                         cumul_threshold);
             break;
         case MPOISSON :
             for(i = 0; i < nb_variable; i++)
                 param_marginals[i] = new DiscreteParametric(POISSON, inf_bound,
                                                                                                         I_DEFAULT,
                                                                                                         parameter[0] + parameter[i+1],
                                                                                                         D_DEFAULT, cumul_threshold);
             break;
         case MMULTINOMIAL :
             for(i = 0; i < nb_variable; i++)
                 param_marginals[i] = new DiscreteParametric(BINOMIAL, inf_bound, sup_bound[i], D_DEFAULT, probability[i], cumul_threshold);
             break;
         default :
             break;
     }

     /*
      * MUNIFORM,
      * */

}

/*****************************************************************
 *
 *  Extend the support of a marginal distribution
 *  in DiscreteMultivariateParametric for a given variable
 *  using the new bound
 *
 **/

void DiscreteMultivariateParametric::extend_param_marginal(unsigned int ivariable, int value_update) const
{
   assert(ivariable < nb_variable);

   assert(param_marginals.size() == nb_variable);

   if (param_marginals[ivariable]->nb_value < value_update)
      param_marginals[ivariable]->computation(value_update, 1.);

}

/*****************************************************************
 *
 *  Return marginal distribution of a given variable
 *  in DiscreteMultivariateParametric
 *
 **/

Distribution* DiscreteMultivariateParametric::get_marginal(unsigned int ivariable) const
{ return extract_marginal(ivariable); }

/*****************************************************************
 *
 *  Return parametric marginal distribution of a given variable
 *  in DiscreteMultivariateParametric
 *  (a new instance is allocated)
 *
 **/

DiscreteParametric* DiscreteMultivariateParametric::extract_marginal(unsigned int ivariable) const
{

   DiscreteParametric* res = NULL;

   assert(ivariable < nb_variable);

   assert(param_marginals.size() == nb_variable);

   res = new DiscreteParametric(*param_marginals[ivariable]);
   return res;
}

/*****************************************************************
 *
 *  Print a DiscreteMultivariateParametric on a single line
 *  using an output stream
 *
 **/

ostream& DiscreteMultivariateParametric::line_write(ostream& os) const
{

   os << STAT_TREES_multivariate_distribution_word[this->ident];

   os << " ;  " << nb_variable << " " << STAT_word[STATW_VARIABLES];

   return os;
}

/*****************************************************************
 *
 *  Print a DiscreteMultivariateParametric into a file
 *  using a StatError object, the path
 *  and a flag on the level of detail
 *
 **/

bool DiscreteMultivariateParametric::ascii_write(StatError& error, const char * path,
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
      ascii_print(out_file, true, exhaustive, false, NULL);
   }

   return status;
}

/*****************************************************************
 *
 *  Print a DiscreteMultivariateParametric in a spreadsheet fashion
 *  using a StatError object and the path
 *
 **/

bool DiscreteMultivariateParametric::spreadsheet_write(StatError& error,
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
      spreadsheet_write(out_file);
   }

   return status;
}

/*****************************************************************
 *
 *  Gnuplot output for DiscreteMultivariateParametric class
 *  using a StatError object, a prefix for the files
 *  and the title of output figures
 *
 **/

bool DiscreteMultivariateParametric::plot_write(StatError& error,
                                                const char * prefix,
                                                const char * title) const
{
   bool status = plot_write(prefix, title);

   error.init();

   if (!status)
     error.update(STAT_error[STATR_FILE_PREFIX]);

   return status;
}

/*****************************************************************
 *
 *  Matplotlib output for DiscreteMultivariateParametric class
 *  using potential data
 *
 **/

MultiPlotSet* DiscreteMultivariateParametric::get_plotable(const Statiskit::Marginal::Multivariate::CountTable *histo) const
{
   // plot_set[variable][view]
   register int i, j;
   const int nb_plot_set = param_marginals.size();
   ostringstream title;
   DiscreteParametricModel *pmodel = NULL;
   FrequencyDistribution *chisto = NULL;
   MultiPlotSet *plot_set = new MultiPlotSet(nb_plot_set), *plot_distr = NULL;
   MultiPlot *mplot = NULL;
   Statiskit::Marginal::Multivariate::CountTable::marginal_type::ptr_type mhisto;

   plot_set->border = "15 lw 0";

   title.str("");
   title << STAT_MULTIVARIATE_label[TREESTATL_DISCRETE_MULTIVARIATE_DISTRIBUTION];
   if (histo != NULL)
      title << " " << STAT_label[STATL_FIT];

   plot_set->title = title.str();

   for(i = 0; i < nb_plot_set; i++)
   {
      // it is required in DiscreteParametricModel->plotable
      // that this->frequencydistribution != NULL
      if (histo != NULL)
      {
         mhisto = Statiskit::marginalizing(*histo, i); // shared_ptr
         chisto = Scalar_stat_histogram_to_frequency_distribution(*mhisto);

         if (param_marginals[i]->nb_value < chisto->nb_value)
            extend_param_marginal(i, chisto->nb_value);

         pmodel = new DiscreteParametricModel(*param_marginals[i], chisto);
         plot_distr = pmodel->get_plotable();
         delete pmodel;
         pmodel = NULL;
         (*plot_set)[i].resize((*plot_distr)[0].size());
         for(j = 0; j < (*plot_distr)[0].size(); j++)
            (*plot_set)[i][j] = (*plot_distr)[0][j];
      }
      else
      {
         plot_distr = param_marginals[i]->get_plotable();
      // mplot = new MultiPlot((*plot_distr)[0]);
         (*plot_set)[i] = (*plot_distr)[0];
      }
      // delete mplot;
      // mplot = NULL;
      title.str("");
      title << STAT_label[STATL_VARIABLE] << " " << i+1;
      (*plot_set)[i].title = title.str();
      if (chisto != NULL)
      {
         delete chisto;
         chisto = NULL;
      }
   } // end for i

   return plot_set;
}

/*****************************************************************
 *
 *  Print a DiscreteMultivariateParametric into an output stream, using
 *  a flag on the level of detail
 *
 **/

ostream& DiscreteMultivariateParametric::ascii_write(ostream& os, bool exhaustive) const
{ return ascii_write(os, exhaustive, false); }

/*****************************************************************
 *
 *  Print a DiscreteMultivariateParametric into an output stream, using
 *  a flag on the level of detail and a flag on writing into a file or not
 *
 **/

ostream& DiscreteMultivariateParametric::ascii_write(ostream& os, bool exhaustive,
                                                     bool file_flag) const
{ return ascii_print(os, file_flag, exhaustive, false, NULL); }

/*****************************************************************
 *
 *  Print a DiscreteMultivariateParametric in a spreadsheet fashion
 *  into an output stream
 *
 **/

ostream& DiscreteMultivariateParametric::spreadsheet_write(ostream& os) const
{
   spreadsheet_print(os);
   if (complement > 0.)
   {
      os << STAT_label[STATL_UNPROPER] << " " << STAT_label[STATL_DISTRIBUTION] << "\t"
         << STAT_label[STATL_COMPLEMENTARY_PROBABILITY] << "\t" << complement << "\t" << endl;
    }

   spreadsheet_characteristic_print(os, true);


   os << STAT_label[STATL_INFORMATION] << "\t" << information_computation() << endl;

   DiscreteMultivariateDistribution::spreadsheet_print(os);
}

/*****************************************************************
 *
 *  Gnuplot output for DiscreteMultivariateParametric class
 *  using a prefix for the files and the title of output figures
 *
 **/

bool DiscreteMultivariateParametric::plot_write(const char * prefix, const char * title) const
{
   bool status;
   unsigned int i;
   int variable= 0;
   DiscreteParametricProcess *pprocess = NULL;
   DiscreteParametric **pobservation = NULL;


   if (param_marginals.empty())
      status = false;
   else
   {
      pobservation = new DiscreteParametric*[nb_variable];
      for(i = 0; i < nb_variable; i++)
         if (param_marginals[i] != NULL)
            pobservation[i] = new DiscreteParametric(*param_marginals[i]);
         else
            pobservation[i] = NULL;
       pprocess = new DiscreteParametricProcess(nb_variable, pobservation);
       status = pprocess->plot_print(prefix, title, 0);
       for(i = 0; i < nb_variable; i++)
          if (pobservation[i] != NULL)
             delete pobservation[i];
       delete [] pobservation;
       delete pprocess;
   }
   return status;
}


/*****************************************************************
 *
 *  Print DiscreteMultivariateParametric parameters using an output stream
 *
 **/

ostream& DiscreteMultivariateParametric::ascii_print(ostream &os, bool comment_flag,
                                                     bool exhaustive, bool file_flag, bool cumul_flag,
                                                     const Statiskit::Marginal::Multivariate::CountTable *histo) const
{
   unsigned int i;
   header_ascii_print(os);

   os << STAT_TREES_multivariate_distribution_word[ident] << "   ";

   if (inf_bound != I_DEFAULT)
   {
      os << STAT_word[STATW_INF_BOUND] << " : " << inf_bound << "   ";
      os << " ; ";
   }
   if (!sup_bound.empty())
   {
      os << STAT_word[STATW_SUP_BOUND] << " : ";
      for(i = 0; i < sup_bound.size(); i++)
         os << sup_bound[i] << "   ";
      os << " ; ";
   }

   if (!parameter.empty())
      os << STAT_word[STATW_PARAMETER] << " : " << parameter << "   ";

   if (!probability.empty())
      os << STAT_word[STATW_PROBABILITY] << " : " << probability;

   os << endl;

   if (exhaustive)
   {
      if (complement > 0.)
      {
         os << STAT_label[STATL_UNPROPER] << " " << STAT_label[STATL_DISTRIBUTION] << " ("
            << STAT_label[STATL_COMPLEMENTARY_PROBABILITY] << ": " << complement << ")" << endl;
       }

      // write marginals, means, etc.
      ascii_characteristic_print(os, true, comment_flag);


      if (file_flag)
         os << "# ";

      if ((total_nb_value > 0) && (total_nb_value < MULTIVARIATE_DISTRIBUTION_MAX_VALUE))
      {
         mass_computation();
         os << STAT_label[STATL_INFORMATION] << ": " << information_computation() << endl;
      }

      if ((total_nb_value > 0) && (total_nb_value < MULTIVARIATE_DISTRIBUTION_MAX_VALUE))
         // write probabilities
         DiscreteMultivariateDistribution::probability_ascii_print(os, comment_flag, cumul_flag);
   }

   return os;
}

/*****************************************************************
 *
 *  Print DiscreteMultivariateParametric parameters
 *  in a spreadsheet fashion
 *
 **/

ostream& DiscreteMultivariateParametric::spreadsheet_print(ostream& os) const
{
   os << STAT_discrete_distribution_word[ident] << "\t";

   if (inf_bound != I_DEFAULT)
      os << STAT_word[STATW_INF_BOUND] << "\t" << inf_bound << "\t";
   if (!sup_bound.empty())
      os << STAT_word[STATW_SUP_BOUND] << "\t" << sup_bound << "\t";

   if (!parameter.empty())
      os << STAT_word[STATW_PARAMETER] << "\t" << parameter << "\t";

   if (probability.empty())
      os << STAT_word[STATW_PROBABILITY] << "\t" << probability;

   os << endl;

   return os;
}

/*****************************************************************
 *
 *  Print the marginal characteristics (mean, variance,...)
 *
 **/

ostream& DiscreteMultivariateParametric::ascii_characteristic_print(ostream &os,
                                                                    bool shape,
                                                                    bool comment_flag) const
{
   unsigned int i;

   if (!param_marginals.empty())
   {
      for(i = 0; i < nb_variable; i++)
         if (param_marginals[i] != NULL)
         {
            if (comment_flag)
               os << "# " ;
            os << STAT_word[STATW_VARIABLE] << ": " << i << endl;
            if (comment_flag)
               os << "# " ;
            param_marginals[i]->ascii_print(os);
            param_marginals[i]->ascii_characteristic_print(os, shape, comment_flag);
            os << endl;
         }
   }
   return os;
}

/*****************************************************************
 *
 *  Print the marginal characteristics (mean, variance,...)
 *  in a spreadsheet fashion
 *
 **/

ostream& DiscreteMultivariateParametric::spreadsheet_characteristic_print(ostream &os,
                                                                          bool shape) const
{
   unsigned int i;

   if (!param_marginals.empty())
   {
      for(i = 0; i < nb_variable; i++)
         if (param_marginals[i] != NULL)
         {
            os << STAT_label[STATL_VARIABLE] << "\t" << i << endl;
            param_marginals[i]->spreadsheet_characteristic_print(os, shape);
            os << endl;
         }
   }
   return os;
}

/*****************************************************************
 *
 *  Variable permutation for DiscreteMultivariateParametric class
 *  (checking or not validity of permutation)
 *
 **/

void DiscreteMultivariateParametric::variable_permutation(StatError& error,
                                                          int* perm, bool check_flag)
{
   unsigned int i;
   bool status = true;
   // each element of check_perm must be used exactly once

   DiscreteMultivariateDistribution::variable_permutation(error, perm, check_flag);
   if (error.get_nb_error() > 0)
      status = false;

   if (status)
   {
      std::vector< DiscreteParametric * >  pparam_marginals(nb_variable);
      for(i = 0; i < nb_variable; i++)
         pparam_marginals[perm[i]] = param_marginals[i];

      for(i = 0; i < nb_variable; i++)
         param_marginals[i] = pparam_marginals[i];

      switch (ident)
      {
         case MNEGATIVE_MULTINOMIAL :
         {
            std::vector<double> pprobability(nb_variable);

            for(i = 0; i < nb_variable; i++)
               pprobability[perm[i]] = probability[i];

            for(i = 0; i < nb_variable; i++)
               probability[i] = probability[i];
            break;
         }
      }
   }
}

/*****************************************************************
 *
 *  Destructor for DiscreteMultivariateParametric class
 *
 **/

void DiscreteMultivariateParametric::remove()
{
   unsigned int i;

   if (!param_marginals.empty())
   {
      for(i = 0; i < param_marginals.size(); i++)
      {
         if (param_marginals[i] != NULL)
         {
            delete param_marginals[i];
            param_marginals[i] = NULL;
         }
      }
   }
   if (rparameter != NULL)
   {
      delete rparameter;
      rparameter = NULL;
   }
}

/*****************************************************************
 *
 *  Copy operator for DiscreteMultivariateParametric class
 *
 **/

void DiscreteMultivariateParametric::copy(const DiscreteMultivariateParametric& dist)
{
   unsigned int i;

   param_marginals.resize(dist.param_marginals.size());

   for(i = 0; i < param_marginals.size(); i++)
      param_marginals[i] = new DiscreteParametric(*dist.param_marginals[i]);

   ident = dist.ident;
   inf_bound = dist.inf_bound;
   sup_bound = dist.sup_bound;
   parameter = dist.parameter;
   probability = dist.probability;

   if ((rparameter != NULL) && (dist.rparameter != NULL))
      *rparameter = *dist.rparameter;
   if ((rparameter == NULL) && (dist.rparameter != NULL))
      rparameter = new std::vector<double>(*dist.rparameter);

}

/*****************************************************************
 *
 *  Parse a DiscreteMultivariateParametric from a file
 *  Parsing begins at INF_BOUND: the number of variables
 *  and distribution types must have been parsed before
 *
 **/

DiscreteMultivariateParametric*
Stat_trees::discrete_multivariate_parametric_parsing(StatError &error,
                                                     std::ifstream &in_file,
                                                     int &line,
                                                     RWCTokenizer &next,
                                                     int iident,
                                                     unsigned int inb_variable,
                                                     double cumul_threshold,
                                                     int min_inf_bound)

{
   RWLocaleSnapshot locale("en");
   RWCString buffer, token;
   size_t position;
   bool status = true, lstatus;
   register int i, var, correction_line = 0;
   long iinf_bound;
   double rparameter = 0; // MNEGATIVE_MULTINOMIAL
   std::vector<unsigned int> isup_bound;
   std::vector<long> lsup_bound;
   std::vector<double> vparameter;
   std::vector<double> vprobability;
   DiscreteMultivariateParametric *dist = NULL;

   if ((iident < MPOISSON) || (iident > MMULTINOMIAL))
   {
      status = false;
      error.update(STAT_parsing[STATP_DISTRIBUTION_NAME], line, 0);
      return dist;
   }

   i = 0;

   if (&next == NULL)
      next = RWCTokenizer(buffer);

   while ((!((token = next()).isNull())) && (i < 4))
   {
      // read inferior bound

      switch (i)
      {
         case 0 :
         {
            // inferior bound keyword
            if (token != STAT_word[STATW_INF_BOUND])
            {
               status = false;
               error.correction_update(STAT_parsing[STATP_PARAMETER_NAME], STAT_word[STATW_INF_BOUND], line , i + 2);
            }
            break;
          }
         case 1 :
         {
            // test separator ":"
            if (token != ":")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 2);
            }
            break;
         }
         // value of inf. bound
         case 2 :
         {
            lstatus = locale.stringToNum(token , &iinf_bound);
            if ((lstatus) && ((iinf_bound < min_inf_bound) || (iinf_bound > MAX_INF_BOUND)))
               lstatus = false;
            break;

         }
         case 3 :
         {
            // test separator ";"
            if (token != ";")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 2);
            }
            break;
         }
      }
      i++;
   } // end while

   // i == 4 at this point, in principle
   if (i != 4)
      status = false;

   if (!status)
      error.update(STAT_parsing[STATP_FORMAT] , line);

   i = 2;
   // test for SUP_BOUND
   if (status & (iident == MMULTINOMIAL))
   {
      if (token != STAT_word[STATW_SUP_BOUND])
      {
         status = false;
         error.correction_update(STAT_parsing[STATP_PARAMETER_NAME], STAT_word[STATW_SUP_BOUND], line, i + 1);
      }
      else
      {
         isup_bound.resize(inb_variable);
         lsup_bound.resize(inb_variable);
         while ((!((token = next()).isNull())) && (i < 4))
         {
            switch ((i - 1) % 4) // blocks "SUP_BOUND : val1  valK ;"
            {
               // test separator ":"
               case 1:
               {
                  if (token != ":")
                  {
                     status = false;
                     error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5);
                  }
                  break;
               }

               case 2:
               {
                  var = 0;
                  lstatus = lstatus && locale.stringToNum(token, &lsup_bound[var]);
                  isup_bound[var] = (unsigned int)(lsup_bound[var]);
                  rparameter = isup_bound[var];
                  var++;
                  while ((!((token = next()).isNull())) && (var < inb_variable))
                  {
                     lstatus = lstatus && locale.stringToNum(token , &lsup_bound[var]);
                     isup_bound[var] = (unsigned int)(lsup_bound[var]);
                     if (rparameter != isup_bound[var])
                     {
                        lstatus = false;
                        error.update(STAT_parsing[STATP_PARAMETER_VALUE], line, i + 5 + var);
                        ostringstream correction_message;
                        correction_message << "!= " << token;
                        error.correction_update(STAT_parsing[STATP_PARAMETER_VALUE], (correction_message.str()).c_str(), line);
                     }
                     var++;
                  }
                  if (token != ";")
                  {
                     status = false;
                     error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5);
                  }
                  break;
               }

               // test separator ";"
               default:
               {
                  status = false;
                  error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5);
                  break;
               }
            } // end switch ((i - 1) % 4)

            if (!lstatus)
            {
               status = false;
               error.update(STAT_parsing[STATP_PARAMETER_VALUE], line, i + 5);
            }
            i++;
         } // end while
         if ((i != 4))
         {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
         }
         i = 2;
      }
   }

   // test for parameter keyword
   if ((status) && ((((iident == MPOISSON) || (iident == MNEGATIVE_MULTINOMIAL)
        || (iident == MMULTINOMIAL)) &&(token != STAT_word[STATW_PARAMETER]))))
   {
      status = false;
      error.correction_update(STAT_parsing[STATP_PARAMETER_NAME], STAT_word[STATW_PARAMETER], line, i + 1);
   }

   if (iident == MPOISSON)
      vparameter.resize(inb_variable + 1);

   if ((iident == MNEGATIVE_MULTINOMIAL) || (iident == MMULTINOMIAL))
   {
      vparameter.resize(1);
      vprobability.resize(inb_variable);
   }

   // test for parameter name
   while (status && (!((token = next()).isNull())))
   {
      switch ((i - 1) % 4) // blocks "param : val1  valK ;"
      {
         case 0 :
         {
            switch ((i - 1) / 4)
            {
               // probability keyword
               case 1 :
               {
                  if (((iident == MNEGATIVE_MULTINOMIAL) || (iident == MMULTINOMIAL))
                      && (token != STAT_word[STATW_PROBABILITY]))
                  {
                     status = false;
                     error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PROBABILITY] , line , i + 5);
                  }
                  break;
               }

               default :
               {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PROBABILITY] , line , i + 5);
               }
            } // end switch ((i - 1) / 4)
            break;
         }

         // test separator ":"
         case 1 :
         {
            if (token != ":")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5);
            }
            break;
         }

         // parameter value
         case 2 :
         {
            switch ((i - 1) / 4)
            {
               // parameter value
               case 0 :
               {
                  if ((iident == MNEGATIVE_MULTINOMIAL) || (iident == MMULTINOMIAL)
                      || (iident == MPOISSON))
                  {
                    lstatus = locale.stringToNum(token , &vparameter[0]);
                    if ((lstatus) && (vparameter[0] <= 0.))
                       lstatus = false;
                    if ((iident == MMULTINOMIAL) &&
                        (rparameter != ((int)vparameter[0]) - (inb_variable-1)*iinf_bound))
                    {
                       lstatus = false;
                       error.update(STAT_parsing[STATP_PARAMETER_VALUE], line, i + 5);
                       ostringstream correction_message;
                       correction_message << " == " << rparameter + (inb_variable-1)*iinf_bound;
                       error.correction_update(STAT_parsing[STATP_PARAMETER_VALUE], (correction_message.str()).c_str(), line);
                    }
                    if (iident == MPOISSON)
                    {
                       var = 1;
                       while ((!((token = next()).isNull())) && (var < inb_variable+1))
                       {
                          lstatus = lstatus && locale.stringToNum(token , &vparameter[var]);
                          if (vparameter[var] <= 0.)
                          {
                             lstatus = false;
                             correction_line = var - 1;
                             if (token != ";")
                             {
                                error.update(STAT_parsing[STATP_PARAMETER_VALUE], line, i + 5 + correction_line + 1);
                                ostringstream correction_message;
                                correction_message << "> 0";
                                error.correction_update(STAT_parsing[STATP_PARAMETER_VALUE], (correction_message.str()).c_str(), line);
                             }
                             else
                             {
                                error.update(STAT_parsing[STATP_NB_PARAMETER], line, i + 5 + correction_line);
                                ostringstream correction_message;
                                correction_message << "== " << inb_variable + 1;
                                error.correction_update(STAT_parsing[STATP_NB_PARAMETER], (correction_message.str()).c_str(), line);
                             }
                          }
                          var++;
                       }
                       if (token != ";")
                       {
                          status = false;
                          error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5 + inb_variable + 1);
                       }
                       i++;
                       // nothing more to read
                       while (status && (!((token = next()).isNull())))
                          i++;
                    }
                  }
                  break;
               }

               // probability value
               case 1 :
               {
                  if ((iident == MNEGATIVE_MULTINOMIAL) || (iident == MMULTINOMIAL))
                  {
                     var = 0;
                     rparameter = 0.;
                     lstatus = lstatus && locale.stringToNum(token , &vprobability[var]);
                     if (!((vprobability[var] > 0.) && (vprobability[var] < 1.)))
                     {
                        lstatus = false;
                        error.update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], line, i + 5 + var);
                        ostringstream correction_message;
                        correction_message << " != " << token;
                        error.correction_update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], (correction_message.str()).c_str(), line);
                     }
                     if (iident == MNEGATIVE_MULTINOMIAL)
                        rparameter += 1 - vprobability[var];
                     if (iident == MMULTINOMIAL)
                        rparameter += vprobability[var];
                     var++;
                     while ((!((token = next()).isNull())) && (var < inb_variable))
                     {
                        lstatus = lstatus && locale.stringToNum(token , &vprobability[var]);
                        if (!((vprobability[var] > 0.) && (vprobability[var] < 1.)))
                        {
                           lstatus = false;
                           error.update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], line, i + 5 + var);
                           ostringstream correction_message;
                           correction_message << "!= " << token;
                           error.correction_update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], (correction_message.str()).c_str(), line);
                        }
                        if (iident == MNEGATIVE_MULTINOMIAL)
                           rparameter += 1 - vprobability[var];
                        if (iident == MMULTINOMIAL)
                           rparameter += vprobability[var];
                        var++;
                     }

                     if (lstatus)
                     {
                        switch (iident)
                        {
                           case MNEGATIVE_MULTINOMIAL:
                           {
                              if ((rparameter <= 0.) || (rparameter >= 1.))
                              {
                                 lstatus = false;
                                 ostringstream correction_message;
                                 correction_message << "in (0, 1)";
                                 error.correction_update(STAT_TREES_parsing[TREESTATP_BAD_ONE_MINUS_PROBABILITIES],
                                                         (correction_message.str()).c_str(), line);
                              }
                              break;
                           }
                           case MMULTINOMIAL:
                           {
                              if (abs(1-rparameter) > 1.-cumul_threshold)
                              {
                                 lstatus = false;
                                 ostringstream correction_message;
                                 correction_message << "in (0, 1)";
                                 error.correction_update(STAT_TREES_parsing[TREESTATP_BAD_ONE_MINUS_PROBABILITIES],
                                                         (correction_message.str()).c_str(), line);
                              }
                              break;
                           }
                        }
                     }
                  }
                  if (token != ";")
                  {
                     status = false;
                     error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5);
                  }
                  i++;
                  break;
               }
            } // end switch (i-1)/4
            break;
         } // end case 2

         // test separator ";"
         case 3 :
         {
            if (token != ";")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5 + correction_line);
            }
            break;
         }
      } // end switch ((i - 1) % 4)

      if (!lstatus)
      {
         status = false;
         error.update(STAT_parsing[STATP_PARAMETER_VALUE], line, i + 5 + correction_line);
      }
      i++;
   } // end while

   if (i > 0)
   {
      if ((iident == MNEGATIVE_MULTINOMIAL) && (i != 9))
      {
         status = false;
         error.update(STAT_parsing[STATP_FORMAT] , line);
      }
      if ((iident == MMULTINOMIAL) && (i != 9))
      {
         // ; PARAMETER : ... ; PROBABILITY : ... ;
         status = false;
         error.update(STAT_parsing[STATP_FORMAT] , line);
      }
      if ((iident == MPOISSON) && (i != 5))
      {
         // ; PARAMETER : ... ;
         status = false;
         error.update(STAT_parsing[STATP_FORMAT] , line, i + 4 + inb_variable);
      }
   }

   if (status)
      dist = new DiscreteMultivariateParametric(inb_variable, iident,
                                                iinf_bound, isup_bound,
                                                vparameter, vprobability,
                                                cumul_threshold);
   return dist;
}

/*****************************************************************
 *
 *  Default constructor of IidDiscreteMultivariateParametric class
 *
 **/

IidDiscreteMultivariateParametric::IidDiscreteMultivariateParametric(int iident, int inb_variable,
                                                                     unsigned int iinf_bound,
                                                                     unsigned int isup_bound,
                                                                     double dparameter,
                                                                     double dprobability,
                                                                     double cumul_threshold)
 : DiscreteMultivariateParametric()
 , uident(iident)
 , marginal(NULL)
 , sup_bound1(isup_bound)
 , parameter1(dparameter)
 , probability1(dprobability)
{
   double cmean, cvariance;
   unsigned int i, coffset;

   marginal = new DiscreteParametric(iident, iinf_bound, isup_bound,
                                     dparameter, dprobability, cumul_threshold);

   ident = MIID;

   coffset = marginal->offset;
   marginal->mean_computation();
   cmean = marginal->mean;
   marginal->variance_computation();
   cvariance = marginal->variance;
   complement = marginal->complement;

   nb_variable = inb_variable;
   nb_value.resize(nb_variable);
   offset.resize(nb_variable);
   total_nb_value = pow(marginal->nb_value, nb_variable);
   mean.resize(nb_variable);
   variance.resize(nb_variable);
   marginals.resize(nb_variable);
   sup_bound.resize(nb_variable);

   if (marginal->parameter == D_DEFAULT)
      parameter.resize(0);
   else
   {
      parameter.resize(1);
      parameter[0] = dparameter;
   }

   if (marginal->probability == D_DEFAULT)
      probability.resize(0);
   else
   {
      probability.resize(1);
      probability[0] = dprobability;
   }

   for(i = 0; i < nb_variable; i++)
   {
      mean[i] = cmean;
      variance[i] = cvariance;
      offset[i] = coffset;
      nb_value[i] = marginal->nb_value;
      marginals[i] = new DiscreteParametric(*marginal);
      sup_bound[i] = isup_bound;
   }

   mass->resize(0);

   val_strides_computation();

}

/*****************************************************************
 *
 *  Copy constructor of IidDiscreteMultivariateParametric class
 *
 **/

IidDiscreteMultivariateParametric::IidDiscreteMultivariateParametric(const IidDiscreteMultivariateParametric& dist)
 : DiscreteMultivariateParametric(dist)
 , uident(MIID)
 , marginal(NULL)
 , sup_bound1(I_DEFAULT)
 , parameter1(D_DEFAULT)
 , probability1(D_DEFAULT)
{ copy(dist); }

/*****************************************************************
 *
 *  Copy an IidDiscreteMultivariateParametric and return a pointer
 *
 **/

IidDiscreteMultivariateParametric* IidDiscreteMultivariateParametric::copy_ptr() const
{
   IidDiscreteMultivariateParametric *res = new IidDiscreteMultivariateParametric(*this);

   return res;
}

/*****************************************************************
 *
 *  Destructor for IidDiscreteMultivariateParametric class
 *
 **/

IidDiscreteMultivariateParametric::~IidDiscreteMultivariateParametric()
{ remove(); }

/*****************************************************************
 *
 *  Compute the number of free parameters in DiscreteMultivariateParametric
 *
 **/

unsigned int IidDiscreteMultivariateParametric::nb_parameter_computation() const
{
   unsigned int bnbparameter;

   assert(marginal != NULL);

   return marginal->nb_parameter_computation();
}

/*****************************************************************
 *
 *  Print an IidDiscreteMultivariateParametric on a single line
 *  using an output stream
 *
 **/

ostream& IidDiscreteMultivariateParametric::line_write(ostream& os) const
{

   os << STAT_TREES_multivariate_distribution_word[this->ident];

   os << " : " << STAT_discrete_distribution_word[this->uident];

   os << " ;  " << nb_variable << " " << STAT_word[STATW_VARIABLES];

   return os;
}

/*****************************************************************
 *
 *  Print IidDiscreteMultivariateParametric parameters using an output stream
 *
 **/

ostream& IidDiscreteMultivariateParametric::ascii_print(ostream &os, bool comment_flag,
                                                        bool exhaustive, bool cumul_flag,
                                                        bool file_flag,
                                                        const Statiskit::Marginal::Multivariate::CountTable *histo) const
{
   unsigned int i;

   header_ascii_print(os);

   os << STAT_TREES_multivariate_distribution_word[this->ident] << endl;

   // os << nb_variable << " " << STAT_word[STATW_VARIABLES] << endl;

   if (marginal != NULL)
   {
      marginal->ascii_print(os);
      marginal->ascii_characteristic_print(os, exhaustive, comment_flag);
   }

   return os;
}

/*****************************************************************
 *
 *  Destructor for IidDiscreteMultivariateParametric class
 *
 **/

void IidDiscreteMultivariateParametric::remove()
{
   if (marginal != NULL)
   {
      delete marginal;
      marginal = NULL;
   }
}

/*****************************************************************
 *
 *  Copy operator for DiscreteMultivariateDistribution class
 *
 **/

void IidDiscreteMultivariateParametric::copy(const IidDiscreteMultivariateParametric& dist)
{
   // DiscreteMultivariateParametric::copy(markov) must be used before
   // (or DiscreteMultivariateParametric::DiscreteMultivariateParametric)
   // as well as remove()
   marginal = new DiscreteParametric(*dist.marginal);
   uident = MIID;
   sup_bound1 = dist.sup_bound1;
   parameter1 = dist.parameter1;
   probability1 = dist.probability1;
}

/*****************************************************************
 *
 *  Parse an IidDiscreteMultivariateParametric from a file
 *
 **/

IidDiscreteMultivariateParametric* Stat_trees::iid_discrete_multivariate_parsing(StatError &error,
                                                                                 std::ifstream &in_file,
                                                                                 int &line,
                                                                                 unsigned int inb_variable,
                                                                                 double cumul_threshold,
                                                                                 int min_inf_bound)

{
   IidDiscreteMultivariateParametric *dist = NULL;
   DiscreteParametric *marginal = NULL;

   marginal = discrete_parametric_parsing(error, in_file, line, UNIFORM, cumul_threshold, min_inf_bound);

   if (marginal != NULL)
   {
      dist = new IidDiscreteMultivariateParametric(marginal->ident, inb_variable,
                                                   marginal->inf_bound, marginal->sup_bound,
                                                   marginal->parameter, marginal->probability,
                                                   cumul_threshold);
      delete marginal;

   }
   return dist;
}

/*****************************************************************
 *
 *  Default constructor of DiscreteMultivariateDistribution class
 *
 **/

MultinomialCompoundDiscreteParametric::MultinomialCompoundDiscreteParametric()
 : DiscreteMultivariateParametric()
 , param_compound(NULL)
 , cumul_threshold()
{}

/*****************************************************************
 *
 *  Default constructor of MultinomialCompoundDiscreteParametric class
 *
 **/

MultinomialCompoundDiscreteParametric::MultinomialCompoundDiscreteParametric(int inb_variable,
                                                                             const std::vector<double>& iprobability,
                                                                             const DiscreteParametric& dist,
                                                                             double icumul_threshold)
    : DiscreteMultivariateParametric(),
    param_compound(NULL)
{ init(inb_variable, iprobability, dist, icumul_threshold);}

/*****************************************************************
 *
 *  Copy constructor of MultinomialCompoundDiscreteParametric class
 *
 **/

MultinomialCompoundDiscreteParametric::MultinomialCompoundDiscreteParametric(const MultinomialCompoundDiscreteParametric& dist)
{ copy(dist); }



void MultinomialCompoundDiscreteParametric::copy(const MultinomialCompoundDiscreteParametric& dist)
{
    unsigned int i;

    param_marginals.resize(dist.param_marginals.size());

    for(i = 0; i < param_marginals.size(); i++)
        if (dist.param_marginals[i] != NULL)
            param_marginals[i] = new DiscreteParametric(*dist.param_marginals[i]);

    param_compound = new DiscreteParametric(*dist.param_compound);

    nb_variable = dist.nb_variable;
    ident = dist.ident;
    inf_bound = 0;
    sup_bound = dist.sup_bound;
    parameter = dist.parameter;
    probability = dist.probability;

    if ((rparameter != NULL) && (dist.rparameter != NULL))
        *rparameter = *dist.rparameter;
    if ((rparameter == NULL) && (dist.rparameter != NULL))
        rparameter = new std::vector<double>(*dist.rparameter);
}
/*****************************************************************
 *
 *  Copy a MultinomialCompoundDiscreteParametric and return a pointer
 *
 **/

MultinomialCompoundDiscreteParametric* MultinomialCompoundDiscreteParametric::copy_ptr() const
{
   MultinomialCompoundDiscreteParametric *res = new MultinomialCompoundDiscreteParametric(*this);

   return res;
}

/*****************************************************************
 *
 *  Destructor for MultinomialCompoundDiscreteParametric class
 *
 **/

MultinomialCompoundDiscreteParametric::~MultinomialCompoundDiscreteParametric()
{ remove(); }

/*
 * Init function
 * */


void MultinomialCompoundDiscreteParametric::init(int inb_variable,
                                                                             const std::vector<double>& iprobability,
                                                                             const DiscreteParametric& dist,
                                                                             double icumul_threshold)
{
      remove();
      unsigned int i, j, k;
    double s;
    ident = MCOMPOUND_MULTINOMIAL;
    nb_variable = inb_variable;
    assert(iprobability.size() == nb_variable);
    probability = iprobability;
    s = 0;
    for(i = 0; i < nb_variable; i++){
        s += probability[i];
    }
    if(s != 1){
        for(i = 0 ; i < nb_variable-1; i++){
          probability[i] /= s;
        }
        s = 0;
        for(i = 0; i < nb_variable-1; i++){
            s += probability[i];
        }
        probability[nb_variable-1] = 1-s;
    }
    inf_bound = 0;
    cumul_threshold = icumul_threshold;
        param_compound = new DiscreteParametric(dist.ident, dist.inf_bound, dist.sup_bound, dist.parameter, dist.probability);
        parameter.resize(0);
        mean = parametric_mean_computation();
    variance = parametric_variance_computation();
    param_marginals.resize(nb_variable);
    offset.resize(nb_variable);
    nb_value.resize(nb_variable);

    switch(param_compound->ident){
        case BINOMIAL : {
            sup_bound.resize(nb_variable,param_compound->sup_bound);
            for(i = 0; i < nb_variable; i++){
                param_marginals[i] = marginals_compound_multinomial(i);
            }
            break;
        }
        default :
            break;
    }
}
/*****************************************************************
 *
 *  Destructor for IidDiscreteMultivariateParametric class
 *
 **/

void MultinomialCompoundDiscreteParametric::remove()
{
    unsigned int i;
    if(!param_marginals.empty()){
        for(i = 0; i < param_marginals.size(); i++){
            if(param_marginals[i] != NULL){
                delete param_marginals[i];
                param_marginals[i] = NULL;
            }
        }
    }
    if(param_compound != NULL){
        delete param_compound;
        param_compound = NULL;
    }
}

/*****************************************************************
 *
 *  Compute the number of free parameters in DiscreteMultivariateParametric
 *
 **/

unsigned int MultinomialCompoundDiscreteParametric::nb_parameter_computation() const
{
   unsigned int bnbparameter;

   assert(param_compound != NULL);

   return param_compound->nb_parameter_computation() + nb_variable - 1;
}

/*
 * Compute parametric marginal distributions
 * */
void MultinomialCompoundDiscreteParametric::compute_param_marginals(double cumul_threshold)
{
    unsigned int i, j , k;
    std::vector<double> pprobability;
    double s;

    for(i = 0; i < nb_variable; i++){
            param_marginals[i] = marginals_compound_multinomial(i);
        }
}

/*****************************************************************
 *
 *  Print MultinomialCompoundDiscreteParametric parameters
 *  using an output stream
 *
 **/

ostream& MultinomialCompoundDiscreteParametric::ascii_print(ostream &os, bool comment_flag,
                                                            bool exhaustive, bool file_flag, bool cumul_flag,
                                                            const Statiskit::Marginal::Multivariate::CountTable *histo) const
{
   unsigned int i;
   FrequencyDistribution *sum_histo = NULL;
   Statiskit::Marginal::Univariate::Table<int> *sum_mv_histo = NULL;

   header_ascii_print(os);

   os << STAT_TREES_multivariate_distribution_word[ident] << "   ";

   os << STAT_word[STATW_INF_BOUND] << " : " << 0 << "   ";
   os << " ; ";


   os << STAT_word[STATW_PROBABILITY] << " : " << probability;

   os << endl;

   if (exhaustive)
   {
      // write marginals, means, etc.
      ascii_characteristic_print(os, true, comment_flag);


      if (file_flag)
         os << "# ";

      if ((total_nb_value > 0) && (total_nb_value < MULTIVARIATE_DISTRIBUTION_MAX_VALUE))
      {
         mass_computation();
         os << STAT_label[STATL_INFORMATION] << ": " << information_computation() << endl;
      }

      if ((total_nb_value > 0) && (total_nb_value < MULTIVARIATE_DISTRIBUTION_MAX_VALUE))
         // write probabilities
         DiscreteMultivariateDistribution::probability_ascii_print(os, comment_flag, cumul_flag);
   }

   os << STAT_MULTIVARIATE_PARSING_word[TREESTATW_COMPOUNDING_DISTRIBUTION]
      << ": " << endl;

   if (param_compound != NULL)
   {
      if (histo != NULL)
      {
         sum_mv_histo = Stat_trees::Stat_pseudo_histogram_get_sum(*histo, Stat_histogram_value_sum_int);
         sum_histo = Stat_trees::Scalar_stat_histogram_to_frequency_distribution(*sum_mv_histo);
         delete sum_mv_histo;
         sum_mv_histo = NULL;
      }

      param_compound->ascii_print(os);
      if (sum_histo != NULL)
      {
         delete sum_histo;
         sum_histo = NULL;
      }
      if (exhaustive)
      {
         param_compound->Distribution::ascii_print(os, comment_flag, cumul_flag,
                                                   true, sum_histo);
         param_compound->ascii_characteristic_print(os, false, comment_flag);
      }
   }

   return os;
}

/*****************************************************************
 *
 *  Parse a MultinomialCompoundDiscreteParametric from a file
 *  Parsing begins at INF_BOUND: the number of variables
 *  and distribution types must have been parsed before
 *
 **/

MultinomialCompoundDiscreteParametric*
Stat_trees::multinomial_compound_parsing(StatError &error,
                                         std::ifstream &in_file,
                                         int &line,
                                         RWCTokenizer &next,
                                         unsigned int inb_variable,
                                         double cumul_threshold,
                                         int min_inf_bound)
{
   RWLocaleSnapshot locale("en");
   RWCString buffer, token;
   size_t position;
   bool status = true, lstatus;
   register int i, var, correction_line = 0;
   long iinf_bound;
   std::vector<double> vprobability;
   double rparameter = 0.;
   MultinomialCompoundDiscreteParametric *dist = NULL;
   DiscreteParametric *dparam_compound = NULL;

   i = 0;

   if (&next == NULL)
      next = RWCTokenizer(buffer);

   // read a single line
   while ((!((token = next()).isNull())) && (i < 4))
   {
      // read inferior bound

      switch (i)
      {
         case 0 :
         {
            // inferior bound keyword
            if (token != STAT_word[STATW_INF_BOUND])
            {
               status = false;
               error.correction_update(STAT_parsing[STATP_PARAMETER_NAME], STAT_word[STATW_INF_BOUND], line , i + 2);
            }
            break;
          }
         case 1 :
         {
            // test separator ":"
            if (token != ":")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 2);
            }
            break;
         }
         // value of inf. bound
         case 2 :
         {
            lstatus = locale.stringToNum(token , &iinf_bound);
            if ((lstatus) && (iinf_bound != 0))
            {
               lstatus = false;
               error.update(STAT_error[STATR_MIN_INF_BOUND], line, i + 5 + var);
               ostringstream correction_message;
               correction_message << "== " << token;
               error.correction_update(STAT_error[STATR_MIN_INF_BOUND], (correction_message.str()).c_str(), line);
            }
            break;

         }
         case 3 :
         {
            // test separator ";"
            if (token != ";")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 2);
            }
            break;
         }
      }
      i++;
   } // end while

   // i == 4 at this point, in principle
   if (i != 4)
      status = false;

   if (!status)
      error.update(STAT_parsing[STATP_FORMAT] , line);

   i = 2;

   vprobability.resize(inb_variable);

   // test for parameter name
   // read a single line
   while (status && (!((token = next()).isNull())))
   {
      switch ((i - 1) % 4) // block "prob : val1  valK ;"
      {
         case 0 :
         {
            if (token != STAT_word[STATW_PROBABILITY])
            {
               status = false;
               error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PROBABILITY] , line , i + 5);
            }
            break;
         }

         // test separator ":"
         case 1 :
         {
            if (token != ":")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5);
            }
            break;
         }

         // probability values
         case 2 :
         {
            var = 0;
            lstatus = lstatus && locale.stringToNum(token , &vprobability[var]);
            if (!((vprobability[var] > 0.) && (vprobability[var] < 1.)))
            {
               lstatus = false;
               error.update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], line, i + 5 + var);
               ostringstream correction_message;
               correction_message << " != " << token;
               error.correction_update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], (correction_message.str()).c_str(), line);
            }
            rparameter += vprobability[var];
            var++;
            while ((!((token = next()).isNull())) && (var < inb_variable))
            {
               lstatus = lstatus && locale.stringToNum(token , &vprobability[var]);
               if (!((vprobability[var] > 0.) && (vprobability[var] < 1.)))
               {
                  lstatus = false;
                  error.update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], line, i + 5 + var);
                  ostringstream correction_message;
                  correction_message << "!= " << token;
                  error.correction_update(STAT_TREES_parsing[TREESTATP_BAD_PROBABILITY], (correction_message.str()).c_str(), line);
               }
               rparameter += vprobability[var];
               var++;
            }

            if (lstatus)
            {
               if (abs(1-rparameter) > 1.-cumul_threshold)
               {
                  lstatus = false;
                  ostringstream correction_message;
                  correction_message << "in (0, 1)";
                  error.correction_update(STAT_TREES_parsing[TREESTATP_BAD_ONE_MINUS_PROBABILITIES],
                                          (correction_message.str()).c_str(), line);
               }
            }
            // test separator ";"
            // may need to be skipped for discrete_parametric_parsing
            if (token != ";")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5);
            }
            i++;
            break;
         }

         // test separator ";"
         // may need to be skipped for discrete_parametric_parsing
         case 3 :
         {
            if (token != ";")
            {
               status = false;
               error.update(STAT_parsing[STATP_SEPARATOR], line, i + 5 + correction_line);
            }
            break;
         }
      } // end switch ((i - 1) % 4)

      if (!lstatus)
         status = false;
      i++;
   } // end while

   if (i != 5)
   {
      // ; PROBABILITY : ... ;
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line, i + 4 + inb_variable);
   }

   // check keyword COMPOUNDING_DISTRIBUTION
   while (buffer.readLine(in_file, false))
   {
      // skip lines with comments
      line++;

      #     ifdef DEBUG
      cout << line << "  " << buffer << endl;
      #     endif
      position = buffer.first('#');
      if (position != RW_NPOS)
         buffer.remove(position);

      i = 0;

      RWCTokenizer next(buffer);

      // read a single line
      while (status && (!((token = next()).isNull())))
      {
         if ((i == 0) && (token != STAT_MULTIVARIATE_PARSING_word[TREESTATW_COMPOUNDING_DISTRIBUTION]))
         {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_MULTIVARIATE_PARSING_word[TREESTATW_COMPOUNDING_DISTRIBUTION],
                                    line , i);
         }
         if ((i == 1) && (token != ":"))
         {
            status = false;
            error.update(STAT_parsing[STATP_SEPARATOR], line, i);
         }
         i++;
      }
      if ((i > 0 ) &&  (i != 2))
      {
         // ; PROBABILITY : ... ;
         status = false;
         error.update(STAT_parsing[STATP_FORMAT] , line);
      }
      break;
   }


   dparam_compound = discrete_parametric_parsing(error , in_file , line ,
                                                 UNIFORM, cumul_threshold);

   if (dparam_compound == NULL)
      status = false;

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

   if (status)
      dist = new MultinomialCompoundDiscreteParametric(inb_variable, vprobability,
                                                       *dparam_compound,
                                                       cumul_threshold);
   return dist;
}
