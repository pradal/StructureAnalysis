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


#include "statiskit/core/data/marginal/multivariate.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"

#include "mv_histogram_tools.h"

#include "multivariate_distribution.h"

using namespace stat_tool;

/*****************************************************************
 *
 *
 * Some notes on how to use statistic::Histograms
 *
 **/

/* Frequency of a given value
   Marginal::Univariate::CountTable::const_iterator it;
   it = hist.cfind(value);
   (*it).second;
*/

/* Number of different elements
   hist.get_nb_different_values();
*/

/* Total number of elements (sample size)
   hist.get_nb_values();
*/

/* Use marginalizing
CountTable must be initialized with
PseudoHistogram(scalars_type(nb, type)), ??
e.g. PseudoHistogram(scalars_type(nb, INTEGER)).
*/

/* Also use the following functions to find or list std::vector values */

/*****************************************************************
 *
 *  Conversion between Statiskit::Marginal::Univariate::CountTable
 *  and Stat_tool::FrequencyDistribution
 *
 **/
FrequencyDistribution*
Stat_trees::Scalar_stat_histogram_to_frequency_distribution(const Statiskit::Marginal::Univariate::Table<int>& hist)
{
   int i;
   const Statiskit::Marginal::Univariate::event_type max_val_event = hist.get_max();
   const int max_val = boost::get<int>(max_val_event);
   Statiskit::Marginal::Univariate::Table<int>::const_iterator it;
   Statiskit::Marginal::Univariate::Table<int>::key_type value;

   FrequencyDistribution *res = new FrequencyDistribution;

   res->init(max_val+1);

   for(i=0; i <= max_val; i++)
   {
      it = hist.cfind(i);
      res->frequency[i] += (int)(*it).second;
   }

   res->nb_element_computation();
   res->nb_value_computation();
   res->offset_computation();
   res->max_computation();
   res->mean_computation();
   res->variance_computation();

   return res;
}

/*****************************************************************
 *
 * Conversion from Statistic::Marginal::Vector::Histogram::key_type
 * to std::vector<unsigned int>
 *
 **/

std::vector<unsigned int>
Stat_trees::Stat_histogram_value(const Statiskit::Marginal::Multivariate::CountTable::key_type& stat_value)
{
   int pos;
   std::vector<unsigned int> res_tree_stat;

   if (stat_value.size() > 0)
   {
      res_tree_stat.resize(stat_value.size());
      for(pos=0; pos < stat_value.size(); pos++)
         res_tree_stat[pos] = boost::get<int>(stat_value[pos]);
   }

   return res_tree_stat;
}

/*****************************************************************
 *
 * Conversion from std::vector<unsigned int>
 * to Statistic::Marginal::Multivariate::CountTable::key_type
 *
 **/

Statiskit::Marginal::Multivariate::CountTable::key_type
Stat_trees::Stat_histogram_value(const std::vector<unsigned int>& value)
{
   int pos;
   Statiskit::Marginal::Multivariate::CountTable::key_type res_stat;

   if (value.size() > 0)
   {
      res_stat.resize(value.size());
      for(pos=0; pos < value.size(); pos++)
         res_stat[pos] = (int)value[pos];
   }

   return res_stat;
}

/*****************************************************************
 *
 * Set of (different) values of a Statistic::Marginal::Vector::Histogram
 *
 **/


std::vector<std::vector<unsigned int> >
Stat_trees::Stat_histogram_get_elements(const Statiskit::Marginal::Multivariate::CountTable& hist)
{
   int v, pos;
   std::vector<Statiskit::Marginal::Multivariate::CountTable::key_type> res_stat;
   std::vector<std::vector<unsigned int> > res_tree_stat;
   std::vector<unsigned int> vec;

   res_stat = Stat_trees::Stat_pseudo_histogram_get_elements<Statiskit::Marginal::Multivariate::CountTable::key_type, int>(hist);

   if (res_stat.size() > 0)
   {
      vec.resize(res_stat[0].size());
      res_tree_stat.resize(res_stat.size());
      for(v=0; v < res_stat.size(); v++)
         res_tree_stat[v] = Stat_histogram_value(res_stat[v]);
   }

   return res_tree_stat;
}

/*****************************************************************
 *
 * Multiset of values (with replicates)
 * of a Statistic::Marginal::Vector::Histogram
 *
 **/


std::vector<std::vector<unsigned int> >
Stat_trees::Stat_histogram_get_elements_with_replicates(const Statiskit::Marginal::Multivariate::CountTable& hist)
{
   int v, r, pos;
   std::vector<Statiskit::Marginal::Multivariate::CountTable::key_type> res_stat;
   std::vector<std::vector<unsigned int> > res_tree_stat;
   std::vector<unsigned int> value;

   res_stat = Stat_trees::Stat_pseudo_histogram_get_elements<Statiskit::Marginal::Multivariate::CountTable::key_type, int>(hist);

   if (res_stat.size() > 0)
   {
      for(v=0; v < res_stat.size(); v++)
      {
         value = Stat_histogram_value(res_stat[v]);
         for(r=0; r < Stat_trees::Stat_histogram_get_frequency(hist, value); r++)
            res_tree_stat.push_back(Stat_histogram_value(res_stat[v]));
      }
   }

   return res_tree_stat;
}


/*****************************************************************
 *
 * Get frequency for Statistic::Statiskit::Marginal::Multivariate::CountTable
 *
 **/

int Stat_trees::Stat_histogram_get_frequency(const Statiskit::Marginal::Multivariate::CountTable& hist,
                                             const std::vector<unsigned int>& value)
{
   int freq = 0, pos;
   Statiskit::Marginal::Multivariate::CountTable::key_type stat_value;
   Statiskit::Marginal::Multivariate::CountTable::const_iterator it;

   if (value.size() > 0)
   {
      stat_value.resize(value.size());
      for(pos=0; pos < stat_value.size(); pos++)
         stat_value[pos] = (int)value[pos];
      it = hist.cfind(stat_value);
      freq = (*it).second;
   }

   return freq;
}


/*****************************************************************
 *
 * Print Statistic::Statiskit::Marginal::Multivariate::CountTable
 *
 **/

std::ostream& Stat_trees::print_stat_histogram(const Statiskit::Marginal::Multivariate::CountTable& hist, std::ostream& os)
{
   int i, freq, var, nb_variable = 0;
   const std::vector<std::vector<unsigned int> > values
         = Stat_trees::Stat_histogram_get_elements(hist);

   if (values.size() > 0)
   {
      nb_variable = values[0].size();

      os << nb_variable << " " << STAT_word[STATW_VARIABLES] << endl;

      os << "frequencies (" << hist.get_cardinal() << " values) : " << endl;

      for(i=0; i < values.size(); i++)
      {
         freq = Stat_trees::Stat_histogram_get_frequency(hist, values[i]);
         os << "value: ";
         for(var=0; var < nb_variable; var++)
            os << values[i][var] << "\t";
         os << "; frequency: " << freq << endl;

      }
   }

   return os;
}

/*****************************************************************
 *
 * Conversion fromMarginal::Multivariate::CountTable
 * to DiscreteMultivariateReestimation
 *
 **/

Stat_trees::DiscreteMultivariateReestimation<int>*
Stat_trees::stat_histogram_to_DiscreteMultivariateReestimation(const Statiskit::Marginal::Multivariate::CountTable& hist)
{
   const std::vector<std::vector<unsigned int> > uvalues = Stat_histogram_get_elements_with_replicates(hist);
   Vectors *values;
   DiscreteMultivariateReestimation<int> *res = NULL;

   values = std_vector_to_vectors<unsigned int>(uvalues);
   assert (values != NULL);

   res = new DiscreteMultivariateReestimation<int>(*values);

   delete values;

   return res;
}

