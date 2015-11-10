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

#ifndef MV_HISTOGRAM_HPP
#define MV_HISTOGRAM_HPP


/*****************************************************************
 *
 * Set of (different) values of a Statiskit::Table
 *
 **/

template<typename V, typename W> std::vector<V>
Stat_pseudo_histogram_get_elements(const Statiskit::Marginal::Table<V, W>& hist)
{
   int i = 0;
   std::vector<V> res;
   typename Statiskit::Marginal::Table<V,W>::const_iterator it, end;

   it = hist.cbegin();
   end = hist.cend();

   res.resize(hist.get_cardinal());
   for (i=0; i < res.size(); i++)
   {
      // it++;
      res[i] = boost::get<V>((*it).first);
      it++;
   }

   return res;
}

/*!
 *
 * \fn template<typename V, typename W> Statiskit::Marginal::Univariate::Table<W>*
 * Stat_pseudo_histogram_get_sum(const Statiskit::Marginal::Multivariate::Table<W>& hist,
 *                               V (*sum) (const Statiskit::Marginal::Multivariate::CountTable::value_type&))
 * \param hist histogram used to compute the sum of its components
 * to get their histogram.
 * \param *sum function used to compute the sum
 * \details: Histogram of sum of components
 *
 **/

template<typename V, typename W> Statiskit::Marginal::Univariate::Table<W>*
Stat_pseudo_histogram_get_sum(const Statiskit::Marginal::Multivariate::Table<W>& hist,
                              V (*sum) (const Statiskit::Marginal::Multivariate::CountTable::key_type&))
{
   int i = 0;
   Statiskit::Marginal::Univariate::Table<W>* res;
   typename Statiskit::Marginal::Multivariate::event_type vec;
   typename Statiskit::Marginal::Multivariate::Table<W>::const_iterator it, end;

   it = hist.cbegin();
   end = hist.cend();

   res = new Statiskit::Marginal::Univariate::Table<W>(hist.get_storage(0));
   // res = new Statiskit::Marginal::Univariate::Table<W>();

   while(it != end)
   {
      vec = (*it).first;
      res->add(sum(vec), (*it).second);
      it++;
   }

   return res;

}


/*****************************************************************
 *
 * Conversion from std::vector<std::vector<> > to Vectors
 *
 **/

template<typename V>
stat_tool::Vectors* std_vector_to_vectors(const std::vector<std::vector<V> >& svec)
{
   const int inb_vector = svec.size();
   int inb_variable = 0, v, i;
   int *itype = NULL, *iidentifier = NULL, **iint_vector = NULL;
   double **ireal_vector = NULL;
   std::vector<V> val;
   stat_tool::Vectors *res = NULL;

   if (inb_vector > 0)
   {
      iidentifier = new int[inb_vector];
      for(v=0; v < inb_vector; v++)
         iidentifier[v] = v;

      inb_variable = svec[0].size();
      itype = new int[inb_variable];
      for(i=0; i < inb_variable; i++)
	itype[i] = stat_tool::INT_VALUE;

      iint_vector = new int*[inb_vector];
      for(v=0; v < inb_vector; v++)
      {
         iint_vector[v] = new int[inb_variable];
         val = svec[v];
         for(i=0; i < inb_variable; i++)
            iint_vector[v][i] = (int)val[i];
      }
   }

   res = new stat_tool::Vectors(inb_vector, iidentifier, inb_variable, itype, iint_vector, ireal_vector);

   if (inb_vector > 0)
   {
      delete [] iidentifier;
      delete [] itype;
      for(v=0; v < inb_vector; v++)
         delete [] iint_vector[v];
      delete [] iint_vector;
   }
   return res;
}

/*****************************************************************
 *
 * Sum of elements in Statiskit::Marginal::Multivariate::CountTable::value_type
 *
 **/

template<typename V> V
Stat_histogram_value_sum(const Statiskit::Marginal::Multivariate::CountTable::key_type& stat_value)
{
   V res = 0;
   int i;

   for(i=0; i < stat_value.size(); i++)
      res += boost::get<V>(stat_value[i]);

   return res;
}

#endif
