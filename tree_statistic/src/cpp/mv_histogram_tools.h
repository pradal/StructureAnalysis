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
 *       $Id: markov_out_tree.h 2722 2010-12-14 14:17:56Z jbdurand $
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

#ifndef MV_HISTOGRAM_H
#define MV_HISTOGRAM_H


namespace Stat_trees
{

/*! \file markov_out_tree.h
    \brief Purpose:
     provide functions for easy use of Statiskit::Marginal::PseudoHistogram
*/


class DiscreteMultivariateParametric;
template <typename Type> class DiscreteMultivariateReestimation;

/****************************************************************
 *
 *  Functions (tools):
 */

/** Conversion between Statiskit::Marginal::Scalar::Histogram
and Stat_tool::FrequencyDistribution */

stat_tool::FrequencyDistribution*
Scalar_stat_histogram_to_frequency_distribution(const Statiskit::Marginal::Univariate::Table<int>& hist);

/** Sum of elements in Statiskit::Marginal::Multivariate::CountTable::key_type */

template<typename V> V
Stat_histogram_value_sum(const Statiskit::Marginal::Multivariate::CountTable::key_type& stat_value);

static int (*Stat_histogram_value_sum_int)(const Statiskit::Marginal::Multivariate::CountTable::key_type& stat_value) = &Stat_histogram_value_sum<int>;
static double (*Stat_histogram_value_sum_double)(const Statiskit::Marginal::Multivariate::CountTable::key_type& stat_value) = &Stat_histogram_value_sum<double>;

/** Histogram of sum of components */

template<typename V, typename W> Statiskit::Marginal::Univariate::Table<W>*
Stat_pseudo_histogram_get_sum(const Statiskit::Marginal::Multivariate::Table<W>& hist,
                              V (*f) (const Statiskit::Marginal::Multivariate::CountTable::key_type&));

/** Set of (different) values of a Marginal::PseudoHistogram */

template<typename V, typename W> std::vector<V>
Stat_pseudo_histogram_get_elements(const Statiskit::Marginal::Table<V,W>& hist);

/** Set of (different) values of a Marginal::Vector::Histogram */

std::vector<std::vector<unsigned int> >
Stat_histogram_get_elements(const Statiskit::Marginal::Multivariate::CountTable& hist);

/** Multiset of values of a Statiskit::Marginal::Vector::Histogram */

std::vector<std::vector<unsigned int> >
Stat_histogram_get_elements_with_replicates(const Statiskit::Marginal::Multivariate::CountTable& hist);

/** Conversion from Statiskit::Marginal::Vector::Histogram::key_type
  to std::vector<unsigned int> */

std::vector<unsigned int>
Stat_histogram_value(const Statiskit::Marginal::Multivariate::CountTable::key_type& stat_value);

/** Conversion from std::vector<unsigned int> to
    Statiskit::Marginal::Multivariate::CountTable::key_type */

Statiskit::Marginal::Multivariate::CountTable::key_type
Stat_histogram_value(const std::vector<unsigned int>& value);

/** Get frequency for Marginal::Multivariate::CountTable */

int Stat_histogram_get_frequency(const Statiskit::Marginal::Multivariate::CountTable& hist,
                                 const std::vector<unsigned int>& value);

/** Print Statiskit::Marginal::Multivariate::CountTable */

std::ostream& print_stat_histogram(const Statiskit::Marginal::Multivariate::CountTable& hist,
                                   std::ostream& os);

/** Conversion from std::vector<std::vector<int> > to Vectors */

template<typename V>
stat_tool::Vectors* std_vector_to_vectors(const std::vector<std::vector<V> >& svec);

/** Conversion from Statiskit::Marginal::Multivariate::CountTable to DiscreteMultivariateReestimation */

DiscreteMultivariateReestimation<int>*
stat_histogram_to_DiscreteMultivariateReestimation(const Statiskit::Marginal::Multivariate::CountTable& hist);


#include "mv_histogram_tools.hpp"


}; // end namespace

#endif
