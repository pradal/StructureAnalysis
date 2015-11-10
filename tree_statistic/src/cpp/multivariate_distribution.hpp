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

#ifndef MULTIVARIATE_DISTRIBUTION_HPP
#define MULTIVARIATE_DISTRIBUTION_HPP

#include <assert.h>

/*****************************************************************
 *
 *  Default constructor of DiscreteMultivariateReestimation class
 *
 **/

template<typename Type>
DiscreteMultivariateReestimation<Type>::DiscreteMultivariateReestimation()
  : nb_variable(0)
  , min_offset(0)
  , max_alloc_nb_value(0)
  , nb_element(0)
  , frequency(NULL)
  , decode(NULL)
  , marginals(NULL)
  , sums(NULL)
{}

/*****************************************************************
 *
 *  Constructor of DiscreteMultivariateReestimation class
 *  using the number of variables, offset and number of allocated values
 *
 **/

template<typename Type>
DiscreteMultivariateReestimation<Type>::DiscreteMultivariateReestimation(unsigned int inb_variable,
                                                                         unsigned int imin_offset,
                                                                         unsigned int imax_alloc_nb_value)
  : nb_variable(0)
  , min_offset(0)
  , max_alloc_nb_value(0)
  , nb_element(0)
  , frequency(NULL)
  , code(NULL)
  , decode(NULL)
  , marginals(NULL)
  , sums(NULL)
{ init(inb_variable, imin_offset, imax_alloc_nb_value); }

/*****************************************************************
 *
 *  Constructor of DiscreteMultivariateReestimation class
 *  using a stat_tool::Vectors instance
 *
 **/

template<typename Type>
DiscreteMultivariateReestimation<Type>::DiscreteMultivariateReestimation(const stat_tool::Vectors& vec)
  : nb_variable(0)
  , min_offset(0)
  , max_alloc_nb_value(0)
  , nb_element(0)
  , frequency(NULL)
  , code(NULL)
  , decode(NULL)
  , marginals(NULL)
  , sums(NULL)
{
   const int nb_vector = vec.get_nb_vector(),
             vec_nb_variable = vec.get_nb_variable();
   unsigned int i, v, j;
   int imin;
   bool *types = NULL, loop = true;
   value val;

   assert(nb_vector > 0);
   nb_variable = 0;

   types = new bool[vec_nb_variable];

   for(i = 0; i < vec_nb_variable; i++)
      if ((vec.get_type(i) == stat_tool::INT_VALUE) || (vec.get_type(i) == stat_tool::STATE))
      {
         nb_variable++;
         types[i] = true;
      }
      else
         types[i] = false;

   assert(nb_variable > 1);

   i = 0;
   loop = true;
   while(loop)
   {
      if (types[i])
      {
         min_offset = (int)vec.get_min_value(i);
         loop = false;
      }
      if (i >= vec_nb_variable)
         loop = false;
      i++;
   }
   if (i < vec_nb_variable)
   {
      loop = true;
      while(loop)
      {
         if (types[i])
         {
            imin = (int)vec.get_min_value(i);
            assert(imin >= 0);
            min_offset = min(min_offset, (unsigned int)imin);
         }
         i++;
         if (i >= vec_nb_variable)
            loop = false;
      }
   }

   max_alloc_nb_value = 0;
   for(i = 0; i < vec_nb_variable; i++)
      max_alloc_nb_value = max(max_alloc_nb_value, (unsigned int)vec.get_max_value(i));

   max_alloc_nb_value = max_alloc_nb_value - min_offset + 1;

   init(nb_variable, min_offset, max_alloc_nb_value);
   // compute histogram

   val.resize(nb_variable);
   for(v = 0; v < nb_vector; v++)
   {
      j = 0;
      for(i = 0; i < vec_nb_variable; i++)
      {
         if (types[i])
            val[j] = vec.get_int_vector(v, i);
            j++;
      }
      update(val, 1);
   }
   frequency_distribution_computation();
   sum_computation();
}

/*****************************************************************
 *
 *  Change properties of DiscreteMultivariateDistribution:
 *  number of variables, offset and number of values
 *
 **/

template<typename Type> void DiscreteMultivariateReestimation<Type>::init(unsigned int inb_variable,
                                                                          unsigned int imin_offset,
                                                                          unsigned int imax_alloc_nb_value)
{
   unsigned int v;

   remove();

   nb_variable = inb_variable;
   min_offset = imin_offset;
   max_alloc_nb_value = imax_alloc_nb_value;
   nb_element = 0;

   // init frequencies
   frequency = new values_dict();

   // init code
   code = new values_code();

   // init decode
   decode = new values_decode();

   // init marginal frequencies
   marginals = new std::vector<stat_tool::Reestimation<Type> *>(nb_variable);
   for(v = 0; v < nb_variable; v++)
      (*marginals)[v] = new stat_tool::Reestimation<Type>(min_offset + max_alloc_nb_value);
   sums = new std::vector<Type>(nb_variable, 0);
}

/*****************************************************************
 *
 *  Copy constructor of DiscreteMultivariateReestimation class
 *
 **/

template<typename Type>
DiscreteMultivariateReestimation<Type>::DiscreteMultivariateReestimation(const DiscreteMultivariateReestimation<Type>& dist)
{ copy(dist); }


/*****************************************************************
 *
 *  Destructor of DiscreteMultivariateReestimation class
 *
 **/

template<typename Type>
DiscreteMultivariateReestimation<Type>::~DiscreteMultivariateReestimation()
{ remove(); }

/*****************************************************************
 *
 *  Get number of variables of DiscreteMultivariateReestimation
 *
 **/

template<typename Type> unsigned int
DiscreteMultivariateReestimation<Type>::get_nb_variable() const
{ return nb_variable; }

/*****************************************************************
 *
 *  Get elements contained in DiscreteMultivariateReestimation
 *  (without their weights)
 *
 **/

template<typename Type> typename std::vector<std::vector<unsigned int> >
DiscreteMultivariateReestimation<Type>::get_elements() const
{
   std::vector<value> res;
   values_decode::iterator it;

   if (frequency != NULL)
      for(it = decode->begin(); it != decode->end(); it++)
         res.push_back((*decode)[(*it).first]);

   return res;
}

/*****************************************************************
 *
 *  Get marginal of DiscreteMultivariateReestimation class
 *  (return a pointer; object should not be deallocated)
 *
 **/

template<typename Type> const stat_tool::Reestimation<Type>*
DiscreteMultivariateReestimation<Type>::get_marginal_ptr(unsigned int ivariable) const
{
   assert(ivariable < nb_variable);
   assert((marginals != NULL) && (marginals->size() == nb_variable));

   return (*marginals)[ivariable];
}

/*****************************************************************
 *
 *  Get marginal stat_tool::FrequencyDistribution for DiscreteMultivariateReestimation class
 *  (a new instance is allocated)
 *
 **/

template<typename Type> stat_tool::FrequencyDistribution*
DiscreteMultivariateReestimation<Type>::get_marginal_frequency(unsigned int ivariable) const
{
   stat_tool::FrequencyDistribution *res = NULL;
   const stat_tool::Reestimation< Type > *ptr_res = get_marginal_ptr(ivariable);

   if (ptr_res != NULL)
   {
      res = new stat_tool::FrequencyDistribution();
      res->copy(*ptr_res);
   }
   return res;
}


/*****************************************************************
 *
 *  Get sums of DiscreteMultivariateReestimation class
 *
 **/
template<typename Type> typename std::vector<Type>
DiscreteMultivariateReestimation<Type>::get_sums() const
{
   assert(sums != NULL);
   return sums;
}

/*****************************************************************
 *
 *  Get offset value of DiscreteMultivariateReestimation
 *
 **/

template<typename Type> unsigned int
DiscreteMultivariateReestimation<Type>::get_offset() const
{ return min_offset; }

/*****************************************************************
 *
 *  Get maximal number of values of DiscreteMultivariateReestimation
 *
 **/

template<typename Type> unsigned int
DiscreteMultivariateReestimation<Type>::get_max_alloc_nb_value() const
{ return max_alloc_nb_value; }

/*****************************************************************
 *
 *  Return the sum of every values in dictionary
 *
 **/

template<typename Type> Type
DiscreteMultivariateReestimation<Type>::get_nb_element() const
{ return nb_element; }


/*****************************************************************
 *
 *  Compute and return the sum of every values in dictionary
 *
 **/

template<typename Type> Type
DiscreteMultivariateReestimation<Type>::nb_element_computation() const
{
   nb_element = 0;
   iterator it;

   if (frequency != NULL)
      for(it = frequency->begin(); it != frequency->end(); it++)
         nb_element += (*frequency)[(*it).second];

   return nb_element;
}


/*****************************************************************
 *
 *  Return number of elements in dictionary and compute
 *  the number of values in marginals for DiscreteMultivariateReestimation
 *
 **/

template<typename Type> unsigned int
DiscreteMultivariateReestimation<Type>::nb_value_computation() const
{
   unsigned int v;

   if (marginals != NULL)
      for(v = 0; v < nb_variable; v++)
      {
         (*marginals)[v]->nb_value_computation();
         (*marginals)[v]->nb_element_computation();
      }

   if (frequency != NULL)
      return frequency->size();
   else
      return 0;
}

/*****************************************************************
 *
 *  Return offset value and compute offsets in marginals
 *  for DiscreteMultivariateReestimation
 *
 **/

template<typename Type> unsigned int
DiscreteMultivariateReestimation<Type>::offset_computation() const
{
   unsigned int v;

   if (marginals != NULL)
      for(v = 0; v < nb_variable; v++)
         (*marginals)[v]->offset_computation();

   return min_offset;
}

/*****************************************************************
 *
 *  Return maximal number of allocated values and compute maxima
 *  in marginals for DiscreteMultivariateReestimation
 *
 **/
template<typename Type> unsigned int
DiscreteMultivariateReestimation<Type>::max_computation() const
{
   unsigned int v;

   if (marginals != NULL)
      for(v = 0; v < nb_variable; v++)
         (*marginals)[v]->max_computation();

   return max_alloc_nb_value;
}

/*****************************************************************
 *
 *  Compute and return marginal sums for DiscreteMultivariateReestimation
 *  Offset is not substracted
 *
 **/

template<typename Type> typename std::vector<Type>
DiscreteMultivariateReestimation<Type>::sum_computation() const
{
   unsigned int var;
   register int i;

   assert((marginals != NULL) && (sums != NULL));

   sums->resize(nb_variable, 0);

   for(var = 0; var < nb_variable; var++)
   {
      if ((*marginals)[var]->nb_element > 0)
      {
         for (i = min_offset; i < (*marginals)[var]->nb_value; i++)
            (*sums)[var] += (*marginals)[var]->frequency[i] * i;
      }
   }

   return *sums;
}

/*****************************************************************
 *
 *  Compute means of DiscreteMultivariateReestimation
 *
 **/

template<typename Type> void
DiscreteMultivariateReestimation<Type>::mean_computation() const
{
   unsigned int v;

   if (marginals != NULL)
      for(v = 0; v < marginals->size(); v++)
         if ((*marginals)[v] != NULL)
            (*marginals)[v]->mean_computation();
}

/*****************************************************************
 *
 *  Compute variances of DiscreteMultivariateReestimation
 *
 **/

template<typename Type> void
DiscreteMultivariateReestimation<Type>::variance_computation(bool biais) const
{
   unsigned int v;

   if (marginals != NULL)
      for(v = 0; v < marginals->size(); v++)
         if ((*marginals)[v] != NULL)
            (*marginals)[v]->variance_computation();
}


/*****************************************************************
 *
 *  Return means of DiscreteMultivariateReestimation
 *
 **/

template<typename Type> std::vector<double>
DiscreteMultivariateReestimation<Type>::get_means() const
{
   std::vector<double> res(0);
   unsigned int v;

   if (marginals != NULL)
   {
      res.resize(nb_variable);
      for(v = 0; v < nb_variable; v++)
         res[v] = (*marginals)[v]->mean;
   }

   return res;
}

/*****************************************************************
 *
 *  Return variances of DiscreteMultivariateReestimation
 *
 **/

template<typename Type> std::vector<double>
DiscreteMultivariateReestimation<Type>::get_variances() const
{
   std::vector<double> res(0);
   unsigned int v;

   if (marginals != NULL)
   {
      res.resize(nb_variable);
      for(v = 0; v < nb_variable; v++)
         res[v] = (*marginals)[v]->variance;
   }

   return res;
}

/*****************************************************************
 *
 *  Update frequency for a given value in DiscreteMultivariateReestimation
 *  Offset is substracted in code, not in marginals
 *
 **/

template<typename Type> void
DiscreteMultivariateReestimation<Type>::update(const value& v, Type freq)
{
   unsigned int ind, val;
   std::pair<bool, unsigned int> tcode;

   if (frequency != NULL)
   {
      assert(v.size() == nb_variable);
      // compute code in dictionary
      tcode = decoder(v);
      for(ind = 0; ind < nb_variable; ind++)
         (*marginals)[ind]->frequency[v[ind]]++;
      if (tcode.first)
         (*frequency)[tcode.second] += freq;
      else
      {
         // first occurrence of v in *this
         frequency->insert(pair<unsigned int, Type>(tcode.second, freq));
         code->insert(pair<value, unsigned int >(v, tcode.second));
         decode->insert(pair<unsigned int, value >(tcode.second, v));
      }
      nb_element += freq;
   }
}

/*****************************************************************
 *
 *  Get frequency for a given value in DiscreteMultivariateReestimation
 *
 **/

template<typename Type> Type
DiscreteMultivariateReestimation<Type>::get_frequency(const value& v) const
{
   std::pair< bool, unsigned int > pcode;
    unsigned int tcode;

   pcode = decoder(v);
   tcode = pcode.second;

   if (pcode.first)
   {
      assert(tcode >= 0);
      return (*frequency)[tcode];
    } else {
        return 0;
    }
}

/*****************************************************************
 *
 *  Print frequencies and characteristics of
 *  DiscreteMultivariateReestimation
 *
 **/

template<typename Type> std::ostream&
DiscreteMultivariateReestimation<Type>::print(std::ostream& os) const
{
   unsigned int var; // variable
   value v;
   values_decode::iterator it;

   os << nb_variable << " " << stat_tool::STAT_word[stat_tool::STATW_VARIABLES] << endl;

   if (frequency != NULL)
   {
      os << endl;
      for(var = 0; var < nb_variable; var++)
      {
         if (marginals != NULL)
         {
            os << stat_tool::STAT_word[stat_tool::STATW_VARIABLE] << " " << var << ": " << endl;
            (*marginals)[var]->ascii_characteristic_print(os);
         }
      }
   }
   os << endl;
   os << "offset : " << min_offset;
   if (max_alloc_nb_value > 0)
      os << "   maximum : " << max_alloc_nb_value + min_offset - 1;
   os << endl;

   if (frequency != NULL)
   {
      os << "frequencies (" << max_alloc_nb_value << " values) : " << endl;
      for(it = decode->begin(); it != decode->end(); it++)
      {
         v = (*decode)[(*it).first]; // vector
         os << v << ": " << (*frequency)[(*it).first] << endl;
      }
   }
   return os;
}

/*****************************************************************
 *
 *  Matplotlib output for DiscreteMultivariateReestimation class
 *
 **/

template<typename Type> stat_tool::MultiPlotSet*
DiscreteMultivariateReestimation<Type>::get_plotable() const
{
   // plot_set[variable][view]
   register int i, j;
   const int nb_plot_set = nb_variable;
   ostringstream title;
   stat_tool::DiscreteParametricModel *pmodel = NULL;
   stat_tool::FrequencyDistribution *chisto = NULL;
   stat_tool::MultiPlotSet *plot_set = new stat_tool::MultiPlotSet(nb_plot_set), *plot_histo = NULL;
   stat_tool::MultiPlot *mplot = NULL;

   plot_set->border = "15 lw 0";

   title.str("");
   //title << stat_tool::STAT_MULTIVARIATE_label[TREEstat_tool::STATL_DISCRETE_MULTIVARIATE_DISTRIBUTION];
   //title << " " << stat_tool::STAT_label[stat_tool::STATL_FREQUENCY_DISTRIBUTION];

   plot_set->title = title.str();

   for(i = 0; i < nb_plot_set; i++)
   {
      // it is required in stat_tool::DiscreteParametricModel->plotable
      // that this->frequencydistribution != NULL
      if (get_marginal_ptr(i) != NULL)
      {
         chisto = new stat_tool::FrequencyDistribution();
         chisto->copy(*get_marginal_ptr(i));

         plot_histo = chisto->get_plotable();
         (*plot_set)[i].resize((*plot_histo)[0].size());
         for(j = 0; j < (*plot_histo)[0].size(); j++)
            (*plot_set)[i][j] = (*plot_histo)[0][j];
      }
      title.str("");
      title << stat_tool::STAT_label[stat_tool::STATL_VARIABLE] << " " << i+1;
      (*plot_set)[i].title = title.str();
      if (chisto != NULL)
      {
         delete chisto;
         chisto = NULL;
      }
   }

   return plot_set;
}

/*****************************************************************
 *
 *  Compute log likelihood of a given parametric model
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::likelihood_computation(const DiscreteMultivariateDistribution& dist) const
{
   double log_likelihood = 0, mass;
   value v;
   values_decode::iterator it;


#  ifdef DEBUG
   std::cout << "like.." << std::endl;
#  endif
   if (frequency != NULL)
   {
      for(it = decode->begin(); it != decode->end(); it++)
      {
         v = (*decode)[(*it).first]; // vector
#        ifdef DEBUG
         std::cout << v << std::endl;
#        endif
         mass = dist.get_mass(v, true); // compute log mass
#        ifdef DEBUG
         cout << "stat_tool::D_INF: " << stat_tool::D_INF << endl;
         cout << "mass > stat_tool::D_INF: " << (mass > stat_tool::D_INF) << endl;
#        endif
         if (mass > stat_tool::D_INF)
         {
            log_likelihood += mass * (*frequency)[(*it).first];
         }
         else
         {
            return stat_tool::D_INF;
         }
      }
   }
   return log_likelihood;
}

/*****************************************************************
 *
 *  Compute log likelihood of a given parametric model
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::likelihood_computation(const MultinomialCompoundDiscreteParametric& dist) const
{
   double log_likelihood = 0, mass;
   value v;
   values_decode::iterator it;

   if (frequency != NULL)
   {
      for(it = decode->begin(); it != decode->end(); it++)
      {
         v = (*decode)[(*it).first]; // vector
         mass = log(dist.get_mass(v)); // compute log mass
         if (mass > stat_tool::D_INF)
            log_likelihood += mass * (*frequency)[(*it).first];
         else
            return stat_tool::D_INF;
      }
   }
   return log_likelihood;
}

/*****************************************************************
 *
 *  Compute stat_tool::FrequencyDistribution nb_element, mean, etc.
 *  in marginals of DiscreteMultivariateReestimation class
 *
 **/

template<typename Type> void
DiscreteMultivariateReestimation<Type>::frequency_distribution_computation() const
{
   unsigned int s;

   if (marginals != NULL)
      for(s = 0; s < marginals->size(); s++)
         if ((*marginals)[s] != NULL)
         {
            (*marginals)[s]->nb_value_computation();
            (*marginals)[s]->offset_computation();
            (*marginals)[s]->nb_element_computation();
            (*marginals)[s]->max_computation();
            (*marginals)[s]->mean_computation();
            (*marginals)[s]->variance_computation();
         }
}

/*****************************************************************
 *
 *  Compute stat_tool::FrequencyDistributionof the sum of elements in each vector
 *  for DiscreteMultivariateReestimation class
 *
 **/

template<typename Type> stat_tool::FrequencyDistribution*
DiscreteMultivariateReestimation<Type>::get_sum_frequency_distribution() const
{
   std::vector<unsigned int> sum_vec(0); // vector of the sums of elements
   std::vector<Type> freq_sum(0); // vector of the frequency of each sum
   unsigned int i, j, nmax;
   iterator it;
   value val;
   stat_tool::FrequencyDistribution* histo = NULL;

   i = 0;
   nmax = 0;

   for(it = begin(); it != end(); it++)
   {
      val = (*decode)[(*it).first]; // get current value of vector
      sum_vec.push_back(0);
      freq_sum.push_back(get_frequency(val));
      for(j = 0; j < nb_variable; j++) // compute sum of elements
         sum_vec[i] += val[j];
      nmax = MAX(nmax, sum_vec[i]); // compute max value of sums
      i++;
   }

   histo = new stat_tool::FrequencyDistribution(nmax+1);
   for(i = 0; i < sum_vec.size(); i++)
      histo->frequency[sum_vec[i]] += freq_sum[i];

   histo->nb_element_computation();
   histo->offset_computation();
   histo->max_computation();
   histo->mean_computation();
   histo->variance_computation();

   return histo;
}

/*****************************************************************
 * Estimate a DiscreteMultivariateDistribution
 * */
template<typename Type> DiscreteMultivariateDistribution* DiscreteMultivariateReestimation<Type>::distribution_estimation() const
{
    /*
     * TODO en mieux
     * */
    DiscreteMultivariateDistribution *dist = new DiscreteMultivariateDistribution();
    dist->nb_variable = nb_variable;
    dist->nb_value.resize(nb_variable);
    dist->total_nb_value = nb_value_computation();
    for(unsigned int i = 0; i < nb_variable; ++i)
        dist->nb_value[i] = (*marginals)[i]->nb_value;
    dist->nb_parameter = dist->total_nb_value-1;
    dist->offset = std::vector<unsigned int>(nb_variable, 0);
    //dist->maxima = std::vector<unsigned int>(nb_variable, max_computation());
    std::vector<value> values = get_elements();
    for(std::vector<value>::iterator it = values.begin(); it != values.end(); ++it)
    {
        dist->set_mass(*it, get_frequency(*it)/((double)nb_element));
        std::cout << dist->get_mass(*it,false) << std::endl;
    }
    return dist;
}

/*****************************************************************
 *
 *  Estimate a DiscreteMultivariateParametric with given type
 *  using a stat_tool::StatError object, the type, the likelihood (to be updated)
 *  minimal possible inferior bound, flag on estimating this bound
 *  and threshold on cumulative distribution function.
 *
 **/

template<typename Type> DiscreteMultivariateParametric*
DiscreteMultivariateReestimation<Type>::parametric_estimation(stat_tool::StatError &error,
                                                              int ident, double &likelihood,
                                                              int min_inf_bound,
                                                              bool min_inf_bound_flag, double cumul_threshold) const
{
   DiscreteMultivariateParametric *pdist = new DiscreteMultivariateParametric();

   // error.init();
   pdist->ident = ident;

   if ((min_inf_bound < 0) || (min_inf_bound > 1) || (min_inf_bound > min_offset))
      error.update(stat_tool::STAT_error[stat_tool::STATR_MIN_INF_BOUND]);
   else
   {
      switch (pdist->ident)
      {
         case MNEGATIVE_MULTINOMIAL:
            likelihood = negative_multinomial_estimation(pdist, min_inf_bound,
                                                        min_inf_bound_flag, cumul_threshold);
            break;
         case MPOISSON:
            likelihood = mpoisson_estimation_ML(pdist, min_inf_bound,
                                                min_inf_bound_flag, cumul_threshold);
            break;
         case MMULTINOMIAL:
            likelihood = mmultinomial_estimation(pdist, min_inf_bound,
                                                 min_inf_bound_flag, cumul_threshold);
            break;
         default:
            likelihood = stat_tool::D_INF;
            break;
      }


      if (likelihood == stat_tool::D_INF)
      {
         if (pdist != NULL)
         {
            delete pdist;
            pdist = NULL;
            error.update(stat_tool::STAT_error[stat_tool::STATR_ESTIMATION_FAILURE]);
         }
         else
         {
            pdist->nb_parameter_update();
            if (!min_inf_bound_flag)
            // inf_bounds do not count as parameters if these were not estimated
              (pdist->nb_parameter) -= nb_variable;
         }
      }
   }
   return pdist;
}

/*****************************************************************
 *
 *  Estimate both a DiscreteMultivariateParametric and its type
 *  using a stat_tool::StatError object, the type, the likelihood (to be updated)
 *  minimal possible inferior bound, flag on estimating this bound
 *  and threshold on cumulative distribution function.
 *
 **/

template<typename Type> DiscreteMultivariateParametric*
DiscreteMultivariateReestimation<Type>::type_parametric_estimation(stat_tool::StatError &error,
                                                                   double &likelihood,
                                                                   int min_inf_bound,
                                                                   bool min_inf_bound_flag,
                                                                   double cumul_threshold) const
{

   double current_likelihood;
   DiscreteMultivariateParametric *bdist = new DiscreteMultivariateParametric(),
                                  *dist = new DiscreteMultivariateParametric();
   MultinomialCompoundDiscreteParametric *bcdist = new MultinomialCompoundDiscreteParametric();

   // error.init();

   if ((min_inf_bound < 0) || (min_inf_bound > 1) || (min_inf_bound > min_offset))
      error.update(stat_tool::STAT_error[stat_tool::STATR_MIN_INF_BOUND]);
   else
   {

      likelihood = negative_multinomial_estimation(dist, min_inf_bound, min_inf_bound_flag, cumul_threshold);

      if (likelihood != stat_tool::D_INF)
         dist->ident = MNEGATIVE_MULTINOMIAL;


      current_likelihood = mpoisson_estimation_ML(bdist, min_inf_bound, min_inf_bound_flag, cumul_threshold);
      if (current_likelihood > likelihood)
      {
         bdist->ident = MPOISSON;
         likelihood = current_likelihood;
         delete dist;
         dist = new DiscreteMultivariateParametric(*bdist);
      }


      current_likelihood = mmultinomial_estimation(bdist, min_inf_bound, min_inf_bound_flag, cumul_threshold);
      if (current_likelihood > likelihood)
      {
         bdist->ident = MMULTINOMIAL;
         likelihood = current_likelihood;
         delete dist;
         dist = new DiscreteMultivariateParametric(*bdist);
      }

            current_likelihood = compound_multinomial_estimation(bcdist, stat_tool::I_DEFAULT, min_inf_bound, min_inf_bound_flag, cumul_threshold);
      if (current_likelihood > likelihood)
      {
         likelihood = current_likelihood;
         delete dist;
         dist = new MultinomialCompoundDiscreteParametric(*bcdist);
      }

      // etc.


      if (bdist != NULL)
      {
         delete bdist;
         bdist = NULL;
      }
            if(bcdist != NULL){
                delete bcdist;
                bcdist = NULL;
            }

      if (likelihood == stat_tool::D_INF)
      {
         if (dist != NULL)
         {
            delete dist;
            dist = NULL;
            error.update(stat_tool::STAT_error[stat_tool::STATR_ESTIMATION_FAILURE]);
         }
         else
         {
            dist->nb_parameter_update();
            if (!min_inf_bound_flag)
            // inf_bounds do not count as parameters if these were not estimated
              (dist->nb_parameter) -= nb_variable;
         }
      }
   }

   if (dist != NULL)
      dist->ascii_write(std::cout, false);

   return dist;

}

/*****************************************************************
 *
 *  Estimate negative multinomial distribution using marginals only
 *  using minimal possible inferior bound, flag on estimating this bound
 *  and threshold on cumulative distribution function.
 *  Return the likelihood value
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_marginal_estimation(DiscreteMultivariateParametric *pdist, int min_inf_bound,
                                                                                 bool min_inf_bound_flag, double cumul_threshold) const
{
   unsigned int var, minf_bound;
   double norm = 1., // normalization of probabilities
          likelihood = stat_tool::D_INF;
   stat_tool::DiscreteParametric *mdist; // marginal distribution
   std::vector<double> probability, parameter;
   std::vector<unsigned int> sup_bound;

   assert(pdist != NULL);

   if ((frequency != NULL) && (marginals != NULL))
   {
      mdist = new stat_tool::DiscreteParametric((int)((min_offset + max_alloc_nb_value + 1) *  stat_tool::SAMPLE_NB_VALUE_COEFF), stat_tool::NEGATIVE_BINOMIAL);
      minf_bound = (*marginals)[0]->nb_value;
      probability.resize(nb_variable);
      parameter.resize(1);
      // estimate initial parameter and offsets from marginals
      for(var = 0; var < nb_variable; var++)
      {
         (*marginals)[var]->negative_binomial_estimation(mdist, min_inf_bound,
                                                         min_inf_bound_flag, cumul_threshold);
         parameter[0] += mdist->parameter / nb_variable;
         probability[var] = mdist->probability;
         minf_bound = MIN(minf_bound, MIN(mdist->inf_bound, (*marginals)[var]->offset));
         norm += (1-mdist->probability);
      }
      if (parameter[0] > 0.)
      {
         // if the marginals could not be estimated, parameter[0] <= 0
         // estimation should be aborted in this case

         // mass and cumul are always deallocated in ~DiscreteParametric
         delete mdist;
         mdist = NULL;
         if (norm > 0)
         {
            // computation of probabilities
            for(var = 0; var < nb_variable; var++)
               probability[var] =   1 - ((1 - probability[var]) / norm);

         pdist->init(nb_variable, MNEGATIVE_MULTINOMIAL, minf_bound, sup_bound,
                     parameter, probability, cumul_threshold);
         likelihood = likelihood_computation(*pdist);
         }
      }
   }
   return likelihood;
}

/*****************************************************************
 *
 *  Estimate negative multinomial distribution from
 *  DiscreteMultivariateReestimation
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_estimation(DiscreteMultivariateParametric *pdist, int min_inf_bound,
                                                                        bool min_inf_bound_flag, double cumul_threshold) const
{
   unsigned int var, nb_element, inf_bound;
   double likelihood = stat_tool::D_INF;
   Type total_counts = 0;
   std::vector<unsigned int> sup_bound;
   std::vector<double> parameter, probability;
   std::vector<Type> mean_counts(nb_variable);

   // Estimate inferior bound and initialize estimation of n
   // using estimation on marginals

   likelihood = negative_multinomial_marginal_estimation(pdist, min_inf_bound, min_inf_bound_flag, cumul_threshold);
   if (likelihood > stat_tool::D_INF)
   {
      nb_element = (*marginals)[0]->nb_element;
      inf_bound = pdist->get_inf_bound();
      sup_bound = pdist->get_sup_bound();
      parameter = pdist->get_parameter();
      probability = pdist->get_probability();
      if (sums == NULL)
         sum_computation();

      mean_counts = *sums;
      for(var = 0; var < nb_variable; var++)
      {
         assert(nb_element == (*marginals)[var]->nb_element);
         mean_counts[var] -= inf_bound * nb_element;
#        ifdef MESSAGE
         if (mean_counts[var] == 0)
            cout << "Estimation error for negative multinomial distribution:"
                 << " sum of counts is null for variable " << var << endl;
#        endif
         total_counts += mean_counts[var];
      }

      // estimate parameter
      parameter[0] = negative_multinomial_parameter_estimation(*pdist, inf_bound, mean_counts, total_counts, parameter[0]);


      // estimate probabilities
      for(var = 0; var < nb_variable; var++)
         probability[var] = 1- (mean_counts[var] /
                            (parameter[0] * nb_element + total_counts));

      pdist->init(nb_variable, MNEGATIVE_MULTINOMIAL, inf_bound, sup_bound,
                  parameter, probability, cumul_threshold);

      likelihood = likelihood_computation(*pdist);
   }

   return likelihood;
}

/*****************************************************************
 *
 *  Optimization of the log likelihood with respect to parameter
 *  in negative multinomial setting from DiscreteMultivariateReestimation
 *  using data, inferior bound, counts, initial parameter value, tolerance
 *  (relative increase of target function) and maximal number of iterations
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_parameter_estimation(const DiscreteMultivariateParametric& dist,
                                                                                  int inf_bound,
                                                                                  const std::vector<Type>& mean_counts,
                                                                                  Type total_counts,
                                                                                  double initial_parameter,
                                                                                  double tol,
                                                                                  unsigned int nb_iter) const
{
   double p_inf, p_sup, val_inf, val_sup, cur_val;
   const double scale = 10.;
   std::pair<double, double> res1, res2;

   cur_val = negative_multinomial_parameter_estimation_score_target(dist, inf_bound, mean_counts, total_counts, initial_parameter);

   val_inf = negative_multinomial_parameter_estimation_score_target(dist, inf_bound, mean_counts, total_counts, initial_parameter / scale);

   val_sup = negative_multinomial_parameter_estimation_score_target(dist, inf_bound, mean_counts, total_counts, initial_parameter * scale);
   val_sup = stat_tool::D_INF;

   if (cur_val > MIN(val_inf, val_sup))
   {
      p_inf = initial_parameter;
      p_sup = initial_parameter * scale;
      res1 = negative_multinomial_parameter_estimation(dist, inf_bound, mean_counts, total_counts,
                                                       p_sup, p_inf, tol, nb_iter/2);
      p_sup = initial_parameter;
      p_inf = initial_parameter / scale;

      res2 = negative_multinomial_parameter_estimation(dist, inf_bound, mean_counts, total_counts,
                                                       p_sup, p_inf, tol, nb_iter/2);
      if (res1.second > res2.second)
         return res1.first;
      else
         return res2.first;
   }
   else
   {
      p_sup = initial_parameter * scale;
      p_inf = initial_parameter / scale;
      res1 = negative_multinomial_parameter_estimation(dist, inf_bound, mean_counts, total_counts,
                                                       p_sup, p_inf, tol, nb_iter);
      return res1.first;
   }

}

/*****************************************************************
 *
 *  Optimization of the log likelihood with respect to parameter
 *  in negative multinomial setting from DiscreteMultivariateReestimation
 *  using data, inferior bound, counts, two bounds for parameter value, tolerance
 *  (relative increase of target function) and maximal number of iterations
 *  Return parameter value and target function
 *
 **/


template<typename Type> std::pair<double, double>
DiscreteMultivariateReestimation<Type>::negative_multinomial_parameter_estimation(const DiscreteMultivariateParametric& pdist,
                                                                                  int inf_bound,
                                                                                  const std::vector<Type>& mean_counts,
                                                                                  Type total_counts,
                                                                                  double parameter_sup,
                                                                                  double parameter_inf,
                                                                                  double tol,
                                                                                  unsigned int nb_iter) const
{
   bool stop = false;
   unsigned int it = 0; // number of iterations
   double p, p_inf, p_sup, val_inf, val_sup, cur_val, prev_val = stat_tool::D_INF;
   const unsigned int min_iter = 10; // minimal number of iterations
   std::pair<double, double> res;

   p_sup = parameter_sup;
   p_inf = parameter_inf;

   cur_val = negative_multinomial_parameter_estimation_score_target(pdist, inf_bound, mean_counts, total_counts, p_sup);
   val_sup = cur_val;
   val_inf = negative_multinomial_parameter_estimation_score_target(pdist, inf_bound, mean_counts, total_counts, p_inf);

   while (!stop)
   {
      it++;
      p = (p_sup + p_inf) / 2;
      prev_val = cur_val;
      cur_val = negative_multinomial_parameter_estimation_score_target(pdist, inf_bound, mean_counts, total_counts, p);

      if (cur_val < MIN(val_inf, val_sup))
         stop = true;
      else
      {
         if (cur_val > val_sup)
         {
            p_sup = p;
            val_sup = cur_val;
         }
         else
         {
            p_inf = p;
            val_inf = cur_val;
         }
      }
      if (it > nb_iter)
         stop = true;
      if ((it > min_iter) && ((cur_val - prev_val) / std::abs(prev_val) < tol))
         stop = true;
   }
   res.first = p;
   res.second = cur_val;

   return res;
}

/*****************************************************************
 *
 *  Compute term of the log likelihood to optimize with respect to parameter
 *  in negative multinomial setting from DiscreteMultivariateReestimation
 *  using data, inferior bound, counts and parameter value (solve likelihood equations)
 *
 **/

template <typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_parameter_estimation_score_target(const DiscreteMultivariateParametric& dist,
                                                                                               int inf_bound,
                                                                                               const std::vector<Type>& mean_counts,
                                                                                               Type total_counts,
                                                                                               double parameter) const
{
   double sqrd_target = negative_multinomial_parameter_estimation_score_value(dist, inf_bound,
                                                                              mean_counts,
                                                                              total_counts,
                                                                              parameter);

   return -(sqrd_target * sqrd_target);
}

/*****************************************************************
 *
 *  Compute term of the log likelihood to cancel with respect to parameter
 *  in negative multinomial setting from DiscreteMultivariateReestimation
 *  using data, inferior bound, counts and parameter value (solve likelihood equations)
 *
 **/

template <typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_parameter_estimation_score_value(const DiscreteMultivariateParametric& dist,
                                                                                              int inf_bound,
                                                                                              const std::vector<Type>& mean_counts,
                                                                                              Type total_counts,
                                                                                              double parameter) const
{
   unsigned int i, sum_val; // sum of values
   // value of likelihood derivative, to be canceled
   double sqrd_target = 0., digamma_val = 0.;
   value val;
   Type freq;
   iterator it;

   assert((frequency != NULL) && (parameter > 0));

   for(it = begin(); it != end(); it++)
   {
      digamma_val = .0;
      freq = (*it).second;
      val = (*decode)[(*it).first];
      sum_val = 0;
      for(i = 0; i < nb_variable; i++)
         sum_val += (val[i] - inf_bound);
      for(i = 0; i < sum_val; i++)
         digamma_val += 1. / (parameter + i);
      sqrd_target += freq * digamma_val;
   }
   sqrd_target -= nb_element * log(1 +  total_counts / (nb_element * parameter));

   return sqrd_target;
}

/*****************************************************************
 *
 *  Estimate multinomial distribution from
 *  DiscreteMultivariateReestimation by likelihood maximization
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::mmultinomial_estimation(DiscreteMultivariateParametric *pdist,
                                                                int min_inf_bound, bool min_inf_bound_flag,
                                                                double tol) const
{
   bool status;
   stat_tool::DiscreteParametric *mdist;
   unsigned int var, minf_bound, n_cur, n_prev;
   iterator it;
   double norm = 0, likelihood = stat_tool::D_INF;
   std::vector<double> probability, parameter, means;
   std::vector<unsigned int> sup_bound, values;

   assert(pdist != NULL);

   if ((frequency != NULL) && (marginals != NULL))
   {
      if (min_inf_bound_flag)
      {
         mdist = new stat_tool::DiscreteParametric((int)((min_offset + max_alloc_nb_value + 1) *  stat_tool::SAMPLE_NB_VALUE_COEFF), stat_tool::BINOMIAL);
         (*marginals)[0]->binomial_estimation(mdist, min_inf_bound,
                                             min_inf_bound_flag);
         minf_bound = mdist->inf_bound;
         delete mdist;
         for(var = 1; var < nb_variable; var++)
         {
            mdist = new stat_tool::DiscreteParametric((int)((min_offset + max_alloc_nb_value + 1) *  stat_tool::SAMPLE_NB_VALUE_COEFF), stat_tool::BINOMIAL);
            (*marginals)[var]->binomial_estimation(mdist, min_inf_bound, min_inf_bound_flag);
            minf_bound = MIN(minf_bound, mdist->inf_bound);
            delete mdist;
         }
         mdist = NULL;
      }
      else
         minf_bound = min_inf_bound;

      it = begin();
      n_prev = 0;
      values = (*decode)[(*it).first];
      for(var = 0; var < nb_variable; var++)
         n_prev += values[var];

      for(it = it++; it != end(); it++)
      {
         values = (*decode)[(*it).first];
         n_cur = 0;
         for(var = 0; var < nb_variable; var++)
            n_cur += values[var];

         // Prevents the program to stop if the user provides
         // incompatible data (likelihood is 0. then),
         // as opposed to assert(n_cur == n_prev)
         if (n_cur != n_prev)
            return stat_tool::D_INF;

         // inutile : condition toujours verifiee
         // n_prev = n_cur;
      }
      parameter.resize(1, n_cur);
      probability.resize(nb_variable, 0);
      if (parameter[0] > 0.)
      {
         // i.e. the number of trials should be > 0
         mean_computation();
         means = get_means();
         for(var = 0; var < nb_variable; var++)
         {
            probability[var] = (means[var]-minf_bound)/n_cur;
            norm += probability[var];
         }
         if (norm != 1.)
            for(var = 0; var < nb_variable; var++)
               probability[var] /= norm;

         sup_bound.resize(4, parameter[0]-(nb_variable-1)*minf_bound);
         pdist->init(nb_variable, MMULTINOMIAL, minf_bound, sup_bound,
                     parameter, probability, tol);
         likelihood = likelihood_computation(*pdist);
      }
   }

   return likelihood;
}

/*****************************************************************
 *
 *  Estimate compound multinomial distribution from
 *  DiscreteMultivariateReestimation
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::compound_multinomial_estimation(MultinomialCompoundDiscreteParametric *pdist, int iident, int min_inf_bound, bool min_inf_bound_flag, double tol) const
{
    double likelihood = stat_tool::D_INF;
    value val;
    std::vector<double> eprobability(nb_variable,0);
    unsigned int  var;
    stat_tool::Reestimation<Type> *histo = get_sum_frequency_distribution();
    stat_tool::DiscreteParametric *cdist = NULL;

        for(var = 0; var < nb_variable; var++){
            eprobability[var] = (*sums)[var];
        }

        if(iident == stat_tool::I_DEFAULT){
            cdist = histo->type_parametric_estimation(min_inf_bound, min_inf_bound_flag, tol);
        } else {
            cdist = new stat_tool::DiscreteParametric((int)((min_offset + max_alloc_nb_value + 1) * stat_tool::SAMPLE_NB_VALUE_COEFF), iident);
            histo->parametric_estimation(cdist, min_inf_bound, min_inf_bound_flag, tol);
        }

    pdist->init(nb_variable, eprobability, *cdist, tol);
    likelihood = likelihood_computation(*pdist);
    delete histo;

    return likelihood;
}

/*****************************************************************
 *
 *  Compute and return empirical covariance matrix from
 *  DiscreteMultivariateReestimation
 *
 **/

template<typename Type> std::vector< std::vector<double> >
DiscreteMultivariateReestimation<Type>::get_covariances() const
{
    iterator it;
    unsigned int i, j, nb_val;
    value values;
    double res;
    std::vector<double> means;
    std::vector<double> vars;
    std::vector< std::vector<double> > VarCovar(nb_variable, std::vector<double>(nb_variable, 0));
    mean_computation();
    variance_computation();
    means = get_means();
    vars = get_variances();
    nb_val = get_nb_element();

    for(it = begin(); it != end(); it++){
        values = (*decode)[(*it).first];
        for(i = 0; i < nb_variable-1; i++){
            for(j = i+1; j < nb_variable; j++){
                res = 1/((double) nb_val) * (values[i]-means[i]) * (values[j]-means[j]);
                VarCovar[i][j] += res;
                VarCovar[j][i] += res;
            }
        }
    }
    for(i = 0; i < nb_variable; i++){
        VarCovar[i][i] = vars[i];
    }

    return VarCovar;
}

/*****************************************************************
 *
 *  Estimate multivariate Poisson distribution from
 *  DiscreteMultivariateReestimation by likelihood maximization
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::mpoisson_estimation_ML(DiscreteMultivariateParametric* pdist, int min_inf_bound,
                                                               bool min_inf_bound_flag, double tol) const
{
    unsigned int var, minf_bound, max_it = 100000, min_it = 10, i = 0, j, k, min;
    double likelihood = stat_tool::D_INF, cur_val, prev_val, s, s2, esp;
    std::vector<double> log_p, parameter, probability;
    std::vector<unsigned int> sup_bound;
    value values;
    iterator it;
    std::vector<double> means(nb_variable, 0);
    bool loop = true;
        DiscreteMultivariateParametric *cdist = NULL;

    likelihood = mpoisson_estimation_M(pdist, min_inf_bound, min_inf_bound_flag, tol);

    if (likelihood > stat_tool::D_INF)
    {
       minf_bound = pdist->inf_bound;
       parameter = pdist->parameter;
             if(parameter[0] == 0)
                 parameter[0] = 0.01;

       while(loop){
           s2 = 0;
           for(it = begin(); it != end(); it++){
               values = (*decode)[(*it).first];
               min = values[0]-minf_bound;
               for(var = 1; var < nb_variable; var++){
                   min = MIN(min, values[var]-minf_bound);
               }
               esp = 0;
               s = 1;
               j = 1;
               log_p.resize(min+1,0);
               while(j <= min){
                   log_p[j] = std::log(parameter[0]) - std::log((double)j);
                   for(k = 0; k < nb_variable; k++){
                       log_p[j] += std::log(values[k]-minf_bound-j+1)- std::log(parameter[k+1]);
                   }
                   log_p[j] += log_p[j-1];
                   s += exp(log_p[j]);
                   j++;
               }
               for(j = 1; j <= min; j++){
                   esp += j*exp(log_p[j])/s;
               }
               s2 += get_frequency(values) * esp;
           }
           parameter[0] = s2/((double)get_nb_element());
           mean_computation();
           means = get_means();
           for(var = 0; var < nb_variable; var++){
               parameter[var+1] = means[var] - minf_bound - parameter[0];
           }
           i++;
           if(i > min_it){
               prev_val = likelihood_computation(*pdist);
           }
                     pdist->init(nb_variable, MPOISSON, minf_bound, sup_bound, parameter, probability, tol);
           if(i > min_it){
               cur_val = likelihood_computation(*pdist);
               if((cur_val-prev_val)/(std::abs(prev_val)) < tol)
                   loop = false;
           }
                     if(i > max_it)
                         loop = false;
       }
       likelihood = cur_val;
    }

    return likelihood;
}

/*****************************************************************
 *
 *  Estimate multivariate Poisson distribution from
 *  DiscreteMultivariateReestimation by Moments
 *
 **/

template<typename Type> double
DiscreteMultivariateReestimation<Type>::mpoisson_estimation_M(DiscreteMultivariateParametric *pdist, int min_inf_bound, bool min_inf_bound_flag, double cumul_threshold) const
{
     bool status = true;
     unsigned int var, var2, minf_bound, nb_val = 1;
     stat_tool::DiscreteParametric *mdist;
     double s, likelihood = stat_tool::D_INF;
     std::vector<double> means(nb_variable,0), parameter(nb_variable+1,0);
     value O; // offset vector
     value values;
     std::vector<double> probability;
     std::vector<unsigned int> sup_bound;
     iterator it;

     assert(pdist != NULL);

     if((frequency != NULL) && (marginals != NULL)){
         mean_computation();
         means = get_means();
         if(min_inf_bound_flag){
             mdist = new stat_tool::DiscreteParametric((int)((min_offset + max_alloc_nb_value + 1) * stat_tool::SAMPLE_NB_VALUE_COEFF), stat_tool::POISSON);
             (*marginals)[0]->poisson_estimation(mdist, min_inf_bound, min_inf_bound_flag, cumul_threshold);
             minf_bound = mdist->inf_bound;
             delete mdist;
             for(var = 1; var < nb_variable; var++){
                 mdist = new stat_tool::DiscreteParametric((int)((min_offset + max_alloc_nb_value + 1) * stat_tool::SAMPLE_NB_VALUE_COEFF), stat_tool::POISSON);
                 (*marginals)[var]->poisson_estimation(mdist, min_inf_bound, min_inf_bound_flag, cumul_threshold);
                 minf_bound = MIN(minf_bound, MAX(mdist->inf_bound,0));
                 delete mdist;
             }
             mdist = NULL;
         } else {
           minf_bound = min_inf_bound;
         }

         for(var = 0; var < nb_variable; var++){
             s += means[var]-minf_bound;
         }
         nb_val = get_nb_element();

         O.resize(nb_variable, minf_bound);
         if(get_frequency(O) != 0){
             parameter[0] = 1/((double)nb_variable - 1) * (s + log(((double)get_frequency(O))/((double)nb_val))); // estimation of the first parameter
         } else {
             std::vector< std::vector<double> > VarCovar;
             VarCovar = get_covariances();
             for(var = 0; var < nb_variable-1; var++){
                 for(var2 = var+1; var2 < nb_variable; var2++){
                     parameter[0] += 2/((double)nb_variable * (double)(nb_variable-1)) * VarCovar[var][var2];
                 }
             }
         }
         if (parameter[0] < 0.)
            parameter[0] = 0;

         for(var = 0; var < nb_variable; var++){
             parameter[var+1] = means[var]-minf_bound-parameter[0]; // esitmation of all parameters
             if (parameter[var+1] <= 0.)
                status = false;
         }

         if (status)
            pdist->init(nb_variable, MPOISSON, minf_bound, sup_bound, parameter, probability, cumul_threshold);
         else
            return stat_tool::D_INF;
     }
     likelihood = likelihood_computation(*pdist);

     return likelihood;
}

template<typename Type> DiscreteMultivariateReestimation<Type>*
DiscreteMultivariateReestimation<Type>::get_bootstrap_distribution(unsigned int nb_ind, double cumul_threshold) const
{
    unsigned int i, j, tot;
    unsigned int dim = nb_value_computation();
    std::vector<double> parameter(1, nb_ind);
    std::vector<unsigned int> sup_bound;
    std::vector<double> probability(dim);
    iterator it;
    i = 0;
    for(it = begin(); it != end(); it++){
        probability[i] = (*it).second;
        i++;
    }
    DiscreteMultivariateParametric *sample_dist = new DiscreteMultivariateParametric(dim, MMULTINOMIAL, 0, sup_bound, parameter, probability);
    std::vector<unsigned int> sample = sample_dist->simulation();
    delete sample_dist;
    offset_computation();
    max_computation();
    DiscreteMultivariateReestimation<Type> *breest = new DiscreteMultivariateReestimation<int>(nb_variable, get_offset(), max_alloc_nb_value-get_offset()+1);
    i = 0;
    for(it = begin(); it != end(); it++){
        if(sample[i] > 0)
          breest->update((*decode)[(*it).first], sample[i]);
        i++;
    }

    return breest;
}

# ifdef DEBUG

/*****************************************************************
 *
 *  Compute term of the log likelihood to optimize with respect to parameter
 *  in negative multinomial setting from DiscreteMultivariateReestimation
 *  using data, inferior bound, counts and parameter value (solve likelihood equations)
 *
 **/

template <typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_parameter_estimation_likelihood_target(const DiscreteMultivariateParametric& dist,
                                                                                                    int inf_bound,
                                                                                                    const std::vector<Type>& mean_counts,
                                                                                                    Type total_counts,
                                                                                                    double parameter) const
{
   unsigned int i, sum_val; // sum of values
   // likelihood  to be maximized
   double target = 0., val_target; // sum associated to value n
   value val;
   Type freq;
   iterator it;

   assert((frequency != NULL) && (parameter > 0));

   for(it = begin(); it != end(); it++)
   {
      val_target = .0;
      freq = (*it).second;
      val = (*decode)[(*it).first];
      sum_val = 0;
      for(i = 0; i < nb_variable; i++)
         sum_val += (val[i] - inf_bound);
      for(i = 0; i < sum_val; i++)
         val_target += log(parameter + i);
      // cout << "Lgamma rel. error:" << std::abs((val_target - (lgamma(parameter + sum_val) - lgamma(parameter))) / val_target) << endl;
      target += val_target * freq;
   }
   target -= parameter * nb_element * log(1 + (total_counts / (nb_element * parameter)));

   for(i = 0; i < nb_variable; i++)
      target += mean_counts[i] * log(mean_counts[i] /(parameter * nb_element + total_counts));

   return target;
}

/*****************************************************************
 *
 *  Compute likelihood curve with respect to parameter for probabilities
 *  satisfying ML equations in negative multinomial setting from DiscreteMultivariateReestimation
 *
 **/

template <typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_parameter_partial_likelihood(const std::vector<Type>& mean_counts,
                                                                                          int iinf_bound,
                                                                                          Type total_counts,
                                                                                          double parameter) const
{
   unsigned int var;
   double res;
   std::vector<unsigned int> sup_bound(nb_variable);
   std::vector<double> probability(nb_variable), vparameter(1);
   DiscreteMultivariateParametric *dist = NULL;

   vparameter[0] = parameter;
   for(var = 0; var < nb_variable; var++)
      probability[var] = 1 - (mean_counts[var] / ((parameter * nb_element) + total_counts));
   dist = new DiscreteMultivariateParametric(nb_variable, MNEGATIVE_MULTINOMIAL,
                                             iinf_bound, sup_bound, vparameter, probability);
   res = likelihood_computation(*dist);
   delete dist;
   dist = NULL;
   return res;
}

/*****************************************************************
 *
 *  Numerical computation of likelihood derivative  with respect to parameter for probabilities
 *  satisfying ML equations in negative multinomial setting from DiscreteMultivariateReestimation
 *
 **/

template <typename Type> double
DiscreteMultivariateReestimation<Type>::negative_multinomial_parameter_partial_likelihood_derivative(const std::vector<Type>& mean_counts,
                                                                                                     int iinf_bound,
                                                                                                     Type total_counts,
                                                                                                     double parameter,
                                                                                                     double step) const
{
   double l1, l2;

   assert(step != 0);

   l1 = negative_multinomial_parameter_partial_likelihood(mean_counts, iinf_bound, total_counts, parameter);
   l2 = negative_multinomial_parameter_partial_likelihood(mean_counts, iinf_bound, total_counts, parameter+step);

   return (l2 - l1) / step;
}
# endif

/*****************************************************************
 *
 *  Destructor for DiscreteMultivariateReestimation class
 *
 **/

template<typename Type> void
DiscreteMultivariateReestimation<Type>::remove()
{
   unsigned int v;

   if (frequency != NULL)
   {
      delete frequency;
      frequency = NULL;
   }

   if (code != NULL)
   {
      delete code;
      code = NULL;
   }

   if (decode != NULL)
   {
      delete decode;
      decode = NULL;
   }

   if (marginals != NULL) // *marginals is a std::vector
   {
      for(v = 0; v < marginals->size(); v++)
         if ((*marginals)[v] != NULL)
         {
            delete (*marginals)[v];
            (*marginals)[v] = NULL;
         }

      delete marginals;
      marginals = NULL;
   }
   if (sums != NULL)
   {
      delete sums; // *sums is a std::vector
      sums = NULL;
   }
}

/*****************************************************************
 *
 *  Copy operator for DiscreteMultivariateReestimation class
 *  Does not perform deallocation
 *
 **/

template<typename Type> void
DiscreteMultivariateReestimation<Type>::copy(const DiscreteMultivariateReestimation& dist)
{
   unsigned int v;

   nb_variable = dist.nb_variable;
   min_offset = dist.min_offset;
   max_alloc_nb_value = dist.max_alloc_nb_value;
   nb_element = dist.nb_element;
   if (dist.frequency != NULL)
      frequency = new values_dict(*dist.frequency);
   else
      frequency = NULL;
   if (dist.code != NULL)
      code = new values_code(*dist.code);
   else
      code = NULL;
   if (dist.decode != NULL)
      decode = new values_decode(*dist.decode);
   else
      decode = NULL;

   if (dist.marginals != NULL)
   {
      marginals = new std::vector< stat_tool::Reestimation<Type> *>(nb_variable);
      for(v = 0; v < nb_variable; v++)
      {
         if ((*dist.marginals)[v] != NULL)
            (*marginals)[v] = new stat_tool::Reestimation<Type>(*(*dist.marginals)[v]);
         else
            (*marginals)[v] = NULL;
      }
   }
   else
      marginals = NULL;

   if (dist.sums != NULL)
   {
      sums = new std::vector<Type>(dist.sums->size());
      for(v = 0; v < sums->size(); v++)
         (*sums)[v] = (*dist.sums)[v];
   }
   else
      sums = NULL;

}

/*****************************************************************
 *
 *  Code value of DiscreteMultivariateReestimation
 *  Return a pair with first element equals true iif element
 *  is already present in self, and second element is the code
 *  (or pow(max_alloc_nb_value, nb_variable) if the code is incorrect)
 *
 **/

template<typename Type> std::pair<bool, unsigned int>
DiscreteMultivariateReestimation<Type>::decoder(const value& v) const
{
   bool success = true;
   unsigned int ind, tcode, val;
   std::pair<bool, unsigned int> res;

     if (frequency != NULL)
   {
      assert(v.size() == nb_variable);
      // compute code in dictionary
      tcode = 0;
      for(ind = 0; ind < nb_variable; ind++)
      {
         val = v[ind]-min_offset;
         if (val >= max_alloc_nb_value)
            success = false;
               tcode += (unsigned int)(val * pow(max_alloc_nb_value, ind));
      }
      if (success)
      {
         res.second = tcode;
         // check whether element (i.e. tcode)
         // is contained in dictionary
         res.first = (decode->count(tcode) > 0);
            } else {
                res.second = pow(max_alloc_nb_value, nb_variable);
                res.first = false;
            }
   }
   else
   {
      res.first = false;
      res.second = pow(max_alloc_nb_value, nb_variable);
   }
   return res;
}

/*****************************************************************
 *
 *  Get an iterator on DiscreteMultivariateReestimation
 *
 **/

template<typename Type> typename
Stat_trees::DiscreteMultivariateReestimation<Type>::iterator
// typename std::map<unsigned int, Type >::iterator
DiscreteMultivariateReestimation<Type>::end() const
{
   assert(frequency != NULL);
   return frequency->end();
}

/*****************************************************************
 *
 *  Get an iterator on DiscreteMultivariateReestimation
 *
 **/

template<typename Type> typename
Stat_trees::DiscreteMultivariateReestimation<Type>::iterator
// typename std::map<unsigned int, Type >::iterator
DiscreteMultivariateReestimation<Type>::begin() const
{
   assert(frequency != NULL);
   return frequency->begin();
}

#endif
