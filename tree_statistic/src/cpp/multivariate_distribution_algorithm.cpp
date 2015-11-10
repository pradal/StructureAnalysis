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
 *       $Id: multivariate_distribution_algorithm.cpp 9554 2010-09-20 17:10:30Z jbdurand $
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
#include "statiskit/core/data/marginal/multivariate.h"

#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"   // definition of DiscreteParametricModel class
#include "stat_tool/markovian.h"   // definition of DiscreteParametricProcess class
#include "stat_tool/stat_label.h"   // definition of DiscreteParametricProcess class
#include "stat_tool/vectors.h"

#include "tree_labels.h"

#include <assert.h>

#include <boost/math/special_functions/gamma.hpp>

#include "mv_histogram_tools.h"   // definition of DiscreteParametricModel class

#include "multivariate_distribution.h"   // definition of DiscreteParametricModel class

#include <algorithm>

extern char* label(const char*);

using namespace stat_tool;
using namespace Stat_trees;
using namespace boost::math;

/*****************************************************************
 *
 *  Compute log likelihood for a model
 *
 **/

double DiscreteMultivariateDistribution::likelihood_computation(const Statiskit::Marginal::Multivariate::CountTable& hist) const
{
   double log_likelihood = 0, mass;
   const std::vector<std::vector<unsigned int> > values =
       Stat_histogram_get_elements(hist);

   unsigned int i;

#  ifdef DEBUG
   std::cout << "like.." << std::endl;
#  endif
   for(i=0; i < values.size(); i++)
   {
#        ifdef DEBUG
      std::cout << values[i] << std::endl;
#        endif
      mass = this->get_mass(values[i], true); // compute log mass
#        ifdef DEBUG
      cout << "D_INF: " << D_INF << endl;
      cout << "mass > D_INF: " << (mass > D_INF) << endl;
#        endif
      if (mass > D_INF)
      {
         log_likelihood += mass * Stat_histogram_get_frequency(hist, values[i]);
      }
      else
      {
         return D_INF;
      }
   }
   return log_likelihood;
}

/*****************************************************************
 *
 *  Compute distribution entropy
 *
 **/

double DiscreteMultivariateDistribution::information_computation() const
{
   unsigned int i;
   double information = 0, cumul = 0;

   assert(mass != NULL);
   assert(mass->size() == total_nb_value);

   for(i = 0; i < total_nb_value; i++)
      if ((*mass)[i] > 0.)
         information += (*mass)[i] * log((*mass)[i]);

   if (complement > 0.)
      information -=  (1. - complement) * log(1. - complement);

   information /= (1. - complement);

   return information;
}

/*****************************************************************
 *
 *  Simulation of one single vector
 *
 **/

DiscreteMultivariateDistribution::index_array DiscreteMultivariateDistribution::simulation() const
{
   index_array res(0);

   return res;
}

/*****************************************************************
 *
 *  Simulation of several vectors
 *
 **/

Vectors* DiscreteMultivariateDistribution::simulation(unsigned int nb_value) const
{
   Vectors *res = NULL;
   return res;
}

/*****************************************************************
 *
 *  Compute probability for every possible index and set mass accordingly
 *  Does nothing if mass already is the right size
 *
 **/

void DiscreteMultivariateParametric::mass_computation() const
{
   const unsigned int eff_nb_value = min(total_nb_value, MULTIVARIATE_DISTRIBUTION_MAX_VALUE+1);
   unsigned int i;
   index_array indices;

   assert(mass != NULL);

   if (mass->size() != eff_nb_value)
   {
      mass->resize(eff_nb_value, D_DEFAULT);
      if (total_nb_value < MULTIVARIATE_DISTRIBUTION_MAX_VALUE)
      {
         // compute every probabiliy value
         indices = loop_init();
         set_mass(*counters, get_mass(*counters));

         for(i = 0; i < total_nb_value-1; i++)
         {
            set_mass(indices, get_mass(indices));
            indices = loop_next();
         }
      }
      else
      {
         // compute first values only
         indices = loop_init(offset);
         set_mass(*counters, get_mass(*counters));

         for(i = 0; i < MULTIVARIATE_DISTRIBUTION_MAX_VALUE; i++)
         {
            set_mass(indices, get_mass(indices));
            indices = loop_next_low(offset);
         }
         set_mass(indices, D_DEFAULT);
      }
   }

}

/*****************************************************************
 *
 *  Simulation of one single vector
 *
 **/
// EDIT
DiscreteMultivariateDistribution::index_array DiscreteMultivariateParametric::simulation() const
{
     unsigned int i = 0;
     unsigned int s = 0;
     DiscreteParametric *dist = NULL;
     index_array res(0);

     switch(ident){
         case MMULTINOMIAL : {
             dist = extract_marginal(i);
             res.push_back(dist->simulation());
             delete dist;
             s += res[i];
             while(i < nb_variable-2){
                 i++;
                 dist = conditional_distribution(i, res);
                 res.push_back(dist->simulation());
                 s += res[i];
                 delete dist;
             }
             dist = NULL;
             res.push_back(parameter[0]-s);
             break;

         }
         case MNEGATIVE_MULTINOMIAL : {
             dist = extract_marginal(i);
             res.push_back(dist->simulation());
             delete dist;
             while(i < nb_variable-1){
                 i++;
                 dist = conditional_distribution(i, res);
                 res.push_back(dist->simulation());
                 delete dist;
             }

             dist = NULL;
             break;
         }
         case MPOISSON : {
             dist = new DiscreteParametric(POISSON, 0, I_DEFAULT, parameter[i], D_DEFAULT);
             unsigned int res0 = dist->simulation();
             delete dist;
             while(i < nb_variable){
                 i++;
                 dist = new DiscreteParametric(POISSON, inf_bound, I_DEFAULT, parameter[i], D_DEFAULT);
                 res.push_back(dist->simulation()+res0);
                 delete dist;
             }
             dist = NULL;
             break;
         }
         default : {
             break;
         }
     }

     return res;
}

/*****************************************************************
 *
 *  Simulation of several vectors
 *
 **/

Vectors* DiscreteMultivariateParametric::simulation(unsigned int nb_value) const
{
   unsigned int i, j;
   int *iidentifier = NULL;
   int **iint_vector = NULL;
   index_array s;
   Vectors *res = NULL;

   assert(nb_value > 0);

   iint_vector = new int*[nb_value];
   iidentifier = new int[nb_value];

   for(i = 0; i < nb_value; i++)
   {
      iint_vector[i] = new int[nb_variable];
      iidentifier[i] = i;
      s = simulation();
      for(j = 0; j < nb_variable; j++)
         iint_vector[i][j] = s[j];
   }

   res = new Vectors(nb_value, iidentifier, nb_variable, iint_vector);

   for(i = 0; i < nb_value; i++)
      delete [] iint_vector[i];
   delete [] iint_vector;
   delete [] iidentifier;

   return res;
}

/*****************************************************************
 *
 *  Compute probability values for negative multinomial distributions
 *
 **/

double DiscreteMultivariateParametric::negative_multinomial_get_mass(const index_array& indices,
                                                                     bool log_computation) const
{
   bool in_support = true;
   unsigned int i;
   index_array dindices = indices; // indices minus offset

   double log_res = parameter[0] * log((*rparameter)[0]),
          s = 0; // sum of dindices

   for(i = 0; i < nb_variable; i++)
   {
      if (dindices[i] < inf_bound)
         in_support = false; // index not in support
      dindices[i] -= inf_bound;
   }

   if (in_support)
   {
      for(i = 0; i < nb_variable; i++)
         s += dindices[i];

      for(i = 0; i < nb_variable; i++)
         if (probability[i] == 1.)
         {
            if (dindices[i] == 0)
               log_res -= lgamma(dindices[i]+1);
            else
            {
               log_res = D_INF;
               break;
            }
         }
         else
            log_res += dindices[i] * log((1-probability[i])) - lgamma(dindices[i]+1);

      if (log_res != D_INF)
         log_res -= log(tgamma_delta_ratio(parameter[0], s));

   }
   else
      log_res = D_INF;
   if (log_computation)
      return log_res;
   else
      if (log_res > D_INF)
         return exp(log_res);
      else
         return .0;

}

/*********************************************************************
 *
 *  Compute probability values for multivariate Poisson distributions
 *
 **/
// EDIT
double DiscreteMultivariateParametric::mpoisson_get_mass(const index_array& indices, bool log_computation) const
{
    unsigned int i,imin;
  double logmass = 0;
  index_array devent = indices;
  unsigned int min = std::numeric_limits<unsigned int>::max();
  for(i = 0; i < nb_variable; i++) {
        if(devent[i]-inf_bound < 0) { // Not in support
            if(log_computation)
                return D_INF;
            else
                return 0.0;
        } else { // Looking for the minimum
            devent[i] -= inf_bound;
            if(devent[i] < min) {
                min = devent[i];
                imin = i;
            }
        }
  }
  if(min == 0) {
        index_array O(nb_variable,0);
        if(std::equal(devent.begin(), devent.end(), O.begin())) { // mass = e^(-sum_i theta_i)
            for(i = 0; i < nb_variable+1; i++)
                logmass -= parameter[i];
            if(log_computation)
                return logmass;
            else
                return std::exp(logmass);
        } else { // P[X=x] = P[X=x-1] prod_{X_i > 0} theta_i/x_i
            for(i = 0; i < nb_variable; i++) {
                if(devent[i] == 0)
                    devent[i] += inf_bound;
                else {
                    logmass += devent[i]*std::log(parameter[i+1]) - boost::math::lgamma(devent[i]+1);
                    devent[i] += inf_bound-1;
                }
            }
            if(log_computation)
                return logmass+get_mass(devent,true);
            else
                return std::exp(logmass)*get_mass(devent,false);
        }
  } else { // x_j P[X=x] = theta_0*P[X=x-1] + theta_j*P[X_{-j}=x_{-j}, X_j = x_j-1]
        index_array event0=devent, eventi=devent;
        for(i = 0; i < nb_variable; i++) {
            event0[i] += inf_bound-1;
            eventi[i] += inf_bound;
            if(i == imin)
                eventi[i]--;
        }
        if(log_computation)
            return std::log((parameter[0]*get_mass(event0,false)+parameter[imin+1]*get_mass(eventi,false))/(eventi[imin]+1));
        else
            return parameter[0]*get_mass(event0,false)+parameter[imin+1]*get_mass(eventi,false)/(eventi[imin]+1);
  }
/*    bool in_support = true, min_found=false;
    unsigned int i, j, min;
    double s, param_min, log_res = 0;
    index_array dindices = indices; // indices minus offset
    index_array oindices(indices.size(),0); // increasing order of indices
    i = 0;
    while(i<nb_variable and in_support){
        if(dindices[i] < inf_bound)
            in_support = false;
            dindices[i] -= inf_bound;
            i++;
    }
    if(in_support){
        index_array O(nb_variable, 0);
        bool equal = false;
        equal = std::equal(dindices.begin(), dindices.end(), O.begin());
        bool loop = false;
        if (!equal){
            min = dindices[0];
            i = 0;
            bool loop = (min != 0);
            while(loop & i < nb_variable){
                min = MIN(min,dindices[i]);
                if(min == 0){
                    loop = false;
                } else {
                    i++;
                }
            }
            s = 0;
                            if(min == 0){ // Use the first recurrence P[X=x]= P[X=x-1]*prod(X_i != 0) theta_i/x_i
                            for(j = 0;j < nb_variable; j++){
                                if(dindices[j] != 0){
                                    log_res += dindices[j]*std::log(parameter[j+1])-lgamma(dindices[j]+1);
                                }
                            }
            } else {
                index_array ddindices = dindices;
                for(i = 0; i < nb_variable; i++){
                    dindices[i]--;
                }
                i = 0;
                while(i < nb_variable && !min_found){
                    if(ddindices[i]==min){
                        min_found = true;
                        ddindices[i]--;
                        param_min = parameter[i+1];
                    }
                    i++;
                }
                for(i = 0; i < nb_variable; i++){
                    ddindices[i] += inf_bound;
                    dindices[i] += inf_bound;
                }
                log_res = log(param_min*get_mass(ddindices)+parameter[0]*get_mass(dindices))-std::log(min);
            }
        } else {
                    for(i = 0; i < nb_variable+1; i++)
            log_res -= parameter[i];
                }
    } else {
        log_res = D_INF;
    }
    if (log_computation){
        return log_res;
    } else {
        if (log_res > D_INF)
            return std::exp(log_res);
        else
            return .0;
    }*/
}

/*********************************************************************
 *
 *  Compute probability values for multinomial distributions
 *
 **/

double DiscreteMultivariateParametric::mmultinomial_get_mass(const index_array& indices, bool log_computation) const
{
   bool in_support = true;
   double log_res = 0;
   unsigned int i, s = 0;
   index_array dindices = indices;
   assert(dindices.size() == nb_variable);

   for(i = 0; i < nb_variable; i++)
   {
       if(dindices[i] < inf_bound)
          in_support = false;
       dindices[i] -= inf_bound;
       s += dindices[i];
   }

   if (s > (parameter[0]-nb_variable*inf_bound)) // Check value of sum
       in_support = false;

   if (in_support)
   {
      log_res = lgamma(s+1);
      for(i = 0; i < nb_variable; i++)
         log_res += log(probability[i])*dindices[i]-lgamma(dindices[i]+1);
   }
   else
      log_res = D_INF;
   if (log_computation)
       return log_res;
   else
       if(log_res > D_INF)
           return std::exp(log_res);
       else
           return 0.;

}

/*****************************************************************
 *
 *  Compute conditional distribution of variable i given variables
 *  0,...,i-1 for negative multinomial distributions
 *
 **/

DiscreteParametric* DiscreteMultivariateParametric::DiscreteMultivariateParametric::negative_multinomial_conditional(unsigned int ivariablei,
                                                                                                                     const index_array& values) const
{
   unsigned int i;
   double param = parameter[0], proba = (*rparameter)[0];
   DiscreteParametric *res = NULL;

   assert((ivariablei > 0) && (ivariablei < nb_variable) && (values.size() == ivariablei));

   for(i = 0; i < ivariablei; i++)
   {
      param += (values[i] - inf_bound);
      proba += (1-probability[i]);
   }


   proba = proba / ((1. - probability[ivariablei]) + proba);

   res = new DiscreteParametric(NEGATIVE_BINOMIAL,
                                inf_bound, I_DEFAULT,
                                param, proba);
   return res;
}

/*
 * Compute conditional distribution of variable i given variables
 * 0, ..., i-1 for multinomial distributions
 * */
// EDIT
DiscreteParametric* DiscreteMultivariateParametric::mmultinomial_conditional(unsigned int ivariablei, const index_array& values) const
{
     unsigned int i;
     unsigned int param = sup_bound[ivariablei];
     double proba = 1;
     DiscreteParametric *res = NULL;

     assert((ivariablei > 0) && (ivariablei < nb_variable) && (values.size() == ivariablei) && (ivariablei < nb_variable));

     for(i = 0; i < ivariablei; i++){
         param -= (values[i]-inf_bound);
         proba -= probability[i];
     }

     proba = probability[ivariablei]/proba;

     res = new DiscreteParametric(BINOMIAL,
                                                                 inf_bound, param,
                                                                 D_DEFAULT, proba);

     return res;
}

/*
 * Compute the distribution of offspring sum for multivariate poisson dsitribution
 * */

DiscreteParametric* DiscreteMultivariateParametric::mpoisson_sum() const
{
    unsigned int i, j, var;
    double s;
    DiscreteParametric *sum = NULL;
    sum = new DiscreteParametric();
    sum->ident = 0;
    sum->inf_bound = nb_variable*inf_bound;
    sum->nb_value = (int)round(inf_bound + (mean[i] - inf_bound + sqrt(variance[i])) * 20.);
    if(sum->nb_value == sum->inf_bound)
        sum->nb_value++;
    sum->offset = sum->inf_bound;
    sum->nb_parameter = 2;
    sum->mass = new double[sum->nb_value];
    sum->cumul = new double[sum->nb_value];
    sum->mean = 0;
    sum->variance = 0;
    double sp = 0;
    for(var = 0; var < nb_variable; var++){
        sp += parameter[var+1];
    }

    for(i = 0; i < sum->nb_value; i++){
        sum->mass[i] = 0;
        if(i > sum->offset){
            sum->mass[i] = 1;
            s = 0;
            for(var = 0; var < nb_variable+1; var++){
                sum->mass[i] *= std::exp(-parameter[var]);
            }
            for(j = 0; j <= round(i/((double)nb_variable)); j++){
                s += pow(sp, i - nb_variable*j)/boost::math::tgamma(i - nb_variable*j +1) * pow(parameter[0], j)/boost::math::tgamma(j+1);
            }
            sum->mass[i] *= s;
            sum->mean += sum->mass[i]*i;
          sum->variance += sum->mass[i]*pow(i,2) - sum->mass[i]*i;
        }
    }
    sum->cumul[0] = sum->mass[0];
    sum->max = sum->mass[0];
    for(i = 1; i < sum->nb_value; i++){
        sum->cumul[i] = sum->cumul[i-1]+sum->mass[i];
        sum->max = MAX(sum->max, sum->mass[i]);
    }

    return sum;
}


/*
 * Compute the directional cdf for multivariate poisson
 * */

DiscreteParametric* DiscreteMultivariateParametric::mpoisson_dcdf() const
{
    unsigned int i, j, k, var, min;
    double s, p;
    DiscreteParametric *dcdf = NULL;
    dcdf = new DiscreteParametric();
    dcdf->ident = 0;
    dcdf->inf_bound = inf_bound;
    dcdf->nb_value = 0;
    for(var = 0; var < nb_variable; var++){
        dcdf->nb_value = MAX(dcdf->nb_value, (int)round(inf_bound + (mean[i] - inf_bound + sqrt(variance[i])) * 20.));
    }
    if(dcdf->nb_value == inf_bound)
        dcdf->nb_value++;
    dcdf->offset = dcdf->inf_bound;
    dcdf->nb_parameter = nb_variable + 1;
    dcdf->cumul = new double[dcdf->nb_value];

    for(i = 0; i < dcdf->nb_value; i++){
        dcdf->cumul[i] = 0;
        if(i > dcdf->offset){
            dcdf->cumul[i] = 0;
            for(j = 0; j <= i; j++){
                p = pow(parameter[0],j)/boost::math::tgamma(j+1);
                for(var = 0; var < nb_variable; var++){
                    s = 0;
                    for(k = 0; k < i - j; k++){
                        s += pow(parameter[var+1], k)/boost::math::tgamma(k+1);
                    }
                    p *= s;
                }
                dcdf->cumul[i] += p;
            }
            for(var = 0; var < nb_variable+1; var++){
                dcdf->cumul[i] *= std::exp(-parameter[var]);
            }
        }
    }

    return dcdf;
}

/*****************************************************************
 *
 *  Compute distribution entropy in IidDiscreteMultivariateParametric
 *
 **/

double IidDiscreteMultivariateParametric::information_computation() const
{
   if (marginal != NULL)
      return (nb_variable * marginal->information_computation());
}


/*****************************************************************
 *
 *  Return the probability associated with a set of indices
 *
 **/

double IidDiscreteMultivariateParametric::get_mass(const index_array& indices) const
{
   double res = 1;
   unsigned int i;

   assert(marginal != NULL);
   assert((indices.size() == nb_variable) && (total_nb_value > 0));

   for(i = 0; i < nb_variable; i++)
   {
      assert(indices[i] < nb_value[i]);
      res *= marginal->mass[indices[i]];
   }
   return res;
}

/*****************************************************************
 *
 *  Simulation
 *
 **/

DiscreteMultivariateDistribution::index_array IidDiscreteMultivariateParametric::simulation() const
{
   unsigned int i;

   index_array res(nb_variable);

   for(i = 0; i < nb_variable; i++)
      res[i] = marginal->simulation();

   return res;
}

/*****************************************************************
 *
 *  Compute conditional distribution of variable i given variables
 *  0,...,i-1 for negative multinomial distributions
 *
 **/

DiscreteParametric* IidDiscreteMultivariateParametric::IidDiscreteMultivariateParametric::conditional_distribution(unsigned int ivariablei,
                                                                               const index_array& values) const
{
   DiscreteParametric *res = new DiscreteParametric(*marginal);

   return res;
}

/*****************************************************************
 *
 *  Return the probability associated with a set of indices
 *
 **/

double MultinomialCompoundDiscreteParametric::get_mass(const index_array& indices) const
{

  unsigned int tot = 0, i;
  double res = 0;

  assert((indices.size() == nb_variable));

  for(i = 0; i < nb_variable; i++){
        tot += indices[i]-inf_bound;
        res += log(probability[i])*(indices[i]-inf_bound)-lgamma(indices[i]-inf_bound+1);
    }
    res += lgamma(tot+1);
    res = exp(res);

    switch(param_compound->ident){
        case BINOMIAL : {
            res *= pdf(binomial(param_compound->sup_bound-param_compound->inf_bound, param_compound->probability), tot-param_compound->inf_bound);
            break;
        }
        default :
          break;
    }

    return res;
}

/*****************************************************************
 *
 *  Simulation
 *
 **/

DiscreteMultivariateDistribution::index_array MultinomialCompoundDiscreteParametric::simulation() const
{
    unsigned int i, ni = 0;
    double p = 1;
    DiscreteParametric *dist = NULL;


    ni = param_compound->simulation(); // Simulate the sum

  index_array res(nb_variable,0);
    for(i = 0; i < nb_variable-1; i++){ // multinomial simulation
        dist = new DiscreteParametric(BINOMIAL, 0, ni, D_DEFAULT, probability[i]/p, cumul_threshold);
        res[i] = dist->simulation();
        delete dist;
        ni -= res[i];
        p -= probability[i];
        res[i] += inf_bound;
    }
    res[nb_variable-1] = ni+inf_bound;
  dist = NULL;

  return res;
}

std::vector<double> MultinomialCompoundDiscreteParametric::parametric_mean_computation() const
{
    unsigned int i;
    std::vector<double> res;
    res.resize(nb_variable);

    for(i = 0; i < nb_variable; i++){
            res[i] = (param_compound->mean-param_compound->inf_bound) * probability[i] + inf_bound;
    }

    return res;
}

std::vector<double> MultinomialCompoundDiscreteParametric::parametric_variance_computation() const
{
    unsigned int i;
    std::vector<double> res;
    res.resize(nb_variable);

    for(i = 0; i < nb_variable; i++){
        res[i] = probability[i] * ((1-probability[i])*(param_compound->mean-param_compound->inf_bound) + probability[i] * param_compound->variance);
    }

    return res;
}

double  MultinomialCompoundDiscreteParametric::parametric_covariance_computation(unsigned int ivariable1, unsigned int ivariable2) const
{
    double res;
    assert((ivariable1 < nb_variable) && (ivariable2 < nb_variable));

    if(ivariable1 == ivariable2){
        std::vector<double> variance;
        variance = parametric_variance_computation();
        res = variance[ivariable1];
    } else {
        res = - probability[ivariable1] * probability[ivariable2] * ((param_compound->mean - param_compound->inf_bound) - param_compound->variance);
    }

    return res;
}

DiscreteParametric* MultinomialCompoundDiscreteParametric::marginals_compound_multinomial(unsigned int ivariablei) const
{
    DiscreteParametric *marginali = NULL;
    double s, max;
    unsigned int i,j;

    switch(param_compound->ident){
        case BINOMIAL : {
            marginali = new DiscreteParametric();
            marginali->ident = 0;
            marginali->inf_bound = inf_bound;
            marginali->sup_bound = sup_bound[ivariablei];
            marginali->probability = probability[ivariablei];
            marginali->nb_value = sup_bound[ivariablei]+1;
            marginali->offset = inf_bound;
            marginali->max = 0;
            marginali->nb_parameter = 3;
            marginali->mass = new double[marginali->nb_value];
            marginali->cumul = new double[marginali->nb_value];
            marginali->mean = mean[ivariablei];
            marginali->variance = variance[ivariablei];
            for(i = 0; i < marginali->nb_value; i++){
                if(i < inf_bound){
                    marginali->mass[i] = 0;
                } else {
                    marginali->mass[i] = 0;
                    for(j = i-inf_bound; j <= sup_bound[ivariablei]; j++){
                        s = boost::math::pdf(boost::math::binomial(j, probability[ivariablei]), i-inf_bound);
                        s *= boost::math::pdf(boost::math::binomial(sup_bound[ivariablei], param_compound->probability), j);
                        marginali->mass[i] += s;
                    }
                }
            }
            marginali->cumul[0] = 0;
            max = 0;
            for(i = 1; i < marginali->nb_value; i++){
                max = MAX(max,marginali->mass[i]);
                marginali->cumul[i] = marginali->cumul[i-1]+marginali->mass[i];
            }
            marginali->max = max;
            break;
        }
        default :
            assert(false);
            break;
    }

    return marginali;
}

/*
 * Compute the directional cdf for coumpound multinomial distribution
 * */

/*DiscreteParametric* MultinomialCompoundDiscreteParametric::compound_multinomial_dcdf() const
{
    unsigned int i, j, k, var, min;
    double s, p;
    DiscreteParametric *dcdf = NULL;
    dcdf = new DiscreteParametric();
    dcdf->ident = 0;
    dcdf->inf_bound = inf_bound;
    dcdf->nb_value = (int)round(param_compound->inf_bound + (param_compound->mean - param_compound->inf_bound + sqrt(param_compound->variance)) * 20.);
    if(dcdf->nb_value == inf_bound)
        dcdf->nb_value++;
    dcdf->offset = (int)round(param_compound->inf_bound/nb_variable);
    dcdf->nb_parameter = nb_variable;
    switch(param_compound->ident){
        case BINOMIAL :
            dcdf->nb_parameter += 2;
            break;
        default :
            break;
    }
    dcdf->cumul = new double[dcdf->nb_value];

    for(i = 0; i < dcdf->nb_value; i++){
        dcdf->cumul[i] = 0;
        if(i > dcdf->offset){
            dcdf->cumul[i] = 0;
            for(j = 0; j <= i; j++){
                p = pow(parameter[0],j)/boost::math::tgamma(j+1);
                for(var = 0; var < nb_variable; var++){
                    s = 0;
                    for(k = 0; k < i - j; k++){
                        s += pow(parameter[var+1], k)/boost::math::tgamma(k+1);
                    }
                    p *= s;
                }
                dcdf->cumul[i] += p;
            }
            for(var = 0; var < nb_variable+1; var++){
                dcdf->cumul[i] *= std::exp(-parameter[var]);
            }
        }
    }

    return dcdf;

}*/
