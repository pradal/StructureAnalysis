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

#ifndef MULTIVARIATE_DISTRIBUTION_H
#define MULTIVARIATE_DISTRIBUTION_H

#include "stat_tool/stat_tools.h"
#include <valarray>
#include <map>
#include "tool/rw_tokenizer.h"

namespace Stat_trees
{

/****************************************************************
 *
 *  Purpose :
 *     provides classes to implement discrete multivariate
 *     distributions
 */

/*! \file multivariate_distribution.h
    \brief Purpose:
     provide the class "DiscreteMultivariateDistribution"
     and derived classes (DiscreteMultivariateParametric,
     IidDiscreteMultivariateParametric,
     MultinomialCompoundMultivariateParametric)
*/

/*************************************************************
 *
 *  Constants:
 */

/// maximal number of computed probability values
/// for multivariate distributions
const unsigned int MULTIVARIATE_DISTRIBUTION_MAX_VALUE = 100000;
const unsigned int MULTIVARIATE_DISTRIBUTION_NB_ITER = 500;
const double MULTIVARIATE_DISTRIBUTION_LIKELIHOOD_DIFF = 1.e-6;

/// types of multivariate distributions
enum {
  DISCRETE_MULTIVARIATE,
  MIID,
  MPOISSON,
  MNEGATIVE_MULTINOMIAL,
  MMULTINOMIAL, // MULTINOMIAL already exists in stat_tool.h
  MCOMPOUND_MULTINOMIAL
};

/*************************************************************
 *
 *  Template methods:
 */

template <class T> std::ostream& operator<<(std::ostream &os, const std::vector<T>& o)
{
   const unsigned int l = o.size();
   unsigned int i;

   for(i = 0; i < l; i++)
      os << o[i] << "   ";
   os << " ; ";

   return os;
}


/****************************************************************
 *
 *  Class definitions:
 */

template <typename Type> class DiscreteMultivariateReestimation;

/**
   \class DiscreteMultivariateDistribution
   \brief An implementation of discrete multivariate distributions
*/

class DiscreteMultivariateDistribution
{  /// Discrete multivariate distributions

public :

   // public types
   typedef std::vector<unsigned int> index_array;

protected :

   /// partial products of the number of values
   /// for the remaining dimensions
   /// (for use with valarrays)
   index_array *val_strides;
   /// a counter to iterate on indices
   index_array *counters;

   /// compute the val_strides
   void val_strides_computation() const;
   void remove();
   void copy(const DiscreteMultivariateDistribution& dist);

public :

   /// number of variables
   unsigned int nb_variable;
   /// number of values per variable (starting from 0, not from offset)
   std::vector<unsigned int> nb_value;
   /// total number of values
   unsigned int total_nb_value;
   // if needed ?
   // std::vector<unsigned int> alloc_nb_value;
   /// for each variable, 1st value with probability > 0
   std::vector<unsigned int> offset;
   // if needed ?
   // std::vector<unsigned double> max;
   // ou
   // double max;
   /// complementary probability (>0 for improper distributions)
   double complement;
   /// means
   std::vector<double> mean;
   /// variances
   std::vector<double> variance;
   /// probability of each value
   /// the probabilities for the first value of every variable
   /// is stored first, then the values of the last variable
   /// increase, and so on
   std::valarray<double> *mass;
   /// marginals
   std::vector<stat_tool::Distribution*> marginals;
   /// number of unknown parameters
   unsigned int nb_parameter;

   /** Initialize an empty DiscreteMultivariateDistribution */
   DiscreteMultivariateDistribution();
   /** Copy constructor */
   DiscreteMultivariateDistribution(const DiscreteMultivariateDistribution& dist);

   /** Copy a DiscreteMultivariateDistribution and return a pointer */
   virtual DiscreteMultivariateDistribution* copy_ptr() const;

   /** Delete a DiscreteMultivariateDistribution */
   virtual ~DiscreteMultivariateDistribution();

   /** Get the probability value for a set of indices */
   virtual double get_mass(const index_array& indices) const;

   /** Get either the probability value for a set of indices,
       or the log probability */
   virtual double get_mass(const index_array& indices, bool log_computation) const;

   /** Set the probability value for a set of indices */
   bool set_mass(const index_array& indices, double value) const;

   /** Get marginal distribution (a new instance is allocated) */
   virtual stat_tool::Distribution* get_marginal(unsigned int ivariable) const;

   /** Init a loop on every possible index */
   index_array loop_init() const;

   /** Init a loop on every possible index, starting at a given index
       (lexicographic order)*/
   index_array loop_init(const index_array& begin_index) const;

   /** Loop on every possible index (lexicographic order) */
   index_array loop_next() const;

   /** Init a loop on every possible index (lowest indices first) */
   // index_array loop_init_low() const;

   /** Init a loop on every possible index, starting at a given index
       (lowest indices first) */
   // index_array loop_init_low(const index_array& begin_index) const;

   /** Loop on every possible index (lowest indices first) */
   index_array loop_next_low() const;

   /** Loop on every possible index (lowest indices first),
       respecting some given minimal values */
   index_array loop_next_low(const index_array& min_index) const;

   /** Likelihood computation */
   virtual double likelihood_computation(const Statiskit::Marginal::Multivariate::CountTable& hist) const;

   /** Entropy computation */
   virtual double information_computation() const;

   /** Simulation of one vector */
   virtual index_array simulation() const;

   /** Simulation of several vectors */
   virtual stat_tool::Vectors* simulation(unsigned int nb_value) const;

   /** Print distribution with full specification */
   virtual std::ostream& ascii_write(std::ostream &os, bool exhaustive = false) const;

   /** Print the nature of distribution and number of variables */
   std::ostream& header_ascii_print(std::ostream &os) const;

   /** Print header and probabilities */
   virtual std::ostream& ascii_print(std::ostream &os, bool comment_flag,
                                     bool exhaustive, bool file_flag, bool cumul_flag,
                                     const Statiskit::Marginal::Multivariate::CountTable *histo = NULL) const;

   /** Print probabilities */
   virtual std::ostream& probability_ascii_print(std::ostream &os, bool comment_flag = false,
                                                 bool cumul_flag = false) const;

   /** Spreadsheet-like output of probabilities */
   virtual std::ostream& spreadsheet_print(std::ostream &os) const;

   /** Matplotlib output of the distribution */
   virtual stat_tool::MultiPlotSet* get_plotable(const Statiskit::Marginal::Multivariate::CountTable* histo = NULL) const;

   /** Permutation of the dimensions in \e self*/
   void variable_permutation(stat_tool::StatError& error, int* perm, bool check_flag = true);

};

DiscreteMultivariateDistribution* discrete_multivariate_ascii_read(stat_tool::StatError &error , const char *path ,
                                                                   double cumul_threshold = stat_tool::CUMUL_THRESHOLD);

DiscreteMultivariateDistribution* discrete_multivariate_parsing(stat_tool::StatError &error, std::ifstream &in_file ,
                                                                int &line, int last_ident = MCOMPOUND_MULTINOMIAL,
                                                                double cumul_threshold = stat_tool::CUMUL_THRESHOLD,
                                                                int min_inf_bound = 0);

DiscreteMultivariateDistribution* discrete_multivariate_probability_parsing(stat_tool::StatError &error, std::ifstream &in_file ,
                                                                            int &line, unsigned int inb_variable,
                                                                            double cumul_threshold = stat_tool::CUMUL_THRESHOLD,
                                                                            int min_inf_bound = 0);

/**
   \class DiscreteMultivariateParametric
   \brief An implementation of parametric discrete multivariate distributions
   \details Particular choices of parameterization:
   \li In negative multinomial distributions, the sum of the complement to one
       of each probability must be between 0 and 1.
   \li In multinomial distributions, the parameter corresponds
       to the sum of vector components.
*/

class DiscreteMultivariateParametric : public stat_tool::StatInterface, public DiscreteMultivariateDistribution
{  // Parametric multivariate discrete distributions
   template<typename DummyType> friend class DiscreteMultivariateReestimation;

protected :

   /// identifier of the distribution type
   int ident;
   /// inf bound (common to each variable)
   unsigned int inf_bound;
   /// sup bound (specific to each variable)
   std::vector<unsigned int> sup_bound;
   /// parameters (Poisson,...)
   std::vector<double> parameter;
   /// probabilities (binomial, ...)
   std::vector<double> probability;
   /// marginals (parametric)
   std::vector<stat_tool::DiscreteParametric*> param_marginals;

   /// reparametrization used for computational purpose
   std::vector<double>* rparameter;

   void copy(const DiscreteMultivariateParametric& dist);
   void remove();

   virtual std::ostream& spreadsheet_write(std::ostream& os) const;
   virtual bool plot_write(const char * prefix, const char * title) const;

   /** Set distribution type and parameters */
   void init(int inb_variable, int iident, unsigned int iinf_bound,
             const std::vector<unsigned int>& isup_bound,
             const std::vector<double>& vparameter,
             const std::vector<double>& vprobability,
             double cumul_threshold = stat_tool::CUMUL_THRESHOLD);

   /** Compute probability values for negative multinomial distributions */
   double negative_multinomial_get_mass(const index_array& indices,
                                        bool log_computation=false) const;

   /** Compute probability values for multinomial */
   double mmultinomial_get_mass(const index_array& indices, bool log_computation=false) const;

   /** Compute probability values for multivariate poisson distributions */
   double mpoisson_get_mass(const index_array& indices, bool log_computation=false) const;

   /** Conditional distribution of variables
      for negative multinomial distributions */
   stat_tool::DiscreteParametric* negative_multinomial_conditional(unsigned int ivariablei,
                                                        const index_array& values) const;

   /** Conditional distribution of variables for multinomial distributions */
   stat_tool::DiscreteParametric* mmultinomial_conditional(unsigned int ivariablei, const index_array& values) const;

   /** Sum distribution for multivariate poisson */
   stat_tool::DiscreteParametric* mpoisson_sum() const;

   /** Directional cdf for multivariate poisson */
   stat_tool::DiscreteParametric* mpoisson_dcdf() const;

   /** Joint distribution for multivariate poisson */
   DiscreteMultivariateParametric* mpoisson_joint(unsigned int ivariable1, unsigned int ivariable2) const;

public :

   /** Initialize an empty DiscreteMultivariateParametric */
   DiscreteMultivariateParametric();
   /** Initialize a DiscreteMultivariateParametric */
   DiscreteMultivariateParametric(int inb_variable, int iident, unsigned int iinf_bound,
                                  const std::vector<unsigned int>& isup_bound,
                                  const std::vector<double>& vparameter,
                                  const std::vector<double>& vprobability,
                                  double cumul_threshold = stat_tool::CUMUL_THRESHOLD);

   /** Copy constructor of DiscreteMultivariateParametric */
   DiscreteMultivariateParametric(const DiscreteMultivariateParametric& dist);

   /** Copy a DiscreteMultivariateParametric and return a pointer */
   virtual DiscreteMultivariateParametric* copy_ptr() const;

   /** Delete a DiscreteMultivariateParametric */
   virtual ~DiscreteMultivariateParametric();

   /** Compute the number of free parameters */
   virtual unsigned int nb_parameter_computation() const;

   /** Update the number of free parameters */
   void nb_parameter_update();

   /** Get inferior bound of distribution */
   unsigned int get_inf_bound() const;

   /** Get superior bounds of distribution */
   std::vector<unsigned int> get_sup_bound() const;

   /** Get distribution parameters */
   std::vector<double> get_parameter() const;

   /** Get distribution probabilities */
   std::vector<double> get_probability() const;

   /** Get the probability value for a set of indices */
   double get_mass(const index_array& indices) const;

   /** Get either the probability value for a set of indices,
       or the log probability */
   virtual double get_mass(const index_array& indices, bool log_computation) const;

   /** Set every probability values */
   void mass_computation() const;

   /** Simulation of one vector*/
   virtual index_array simulation() const;

   /** Simulation of several vectors */
   virtual stat_tool::Vectors* simulation(unsigned int nb_value) const;

   /** Analytic computation of expectation */
   virtual std::vector<double> parametric_mean_computation() const;
   /** Analytic computation of variance */
   virtual std::vector<double> parametric_variance_computation() const;
   /** Analytic computation of covariance between two variables */
   virtual double parametric_covariance_computation(unsigned int ivariable1, unsigned int ivariable2) const;

   /** Conditional distribution of variable \e i given variables 0, ..., \e i-1 */
   virtual stat_tool::DiscreteParametric* conditional_distribution(unsigned int ivariablei,
                                                        const index_array& values) const;

   /** Compute parametric marginal distributions */
   virtual void compute_param_marginals(double cumul_threshold = stat_tool::CUMUL_THRESHOLD);

   /** Extend the support of a marginal distribution */
   virtual void extend_param_marginal(unsigned int ivariable, int value_update) const;

   /** Get marginal distribution (a new instance is allocated) */
   virtual stat_tool::Distribution* get_marginal(unsigned int ivariable) const;

   /** Return a parametric marginal distribution (a new instance is allocated) */
   stat_tool::DiscreteParametric* extract_marginal(unsigned int ivariable) const;

   /** Print distribution on a single line */
   virtual std::ostream& line_write(std::ostream &os) const;

   /** Print distribution with full specification */
   // method inherited from StatInterface, pure virtual if not redefined
   virtual std::ostream& ascii_write(std::ostream &os, bool exhaustive = false) const;

   /** Print distribution with full specification */
   virtual std::ostream& ascii_write(std::ostream &os, bool exhaustive, bool file_flag) const;
   /** Print distribution with full specification into a file */
   virtual bool ascii_write(stat_tool::StatError &error, const char *path, bool exhaustive = false) const;
   /** Spreadsheet-like output of the distribution into a file */
   virtual bool spreadsheet_write(stat_tool::StatError &error, const char *path) const;
   /** Gnuplot output of the distribution */
   virtual bool plot_write(stat_tool::StatError &error, const char *prefix, const char *title = NULL) const;
   /** Matplotlib output of the distribution */
   virtual stat_tool::MultiPlotSet* get_plotable(const Statiskit::Marginal::Multivariate::CountTable *histo = NULL) const;

   /** Print distribution parameters */
   virtual std::ostream& ascii_print(std::ostream &os, bool comment_flag,
                                     bool exhaustive, bool file_flag, bool cumul_flag,
                                     const Statiskit::Marginal::Multivariate::CountTable *histo = NULL) const;

   /** Spreadsheet-like output for parameters */
   virtual std::ostream& spreadsheet_print(std::ostream &os) const;

   /** Print distribution characteristics and marginals */
   virtual std::ostream& ascii_characteristic_print(std::ostream &os, bool shape = false,
                                                    bool comment_flag = false) const;

   /** Spreadsheet-like output for characteristics and marginals */
   virtual std::ostream& spreadsheet_characteristic_print(std::ostream &os, bool shape = false) const;

   /** Permutation of the dimensions in \e self*/
   void variable_permutation(stat_tool::StatError& error, int* perm, bool check_flag = false);
};

/** Parse a DiscreteMultivariateParametric from a file */

DiscreteMultivariateParametric*
discrete_multivariate_parametric_parsing(stat_tool::StatError &error,
                                         std::ifstream &in_file,
                                         int &line,
                                         RWCTokenizer &next,
                                         int iident,
                                         unsigned int inb_variable,
                                         double cumul_threshold = stat_tool::CUMUL_THRESHOLD,
                                         int min_inf_bound = 0);

class IidDiscreteMultivariateParametric : public DiscreteMultivariateParametric
{  // Iid Parametric discrete multivariate distributions

protected :

   /// identifier for discrete univariate parametric distribution
   unsigned int uident;
   /// common sup bound
   unsigned int sup_bound1;
   /// common parameter
   double parameter1;
   /// common probability
   double probability1;
   // marginal
   stat_tool::DiscreteParametric* marginal;

   void remove();
   void copy(const IidDiscreteMultivariateParametric& dist);

public :

   /** Initialize a DiscreteMultivariateParametric */
   IidDiscreteMultivariateParametric(int iident, int inb_variable, unsigned int iinf_bound,
                                     unsigned int isup_bound,
                                     double dparameter,
                                     double dprobability,
                                     double cumul_threshold = stat_tool::CUMUL_THRESHOLD);
   /** Copy constructor */
   IidDiscreteMultivariateParametric(const IidDiscreteMultivariateParametric& dist);

   /** Copy an IidDiscreteMultivariateParametric and return a pointer */
   virtual IidDiscreteMultivariateParametric* copy_ptr() const;

   /** Delete a DiscreteMultivariateParametric */
   virtual ~IidDiscreteMultivariateParametric();

   /** Compute the number of free parameters */
   unsigned int nb_parameter_computation() const;

   /** Print distribution on a single line */
   std::ostream& line_write(std::ostream &os) const;

   // /** Print distribution with full specification */
   // std::ostream& ascii_write(std::ostream &os, bool exhaustive = false) const;


   /** Print distribution parameters */
   std::ostream& ascii_print(std::ostream &os, bool comment_flag,
                             bool exhaustive, bool file_flag, bool cumul_flag,
                             const Statiskit::Marginal::Multivariate::CountTable *histo = NULL) const;

  /** Entropy computation */
   double information_computation() const;

   /** Get the probability value for a set of indices */
   double get_mass(const index_array& indices) const;

   /** Simulation of one vector */
   virtual index_array simulation() const;


   /** Conditional distribution of variable \e i given variables 0, ..., \e i-1 */
   stat_tool::DiscreteParametric* conditional_distribution(unsigned int ivariablei,
                                                const index_array& values) const;

};

IidDiscreteMultivariateParametric* iid_discrete_multivariate_parsing(stat_tool::StatError &error,
                                                                     std::ifstream &in_file,
                                                                     int &line,
                                                                     unsigned int inb_variable,
                                                                     double cumul_threshold = stat_tool::CUMUL_THRESHOLD,
                                                                     int min_inf_bound = 0);

class MultinomialCompoundDiscreteParametric : public DiscreteMultivariateParametric
{
    // Multinomial distribution with random number of replications

protected :

   /// compound distribution
   stat_tool::DiscreteParametric *param_compound;
   double cumul_threshold;
   /// Destructor
   void remove();
   /// Copy operator
   void copy(const MultinomialCompoundDiscreteParametric& dist);

public :

   /** Initialize an empty MultinomialCompoundDiscreteParametric */
   MultinomialCompoundDiscreteParametric();

   /** Initialize a MultinomialCompoundDiscreteParametric */
   MultinomialCompoundDiscreteParametric(int inb_variable,
                                         const std::vector<double>& iprobability,
                                         const stat_tool::DiscreteParametric& dist,
                                         double icumul_threshold = stat_tool::CUMUL_THRESHOLD);

   /** Copy constructor */
   MultinomialCompoundDiscreteParametric(const MultinomialCompoundDiscreteParametric& dist);

   /** Copy a MultinomialCompoundDiscreteParametric and return a pointer */
   virtual MultinomialCompoundDiscreteParametric* copy_ptr() const;

   /** Delete a MultinomialCompoundDiscreteParametric*/
   virtual ~MultinomialCompoundDiscreteParametric();

    void init(int inb_variable,
                                        const std::vector<double>& iprobability,
                                        const stat_tool::DiscreteParametric& dist,
                                        double icumul_threshold = stat_tool::CUMUL_THRESHOLD);

   /** Compute the number of free parameters */
   unsigned int nb_parameter_computation() const;

   /** Get either the probability value for a set of indices, or the log probability */
   double get_mass(const index_array& indices) const;

   /** Simulation of one vector */
   virtual index_array simulation() const;

   /** Analytic computation of expectation */
   std::vector<double> parametric_mean_computation() const;
   /** Analytic computation of variance */
   std::vector<double> parametric_variance_computation() const;
   /** Analytic computation of covariance between two variables */
   double parametric_covariance_computation(unsigned int ivariable1, unsigned int ivariable2) const;

   /** Print distribution parameters */
   virtual std::ostream& ascii_print(std::ostream &os, bool comment_flag,
                                     bool exhaustive, bool file_flag, bool cumul_flag,
                                     const Statiskit::Marginal::Multivariate::CountTable *histo = NULL) const;

    /** Compute parametric marginal distributions */
    void compute_param_marginals(double cumul_threshold = stat_tool::CUMUL_THRESHOLD);

  /** Compute the ivariabliei-th parametric marginal distribution */
    stat_tool::DiscreteParametric* marginals_compound_multinomial(unsigned int ivariablei) const;

};

/** Parse a DiscreteMultivariateParametric from a file */

MultinomialCompoundDiscreteParametric*
multinomial_compound_parsing(stat_tool::StatError &error,
                             std::ifstream &in_file,
                             int &line,
                             RWCTokenizer &next,
                             unsigned int inb_variable,
                             double cumul_threshold = stat_tool::CUMUL_THRESHOLD,
                             int min_inf_bound = 0);

/**
   \class DiscreteMultivariateReestimation
   \brief Datasets associated with discrete multivariate distributions
   and empirical distribution with integral or floating frequencies
   (EM estimation)
*/


template <typename Type>
class DiscreteMultivariateReestimation
{  /// Datasets associated with Discrete multivariate distributions
   /// and empirical distribution with integral or floating frequencies
   /// (EM estimation)

public :

   // public types
   /// values of observation vectors
   typedef std::vector<unsigned int> value;
   /// number of observation vectors for each value
   /// key is the value of the vector coded in base max_nb_value
   typedef std::map<unsigned int, Type > values_dict;
   /// mapping between values and code
   typedef std::map<value, unsigned int > values_code;
   /// mapping between code and values
   typedef std::map<unsigned int, value > values_decode;
   /// iterator on codes
   typedef typename values_dict::iterator iterator;

protected :
   /// number of variables
   unsigned int nb_variable;
   /// minimum value of offset
   unsigned int min_offset;
   /// maximum number of allocated values (maximum over each variable)
   unsigned int max_alloc_nb_value;
   /// sum of every element in dictionary (equivalent to sample size)
   Type nb_element;
   /// number of observation vectors for each value
   /// key is the value of the vector (minus offset) coded
   /// in base max_nb_value
   values_dict *frequency;
   /// mapping between values (offset not substracted) and code
   values_code *code;
   /// mapping between code (offset not substracted) and values
   values_decode *decode;
   /// marginal frequencies and means
   std::vector<stat_tool::Reestimation<Type> *> *marginals;
   /// marginal sums (frequencies * values)
   std::vector<Type> *sums;
   /// Destructor
   void remove();
   /// Copy operator
   void copy(const DiscreteMultivariateReestimation& dist);

   /// return the code of v in code if v is present
   /// return I_DEFAULT otherwise
   std::pair<bool, unsigned int> decoder(const value& v) const;

   /** Iterator on codes */
   iterator begin() const;

   /** Iterator on codes */
   iterator end() const;

public :

   /** Initialize an empty DiscreteMultivariateDistribution */
   DiscreteMultivariateReestimation();

   /** Initialize an empty DiscreteMultivariateDistribution
       from the number of variables, the smallest possible value,
       and the difference between largest and smallest possible values */
   DiscreteMultivariateReestimation(unsigned int inb_variable, unsigned int imin_offset,
                                    unsigned int imax_nb_value);

   /** Initialize a DiscreteMultivariateDistribution from stat_tool::Vectors*/
   DiscreteMultivariateReestimation(const stat_tool::Vectors& vec);

   /** Change properties of DiscreteMultivariateDistribution */
   void init(unsigned int inb_variable, unsigned int imin_offset,
             unsigned int imax_nb_value);

   /** Copy constructor */
   DiscreteMultivariateReestimation(const DiscreteMultivariateReestimation<Type>& dist);

   /** Delete a DiscreteMultivariateReestimation */
   ~DiscreteMultivariateReestimation();

   /** Get number of variables */
   unsigned int get_nb_variable() const;

   /** Get elements present in \e self */
   std::vector<value> get_elements() const;

   /** Get marginal (return a pointer; object should not be deallocated) */
   const stat_tool::Reestimation< Type >* get_marginal_ptr(unsigned int ivariable) const;

   /** Get marginal (a new instance is allocated) */
   stat_tool::Reestimation< Type >* get_marginal(unsigned int ivariable) const;

   /** Get marginal FrequencyDistribution (a new instance is allocated) */
   stat_tool::FrequencyDistribution* get_marginal_frequency(unsigned int ivariable) const;

   /** Get sums */
   std::vector<Type> get_sums() const;

   /** Get offset value */
   unsigned int get_offset() const;

   /** Get maximal number of values */
   unsigned int get_max_alloc_nb_value() const;

   /** Return number of elements with different keys in dictionary
       and compute the number of values in marginals */
   unsigned int nb_value_computation() const;

   /** Return the sum of every values in dictionary */
   Type get_nb_element() const;

   /** Compute and return the sum of every values in dictionary */
   Type nb_element_computation() const;

   /** Return offset value and compute offsets in marginals */
   unsigned int offset_computation() const;

   /** Return maximal value and compute maxima in marginals */
   unsigned int max_computation() const;

   /** Compute and return sums */
   std::vector<Type> sum_computation() const;

   /** Compute means */
   void mean_computation() const;

   /** Compute variances */
   void variance_computation(bool biais = false) const;

   /** Get means */
   std::vector<double> get_means() const;

   /** Get variances */
   std::vector<double> get_variances() const;

   /** Update frequency for a given value */
   void update(const value& v, Type freq);

   /** Get frequency for a given value */
   Type get_frequency(const value& v) const;

   /** Print frequencies and characteristics */
   std::ostream& print(std::ostream& os) const;

   /** Matplotlib output of histograms */
   stat_tool::MultiPlotSet* get_plotable() const;

   /** Log likelihood computation */
   double likelihood_computation(const DiscreteMultivariateDistribution& dist) const;

      /** Log likelihood computation */
   double likelihood_computation(const MultinomialCompoundDiscreteParametric& dist) const;

   /** Compute FrequencyDistribution nb_element, mean, etc.
       for marginals */
   void frequency_distribution_computation() const;

   /** Compute histogram of the sum of elements in each vector */
   stat_tool::FrequencyDistribution* get_sum_frequency_distribution() const;

   /** Estimation of a DiscreteMultivariateDistribution */
   DiscreteMultivariateDistribution* distribution_estimation() const;

   /** Estimation of a DiscreteMultivariateParametric with given type */
   DiscreteMultivariateParametric* parametric_estimation(stat_tool::StatError &error,
                                                         int ident, double &likelihood,
                                                         int min_inf_bound = 0,
                                                         bool min_inf_bound_flag = true,
                                                         double cumul_threshold = stat_tool::CUMUL_THRESHOLD) const;

   /** Estimation of both a DiscreteMultivariateParametric and its type */
   DiscreteMultivariateParametric* type_parametric_estimation(stat_tool::StatError &error,
                                                              double &likelihood,
                                                              int min_inf_bound = 0,
                                                              bool min_inf_bound_flag = true,
                                                              double cumul_threshold = stat_tool::CUMUL_THRESHOLD) const;

   /** Estimation of negative multinomial distributions using marginals only */
   double negative_multinomial_marginal_estimation(DiscreteMultivariateParametric *pdist, int min_inf_bound,
                                                   bool min_inf_bound_flag, double cumul_threshold) const;

   /** Maximum likelihood estimation of negative multinomial distributions */
   double negative_multinomial_estimation(DiscreteMultivariateParametric *pdist, int min_inf_bound,
                                          bool min_inf_bound_flag, double cumul_threshold) const;

   /** Optimization of the log likelihood with respect to parameter
       in negative multinomial setting.
       Should rather use GSL optim, LOPTI or equivalent */

   double negative_multinomial_parameter_estimation(const DiscreteMultivariateParametric& pdist,
                                                    int iinf_bound,
                                                    const std::vector<Type>& mean_counts,
                                                    Type total_counts,
                                                    double initial_parameter,
                                                    double tol = MULTIVARIATE_DISTRIBUTION_LIKELIHOOD_DIFF,
                                                    unsigned int nb_iter = MULTIVARIATE_DISTRIBUTION_NB_ITER) const;

   /** Optimization of the log likelihood with respect to parameter
       in negative multinomial setting.
       Bounds for parameter value are specified */

   std::pair<double, double> negative_multinomial_parameter_estimation(const DiscreteMultivariateParametric& pdist,
                                                                       int iinf_bound,
                                                                       const std::vector<Type>& mean_counts,
                                                                       Type total_counts,
                                                                       double parameter_sup,
                                                                       double parameter_inf,
                                                                       double tol,
                                                                       unsigned int nb_iter) const;

   /** Term of the score log likelihood to maximize with respect to parameter
       in negative multinomial setting. */

   double negative_multinomial_parameter_estimation_score_target(const DiscreteMultivariateParametric& pdist,
                                                                 int iinf_bound,
                                                                 const std::vector<Type>& mean_counts,
                                                                 Type total_counts,
                                                                 double parameter) const;

   /** Term of the score log likelihood to cancel with respect to parameter
       in negative multinomial setting. */

   double negative_multinomial_parameter_estimation_score_value(const DiscreteMultivariateParametric& pdist,
                                                                int iinf_bound,
                                                                const std::vector<Type>& mean_counts,
                                                                Type total_counts,
                                                                double parameter) const;

   /** Maximum likelihood  estimation of compound multinomial distributions */
   double compound_multinomial_estimation(MultinomialCompoundDiscreteParametric *pdist,
                                          int iident, int min_inf_bound, bool min_inf_bound_flag,
                                          double tol) const;

   /** Moments estimation method for multivariate Poisson distributions */
   double mpoisson_estimation_M(DiscreteMultivariateParametric *pdist, int min_inf_bound, bool
                                min_inf_bound_flag, double cumul_threshold) const;

   /** Maximum likelihood estimation of multivariate Poisson distributions */
   double mpoisson_estimation_ML(DiscreteMultivariateParametric* pdist,
                                 int min_inf_bound, bool min_inf_bound_flag,
                                 double tol) const;

   /** Maximum likelihood estimation of multinomial distributions */
   double mmultinomial_estimation(DiscreteMultivariateParametric *pdist, int min_inf_bound,
                                  bool min_inf_bound_flag, double tol) const;

   /** Get covariances */
   std::vector< std::vector<double> > get_covariances() const;

   /** Permutation of the dimensions in \e self*/
   void variable_permutation(stat_tool::StatError& error, int* perm);

     /** Get a bootstrap distribution from a DiscreteNultivariateReestimation */
     DiscreteMultivariateReestimation<Type>* get_bootstrap_distribution(unsigned int nb_ind, double cumul_threshold = stat_tool::CUMUL_THRESHOLD) const;

   # ifdef DEBUG
   /** Term of the log likelihood to maximize with respect to parameter
       in negative multinomial setting. */

   double negative_multinomial_parameter_estimation_likelihood_target(const DiscreteMultivariateParametric& pdist,
                                                                      int iinf_bound,
                                                                      const std::vector<Type>& mean_counts,
                                                                      Type total_counts,
                                                                      double parameter) const;

   /** Compute likelihood curve with respect to parameter for probabilities
       satisfying ML equations. */
   double negative_multinomial_parameter_partial_likelihood(const std::vector<Type>& mean_counts,
                                                            int iinf_bound,
                                                            Type total_counts,
                                                            double parameter) const;

   /** Numerical computation of likelihood derivative  with respect to parameter for probabilities
       satisfying ML equations.*/

   double negative_multinomial_parameter_partial_likelihood_derivative(const std::vector<Type>& mean_counts,
                                                                       int iinf_bound,
                                                                       Type total_counts,
                                                                       double parameter,
                                                                       double step) const;

   # endif
};

#include "multivariate_distribution.hpp"

};


#endif
