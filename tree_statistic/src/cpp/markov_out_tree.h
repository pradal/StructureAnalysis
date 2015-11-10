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

#ifndef MARKOV_OUT_TREE_H
#define MARKOV_OUT_TREE_H


namespace Stat_trees
{

/*! \file markov_out_tree.h
    \brief Purpose:
     provide the class "MarkovOutTree", which corresponds
     to Markov out trees with dependent children states,
     given their parent, and generation distributions,
     under different hypotheses on order of children
     (including multitype branching processes)
*/


/****************************************************************
 *
 *  Class definitions :
 */

class MarkovOutTreeData;
template <typename Type> class MarkovOutTreeReestimation;

/**
   \class GenerationProcess
   \brief A family of DiscreteMultivariateDistributions, mainly used
          for the representation of generation distributions in
          multitype branching processes
*/

class GenerationProcess
{  // a family of DiscreteMultivariateDistributions, mainly used
   // for the representation of generation distributions in
   // multitype branching processes

   friend class MarkovOutTree;

   friend GenerationProcess* generation_process_parsing(StatError &error,
                                                        std::ifstream &in_file,
                                                        int &line,
                                                        int nb_state,
                                                        double cumul_threshold);

public :
   // public types

   /// configuration of unordered children states
   /// or factors
   /// Configurations of factors are in the following order:
   /// Parent state, ordered states, other factors
   typedef std::vector<unsigned int> index_array;

   typedef Statiskit::Marginal::Multivariate::CountTable mvHistogram;

protected :

   /// number of states or types
   int nb_state;
   /// number of ordered children involved in the branching processes
   /// as factors (due to dependencies regarding these states, excluding parent state)
   unsigned int nb_children_branching;
   /// number of other factors involved in the branching process
   unsigned int nb_factors;
   /// number of possible values for each factor
   index_array *factor_values;
   /// number of generation processes
   /// (combining all possible values of every factor)
   unsigned int nb_generation_values;
   /// number of factors for generation processes
   /// (including parent state)
   unsigned int nb_generation;
   /// base to code and decode configuration of factors
   index_array factor_base;
   /// generation distributions: generation[state]
   DiscreteMultivariateDistribution **generation;

   /** Copy operator */
   void copy(const GenerationProcess& process);

   /** Destructor */
   void remove();

   std::ostream& ascii_print(std::ostream &os,
                             mvHistogram **empirical_observation,
                             bool exhaustive, bool file_flag) const;

   std::ostream& spreadsheet_print(std::ostream &os, int process,
                                   mvHistogram **empirical_observation = NULL) const;

   /** Gnuplot output of \e self*/
   bool plot_print(const char * prefix, const char * title, int process,
                   mvHistogram **empirical_observation = NULL) const;

public :

   /** Create a GenerationProcess using the number of (ordered) children states,
       of factors and parent states determining the total number inb_generation
       of branching processes, and the generation distributions */
   GenerationProcess(int inb_state = 0,
                     unsigned int inb_children_branching = 0,
                     unsigned int inb_factors = 0, const index_array *ifactor_values = NULL,
                     unsigned int inb_generation = 0,
                     DiscreteMultivariateDistribution * const * pgeneration = NULL);
   /** Copy constructor of GenerationProcess class */
   GenerationProcess(const GenerationProcess& process);
   /** Destructor for GenerationProcess class */
   ~GenerationProcess();
   /** Assignement operator for GenerationProcess class */
   GenerationProcess& operator=(const GenerationProcess& process);

   /** Code a configuration of factors with an unsigned int */
   unsigned int code(const index_array& fact) const;
   /** Decode a configuration of factors coded with an unsigned int */
   index_array decode(unsigned int code) const;

   /** get number of ordered children involved in the branching processes
       as factors (due to dependencies regarding these states,
       excluding parent state) */
   unsigned int get_nb_children_branching() const;

   /** get number of other factors involved in the branching process */
   unsigned int get_nb_factors() const;

   /** get number of possible values for each factor */
   index_array get_factor_values() const;

   /* get number of generation processes
      (combining all possible values of every factor) */
   unsigned int get_nb_generation_values() const;

   /** get number of factors for generation processes
       (including parent state) */
   unsigned int get_nb_generation() const;

   /** Matplotlib output of generation process */
   MultiPlotSet* get_plotable(Statiskit::Marginal::Multivariate::CountTable* const * histo = NULL) const;

   /** Return the number of views (i.e. the size) in Matplotlib output */
   unsigned int nb_plot_set_computation(Statiskit::Marginal::Multivariate::CountTable* const * histo = NULL) const;

   /** Write views into a MultiPlotSet object */
   void plotable_write(MultiPlotSet &plot_set, int &index, GenerationProcess::mvHistogram* const * histo = NULL) const;

   /** Permutation of the states of \e self*/
   void state_permutation(StatError& error, int* perm) const;
};

/** Read a generation process from a file */
GenerationProcess* generation_process_parsing(StatError &error,
                                              std::ifstream &in_file,
                                              int &line,
                                              int nb_state,
                                              double cumul_threshold = CUMUL_THRESHOLD);
/**
   \class MarkovOutTree
   \brief a Markov out-tree model, with dependent
    children states, given their parent

*/

class MarkovOutTree : public StatInterface, protected Chain
{  // a Markov out-tree model, with dependent
   // children states, given their parent

   friend class MarkovOutTreeData;

   friend MarkovOutTree* markov_out_tree_parsing(StatError& error,
                                                 std::ifstream &in_file,
                                                 int &line);

   friend MarkovOutTree* markov_out_tree_ascii_read(StatError& error,
                                                    const char * path,
                                                    int depth,
                                                    double cumul_threshold);

   friend ostream& operator<<(ostream &os, const MarkovOutTree& markov);

public :
   // public types

   // configuration of unordered children states
   // or factors
   // typedef std::vector<unsigned int> index_array;
   typedef double** double_array_2d;
   typedef double*** double_array_3d;
   typedef double**** double_array_4d;

protected :

   /// data associated with the model
   MarkovOutTreeData *markov_data;
   /// number of output processes
   unsigned int nb_output_process;

   // An XProcess is a family of distributions indexed
   // by states. XProcess ** represent families of XProcess*,
   // one for each variable.

   /// output processes modeled as (discrete) nonparametric
   CategoricalTreeProcess **nonparametric_process;
   /// output processes modeled as discrete parametric
   DiscreteParametricProcess **discrete_parametric_process;
   /// output processes modeled as continuous (parametric)
   ContinuousParametricProcess **continuous_parametric_process;
   /// rank (respective with partial order of children)
   /// beyond which children are unordered
   /// (including child at rank unordered_rank)
   unsigned int unordered_rank;
   /// number of variable order Markov chain models on ordered children
   unsigned int nb_vomc;
   /// variable order Markov chain models on ordered children
   sequence_analysis::VariableOrderMarkov **vomc_models;
   /// generation processes
   GenerationProcess *generation_distributions;

   /** Copy operator */
   void copy(const MarkovOutTree& markov, bool data_flag = true);
   /** Destructor */
   void remove();

   /** Copy operator for state process */
   void state_copy(const MarkovOutTree& markov);

   /** Initialization of output processes */
   void output_process_init(const int* nb_value);

   /** Copy operator for output process */
   void observation_copy(int inb_output_process,
                         CategoricalProcess * const * npprocess,
                         DiscreteParametricProcess * const * dpprocess,
                         ContinuousParametricProcess * const * cpprocess);

   /** generic function for ascii output, into a file or on screen,
      using potential empirical data at hand */
   ostream& ascii_write(ostream& os, const MarkovOutTreeData* otrees,
                        bool exhaustive = false, bool file_flag = false,
                        bool hidden = false) const;

   ostream& spreadsheet_write(ostream& os, const MarkovOutTreeData * tree,
                              const Test * test = NULL) const;

   bool plot_write(const char * prefix, const char * title,
                   const MarkovOutTreeData * otrees) const;

   /** Matplotlib output of the MarkovOutTree model using data*/
   MultiPlotSet* get_plotable(const MarkovOutTreeData * otrees) const;

   int nb_parameter_computation(double min_probability = 0.) const;

   double* memory_computation() const;

   void state_no_occurrence_probability(int state, double increment = sequence_analysis::LEAVE_INCREMENT);
   void state_first_occurrence_root_distribution(int state, int min_nb_value = 1,
                                                 double cumul_threshold = CUMUL_THRESHOLD);
   void state_first_occurrence_leaves_distribution(int state, int min_nb_value = 1,
                                                   double cumul_threshold = CUMUL_THRESHOLD);
   void state_leave_probability(const double * memory, int state,
                                double increment = sequence_analysis::LEAVE_INCREMENT);
   void state_sojourn_size_distribution(const double * memory, int state,
                                        int min_nb_value = 1,
                                        double cumul_threshold = CUMUL_THRESHOLD);
   void state_nb_pattern_mixture(int state, char pattern);
   void output_no_occurrence_probability(int variable, int output,
                                         double increment = sequence_analysis::LEAVE_INCREMENT);
   void output_first_occurrence_root_distribution(int variable, int output,
                                                  int min_nb_value = 1,
                                                  double cumul_threshold = CUMUL_THRESHOLD);
   void output_first_occurrence_leaves_distribution(int variable, int output,
                                                    int min_nb_value = 1,
                                                    double cumul_threshold = CUMUL_THRESHOLD);
   void output_leave_probability(const double * memory,
                                 int variable, int output,
                                 double increment = sequence_analysis::LEAVE_INCREMENT);
   void output_sojourn_size_distribution(const double * memory, int variable,
                                         int output, int min_nb_value = 1 ,
                                         double cumul_threshold = CUMUL_THRESHOLD);
   void output_nb_zones_mixture(int variable, int output);
   void output_nb_occurrences_mixture(int variable, int output);

   void state_marginal_distribution(const MarkovOutTreeData& trees,
                                    double_array_3d& state_marginal,
                                    int index = I_DEFAULT) const;

   double** state_marginal_distribution(const Trees& trees,
                                        int index) const;

   // methods related to likelihood computation

   double likelihood_correction(const MarkovOutTreeData& trees) const;

   /** Compute output conditional distributions */
   void output_conditional_distribution(const MarkovOutTreeData& trees,
                                        double_array_3d& output_cond,
                                        bool logcomputation = false,
                                        int index = I_DEFAULT) const;

   /** Compute part of the loglikelihood related initial states */
   double initial_state_likelihood_computation(const MarkovOutTreeData& trees) const;

   /** Compute part of the loglikelihood related to output processes */
   double output_likelihood_computation(const MarkovOutTreeData& trees) const;

   /** Compute part of the loglikelihood related to ordered children */
   double ordered_model_likelihood_computation(const MarkovOutTreeData& trees) const;

   /** Compute part of loglikelihood related to unordered children
       without accounting for possible permutations of subtrees */
   double unordered_model_partial_likelihood_computation(const MarkovOutTreeData& trees) const;

   /** Compute part of loglikelihood related to unordered children
       accounting for possible permutations of subtrees */
   double unordered_model_likelihood_computation(const MarkovOutTreeData& trees) const;

   // methods related to parameter estimation

   /** Estimate variable order Markov models */
   void vomc_estimation(StatError &error, const MarkovOutTreeData& trees);

   /** Estimate generation distributions */
   void generation_process_estimation(StatError &error,
                                      const MarkovOutTreeData& trees,
                                      const std::vector<int>& generation_types,
                                      unsigned int generation_min_inf_bound = 0,
                                      bool generation_min_inf_bound_flag = false);

   /** Estimate output processes */
   void output_process_estimation(StatError &error, const MarkovOutTreeData& trees,
                                  const std::vector<bool> &skipped_variables,
                                  bool common_dispersion = false);

   /** Estimate initial distribution and update vomcs if necessary
       (should be called after vomc_estimation) */
   void initial_distribution_estimation(const MarkovOutTreeData& trees);


public :

   /** Create a MarkovOutTree with empty attributes */
   MarkovOutTree();

   /** Create a MarkovOutTree from the number of states,
       the number and values of observation processes, the number
       of ordered children, the number and values of Markov models
       for ordered children, and the branching processes */
   MarkovOutTree(int inb_state, int inb_output_process, CategoricalProcess * const * npprocess = NULL,
                 DiscreteParametricProcess * const * dpprocess = NULL, ContinuousParametricProcess * const * cprocess = NULL,
                 unsigned int iunordered_rank = 0, unsigned int inb_vomc = 0,
                 sequence_analysis::VariableOrderMarkov * const * pvomc_models = NULL,
                 const GenerationProcess *pgeneration = NULL);

   /** Create a  MarkovOutTree from Markov chain,
       the number and values of observation processes, the number
       of ordered children, the number and values of Markov models
       for ordered children, and the branching processes */

   MarkovOutTree(const Chain& markov, int inb_output_process, CategoricalProcess * const * npprocess = NULL,
                 DiscreteParametricProcess * const * dpprocess = NULL, ContinuousParametricProcess *const * cprocess = NULL,
                 unsigned int iunordered_rank = 0, unsigned int inb_vomc = 0,
                 sequence_analysis::VariableOrderMarkov * const *  pvomc_models = NULL,
                 const GenerationProcess *pgeneration = NULL);

   /** Create a MarkovOutTree the number and values of observation processes,
       copying the remaining fields from a given MarkovOutTree */
   MarkovOutTree(const MarkovOutTree * markov, int inb_output_process,
                 CategoricalProcess * const * nonparametric_observation,
                 DiscreteParametricProcess * const * discrete_parametric_observation,
                 ContinuousParametricProcess * const * continuous_parametric_observation,
                 int depth = I_DEFAULT_TREE_SIZE);

   /** Copy constructor of MarkovOutTree class */
   MarkovOutTree(const MarkovOutTree& markov, bool data_flag = true);

   /** Destructor for MarkovOutTree class */
   virtual ~MarkovOutTree();

   /** Return the data part of a MarkovOutTree,
       keeping a reference on \e self */
   MarkovOutTreeData* extract_data(StatError& error) const;
   /** Return the number of states */
   unsigned int get_nb_state() const;
   /** Return the number of output processes */
   unsigned int get_nb_output_process() const;
   /** Return the number of nonparametric output processes */
   unsigned int get_nb_nonparametric_output_process() const;
   /** Return the number of discrete parametric output processes */
   unsigned int get_nb_discrete_parametric_output_process() const;
   /** Return the number of continuous (parametric) output processes */
   unsigned int get_nb_continuous_output_process() const;
   /** Print MarkovOutTree on a single line */
   ostream& line_write(ostream& os) const;
   /** Print model with full specification */
   // method inherited from StatInterface, pure virtual if not redefined
   ostream& ascii_write(ostream& os, bool exhaustive = false) const;
   /** Print distribution into a file */
   bool ascii_write(StatError& error, const char * path,
                    bool exhaustive = false) const;

   /** Spreadsheet-like output of the MarkovOutTree model
       into a file */
   bool spreadsheet_write(StatError& error, const char * path) const;

   /** Return the number of views (i.e. the size) in Matplotlib output
       due to generation process */
   unsigned int nb_generation_plot_set_computation() const;
   /** Gnuplot output of the MarkovOutTree model*/
   bool plot_write(StatError& error, const char * prefix,
                   const char * title = NULL) const;
   /** Matplotlib output of the MarkovOutTree model */
   MultiPlotSet* get_plotable() const;

   /** Compute characteristic distributions of \e self */
   void characteristic_computation(int depth, bool counting_flag,
                                   int variable = I_DEFAULT);

   /** Compute characteristic distributions of \e self
       (using a dataset) */
   void characteristic_computation(const MarkovOutTreeData& tree,
                                   bool counting_flag,
                                   int variable = I_DEFAULT,
                                   bool depth_flag = true);

   /** Get the number of ordered children involved in the branching processes
       as factors (due to dependencies regarding these states,
       excluding parent state) */
   unsigned int get_nb_children_branching() const;

   /** Get the number of other factors involved in the branching process */
   unsigned int get_nb_factors() const;

   /** Get the number of possible values for each factor */
   GenerationProcess::index_array get_factor_values() const;

   /* Get the number of generation processes
      (combining all possible values of every factor) */
   unsigned int get_nb_generation_values() const;

   /** Get the number of factors for generation processes
       (including parent state) */
   unsigned int get_nb_generation() const;

   /** Get the number of variable order Markov chains */
   unsigned int get_nb_vomc() const;

   /** Get the number of ordered children */
   unsigned int get_nb_ordered_children() const;

   /** Get the min offsets for every generation process */
   std::vector<unsigned int> get_generation_min_offsets() const;

   /** Return the marginal distributions associated with a generation process
       (a new instance is allocated) */
   std::vector<Distribution*> get_distribution(StatError &error,
                                               const GenerationProcess::index_array &fact) const;

   /** Return generation process associated with a given factor
        (return a pointer; object should not be deallocated) */
   const DiscreteMultivariateDistribution* get_generation_ptr(StatError &error,
                                                              const GenerationProcess::index_array &fact) const;

   /** Permutation of the states of \e self*/
   void state_permutation(StatError& error, int* perm) const;

   /** Log likelihood computation for a single tree */
   double likelihood_computation(const Trees& trees, unsigned int state_variable,
                                 const std::vector<unsigned int>* factor_indices=NULL,
                                 int index = I_DEFAULT) const;

   /** Log likelihood computation for a whole forest */
   double likelihood_computation(const MarkovOutTreeData& trees,
                                 const std::vector<unsigned int>* factor_indices=NULL,
                                 int index=I_DEFAULT) const;

   /** Log likelihood computation for a single tree,
       accounting or not for possible permutations of subtrees */
   double likelihood_computation(const MarkovOutTreeData& trees, bool partial,
                                 const std::vector<unsigned int>* factor_indices,
                                 int index) const;

   /** Simulation of trees using specified tree depths */
   MarkovOutTreeData* simulation(StatError& error,
                                 const std::vector<unsigned int>& depths,
                                 const std::vector<unsigned int>& factors,
                                 bool counting_flag = true,
                                 bool divergence_flag = false) const;

   /** Simulation of trees using an empirical distribution for
       tree depths */
   MarkovOutTreeData* simulation(StatError& error, const FrequencyDistribution& idepth,
                                 bool depth_flag = true,
                                 bool counting_flag = true,
                                 bool divergence_flag = false) const;

   /** Simulation of trees using maximal values for tree depth or size */
   MarkovOutTreeData* simulation(StatError& error, int inb_trees,
                                 int imax_depth, int imax_size,
                                 bool depth_flag = true,
                                 bool counting_flag = true) const;

   /** Compute an adaptive penalty for \e self used for model selection */
   double penalty_computation(double min_probability = 0.) const;

   // methods below should be used for testing purposes only...
   void get_state_marginal_distribution(const MarkovOutTreeData& trees,
                                        double_array_3d& state_marginal) const;
   void get_output_conditional_distribution(const MarkovOutTreeData& trees,
                                            double_array_3d& output_cond) const;

};

/** Parse a MarkovOutTree from an input stream */
MarkovOutTree* markov_out_tree_parsing(StatError& error,
                                       std::ifstream &in_file,
                                       int &line);

/** Parse a MarkovOutTree from a file */
MarkovOutTree* markov_out_tree_ascii_read(StatError& error,
                                          const char * path,
                                          int size = I_DEFAULT_TREE_SIZE,
                                          double cumul_threshold = sequence_analysis::OCCUPANCY_THRESHOLD);
/**
   \class MarkovOutTreeData
   \brief Set of tree-structured data potentially associated with
          a MarkovOutTree statistical model (Hidden or not)
*/

class MarkovOutTreeData : public Trees
{
   // friend classes
   friend class MarkovOutTree;
   template <typename DummyType> friend class MarkovOutTreeReestimation;

   friend std::ostream& operator<<(std::ostream &os, const MarkovOutTreeData& trees);

public :

   typedef Trees::tree_type tree_type;
   typedef tree_type::value value;

   typedef Typed_edge_one_int_tree::tree_type state_tree_type;
   typedef state_tree_type::value state_value;

   typedef tree_type::key key;
   typedef tree_type::vertices_size_type vertices_size_type;
   typedef tree_type::children_iterator children_iterator;
   typedef tree_type::vertex_iterator vertex_iterator;

   typedef FrequencyDistribution*** ptFrequencyDistribution_array_2d;
   typedef TreeCharacteristics::ptFrequencyDistribution_array FrequencyDistribution_array;
   typedef Histogram*** ptHistogram_array_2d;
   typedef TreeCharacteristics::ptHistogram_array ptHistogram_array;
   typedef Typed_edge_one_int_tree** ptOne_int_tree_set;
   typedef std::vector<HiddenMarkovTreeData*> pt_hmtd_vector;
   typedef double** double_array_2d;
   typedef std::map<int, bool> virtual_vdic;
   typedef GenerationProcess::index_array index_array;
   typedef Statiskit::Marginal::Multivariate::CountTable mvHistogram;
   typedef Statiskit::Marginal::Univariate::CountTable uHistogram;
   typedef Statiskit::Marginal::Univariate::CountTable::ptr_type uHistogram_ptr_type;

protected :

   /// the Markov out tree model
   MarkovOutTree* markov;
   /// initial state and transitions
   ChainData* chain_data;
   /// data part for sequence_analysis::VariableOrderMarkov model
   sequence_analysis::VariableOrderMarkovData** vomc_data;
   // multivariate histogram for joint state empirical distribution
   // DiscreteMultivariateReestimation<int> **empirical_observation;
   /// FrequencyDistribution of tree depths
   FrequencyDistribution *hdepth;
   /// likelihood for the given trees
   double likelihood;
   /// completed likelihood for the given trees (hidden Markov case)
   double hidden_likelihood;
   /// number of states
   int _nb_states;

   /// state trees
   ptOne_int_tree_set state_trees;

   /// dictionaries of virtual leaf vertices
   /// i.e. with censored number of children or missing vertices
   virtual_vdic** virt_vertices;

   /// FrequencyDistribution corresponding to the conditional observation distribution,
   /// depending on the considered observed (integral) variable
   /// and on the value of the state variable
   ptFrequencyDistribution_array_2d observation_distribution;
   /// Histogram corresponding to the conditional observation distribution,
   /// depending on the considered observed (continuous) variable
   /// and on the value of the state variable
   ptHistogram_array_2d observation_histogram;

    /// frequency distribution of the characteristic quantities
    /// for the state variable
   TreeCharacteristics *state_characteristics;

   /// frequency distributions for joint distribution of children states
    MarkovOutTreeReestimation<int> *mot_reestimation;
    // mvHistogram *mot_reestimation;

   /** Copy operator */
   void copy(const MarkovOutTreeData& otrees, bool model_flag = true,
             bool characteristic_flag = true);

   /** Destructor. The markov part is deleted if any. */
   void remove();

   /** compute the number of states if unknown */
   void nb_state_computation();

   /** Compute min and max values of variables */
   virtual void min_max_value_computation();
   /** Compute maximal number of vertices */
   virtual void max_size_computation();
   /** Compute maximal depth */
   virtual void max_depth_computation();
   /** Compute maximal number of children */
   virtual int max_nb_children_computation();
   /** Compute total number of vertices */
   virtual int cumul_size_computation() const;
   /** Compute histogram of sizes */
   virtual void build_size_frequency_distribution();
   /** Compute total number of children */
   virtual int cumul_nb_children_computation() const;
   /** Compute histogram for the number of children,
       taking into account virtual vertices */
   virtual void build_nb_children_frequency_distribution();

   /** Print a summary of the data using specified format */
   ostream& ascii_write(ostream& os, bool exhaustive,
                        bool comment_flag) const;

   /** Return the number of views (i.e. the size) in Matplotlib output
       due to generation processes */
   unsigned int generation_nb_plot_set_computation() const;

   /** Write views due to generation processes into a MultiPlotSet object */
   void generation_plotable_write(MultiPlotSet &plot_set, int &index) const;

   /** Write views due to discrete output processes
       into a MultiPlotSet object, for a given variable */
   void discrete_output_plotable_write(MultiPlotSet &plot_set,
                                       int &index, int var) const;

   /** Compute the frequency distributions for the conditional distribution
       for one observed variable and for each state */
   void observation_frequency_distribution_computation(int ivariable);

   /** Compute the frequency distributions of the characteristic quantities
       for the state variable */
   void build_state_characteristics();

   /** Compute the frequency distributions of intial states */
   void build_chain_reestimation(char itype, int inb_state);

   /** Compute the frequency distributions of the joint distribution
      of children states */
   void build_markov_reestimation(unsigned int inb_children_branching = 0,
                                  unsigned int inb_factors = 0,
                                  const GenerationProcess::index_array *ifactor_values = NULL,
                                  unsigned int inb_generation = 0,
                                  unsigned int iunordered_rank = 0,
                                  unsigned int inb_vomc = 0,
                                  const std::vector<unsigned int>* factor_indices = NULL);

   /** Permutation of the states of \e self */
   void state_permutation(int* perm);

   // operations on state trees
   /** Select subset of state trees from their indices */
   ptOne_int_tree_set state_trees_select_individual(StatError& error, int inb_trees,
                                                   int_array iidentifier,
                                                   int_array& index,
                                                   bool keep = true) const;

public :

   /** Constructor from the number of integral and float variables,
       and the number of trees (composed of a single vertex) */
   MarkovOutTreeData(unsigned int inb_integral, unsigned int inb_float,
                     unsigned int inb_trees);

   /** Copy constructor */
   MarkovOutTreeData(const MarkovOutTreeData& trees,
                     bool model_flag,
                     bool characteristic_flag);

   /** Constructor by conversion from a set of trees */
   MarkovOutTreeData(int inb_trees,
                     const int * itype,
                     Default_tree **otrees);

   /** Constructor by conversion from a Trees instance */
   MarkovOutTreeData(const Trees& otrees, int istate_variable = I_DEFAULT);

   /** Destructor */
   virtual ~MarkovOutTreeData();

   /** Return the size of a given tree of \e self */
   virtual unsigned int get_size(int index) const;

   // Return total number of vertices
   // virtual unsigned int get_total_size() const;

   /** Allocate the FrequencyDistributions associated with the output distributions */
   void create_observation_frequency_distribution(int nb_state);
   /** Compute the FrequencyDistributions associated with the output distributions */
   void observation_frequency_distribution_computation();
   /** Allocate and compute the FrequencyDistributions associated with the output distributions */
   void build_observation_frequency_distribution();
   /** Allocate and compute the FrequencyDistributions associated with state process */
   void build_state_frequency_distribution(StatError &error,
                                           unsigned int inb_children_branching = 0,
                                           unsigned int inb_factors = 0,
                                           const GenerationProcess::index_array *ifactor_values = NULL,
                                           unsigned int inb_generation = 0,
                                           unsigned int iunordered_rank = 0,
                                           unsigned int inb_vomc = 0,
                                           const std::vector<unsigned int>* factor_indices = NULL);

   // access to class members
   /** Return number of states */
   int get_nb_states() const;
   /** Return model associated with the data (a new instance is allocated) */
   MarkovOutTree* get_markov() const;
   /** Return model associated with the data (return a pointer; object should not be deallocated) */
   MarkovOutTree* get_markov_ptr() const;

   /** Return the set of state trees (a new instance is allocated) */
   ptOne_int_tree_set get_state_trees() const;
   /** Return the set of state trees (return a pointer; object should not be deallocated) */
   ptOne_int_tree_set get_state_trees_ptr() const;

   /** Return a given state tree (a new instance is allocated) */
   Typed_edge_one_int_tree* get_state_tree(int itree) const;
   /** Return a given state tree (return a pointer; object should not be deallocated) */
   Typed_edge_one_int_tree* get_state_tree_ptr(int itree) const;

   /** Return a MarkovOutTreeData containing the states
       as a variable */
   MarkovOutTreeData* get_state_markov_out_tree_data() const;

   /** Return true iif a vertex is virtual */
   bool is_virtual(int itree, int ivertex) const;

   /** Return dictionary of virtual vertices for a given tree
       (a new instance is allocated) */
   virtual_vdic* get_virtual_vertices(int itree) const;

   /** Return dictionary of virtual vertices for a given tree
       (return a pointer; object should not be deallocated) */
   const virtual_vdic * get_virtual_vertices_ptr(int itree) const;

   /** Set virtual status of a vertex */
   void set_virtual_vertex(int itree, int ivertex, bool bvirtual = true) const;

   /** Update histogram for the number of children
       and minimal / maximal values of variables,
       taking into account virtual vertices */
   void update_frequency_distributions();

   /** Get every possible combination of factors
       for generation processes (including parent state).
       Parent state comes first (?) if included */
   std::vector<GenerationProcess::index_array> get_factor_combinations() const;

   /** Return the marginal distributions associated with a generation process
       (a new instance is allocated) */
   std::vector<DiscreteDistributionData*> get_distribution_data(StatError &error,
                                                                const index_array &fact) const;
   /** Return observed vectors associated with a generation process
       together with their weights */
   std::vector<std::pair<std::vector<unsigned int>, unsigned int> >
      get_joint_distribution(StatError &error,
                             const index_array &fact) const;

   /** Print summary of the data */
   // method inherited from Trees, itself inherited from StatInterface
   ostream& ascii_write(ostream& os, bool exhaustive = false) const;

   // Print data into a file
   // Not necessary: already in Trees
   /* bool ascii_write(StatError& error, const char * path,
                    bool exhaustive = false) const;*/

   /** Matplotlib output of MarkovOutTreeData */
   MultiPlotSet* get_plotable() const;

   /** Update the frequency distributions of the joint distribution
      of children states */
   void update_markov_reestimation(StatError& error,
                                   unsigned int inb_vomc = 0,
                                   unsigned int inb_ordered_children = 0,
                                   unsigned int inb_children_branching = 0,
                                   const std::vector<unsigned int>& factor_variables = std::vector<unsigned int>(),
                                   bool parent_dependent = true);

   // operations on trees
   /** Select subset of trees from their indices */
   MarkovOutTreeData* select_individual(StatError& error, int inb_trees,
                                        int_array iidentifier, bool keep = true) const;


   /** Check whether a given family of MarkovOutTrees is valid
       on current data set */
   bool markov_out_tree_check_estimation(StatError& error,
                                         unsigned int inb_vomc,
                                         unsigned int inb_ordered_children,
                                         unsigned int inb_children_branching,
                                         const std::vector<unsigned int>& factor_variables,
                                         bool parent_dependent,
                                         const std::vector<int>& generation_types,
                                         unsigned int& inb_generation,
                                         unsigned int& inb_generation_values,
                                         std::vector<int>& generation_types_estimation,
                                         const GenerationProcess::index_array* ifactor_values,
                                         std::vector<bool>& skipped_variables) const;
   // model estimation
   /** Estimate a MarkovOutTree */
   MarkovOutTree* markov_out_tree_estimation(StatError& error, unsigned int inb_vomc,
                                             unsigned int inb_ordered_children,
                                             unsigned int inb_children_branching,
                                             const std::vector<unsigned int>& factor_variables = std::vector<unsigned int>(),
                                             bool parent_dependent = true,
                                             const std::vector<int>& generation_types = std::vector<int>(),
                                             unsigned int generation_min_inf_bound = 0,
                                             bool generation_min_inf_bound_flag = true) const;

};

/**
   \class MarkovOutTreeReestimation
   \brief Data structure associated with MarkovOutTree models
   for exploratory analysis and estimation
*/

template <typename Type>
class MarkovOutTreeReestimation
{  /// Data structure associated with MarkovOutTree models
   /// for exploratory analysis and estimation

public :

   // public types

   /// configuration of unordered children states
   /// or factors
   typedef GenerationProcess::index_array index_array;

   /// Configurations of factors are in the following order:
   /// Parent state, ordered states, other factors

   typedef Statiskit::Marginal::Multivariate::CountTable mvHistogram;
   typedef Statiskit::Marginal::Univariate::CountTable uHistogram;
   // public types of DiscreteMultivariateReestimation<Type> may be required
   /* // public types
   /// number of observation vectors for each value
   /// key is the value of the vector coded in base max_nb_value
   typedef std::map<unsigned int, Type > values_dict;
   /// mapping between values and code
   typedef std::map<value, unsigned int > values_code;
   /// mapping between code and values
   typedef std::map<unsigned int, value > values_decode;
   /// iterator on values
   typedef typename values_dict::iterator iterator; */


protected :

   /// number of states
   int nb_state;
   /// number of ordered children involved in the branching processes
   /// as factors (due to dependencies regarding these states, excluding parent state)
   unsigned int nb_children_branching;
   /// number of other factors involved in the branching process
   unsigned int nb_factors;
   /// number of possible values for each factor
   index_array *factor_values;
   /// number of generation processes
   /// (combining all possible values of every factor)
   unsigned int nb_generation_values;
   /// number of factors for generation processes
   /// (including parent state)
   unsigned int nb_generation;
   /// base to code and decode configuration of factors
   index_array factor_base;
   /// rank (respective with partial order of children) beyond which
   /// children are unordered
   unsigned int unordered_rank;
   /// number of variable order Markov chain models on ordered children
   unsigned int nb_vomc;
   /// sequences of orderered children (one for each parent state)
   sequence_analysis::MarkovianSequences **seq;
   /// empirical generation distributions: generation_reestim[code]
   /// where code denotes an ordered set of factors
   // DiscreteMultivariateReestimation<Type> **generation_reestim;
   mvHistogram **generation_reestim;

   /** Copy operator */
   void copy(const MarkovOutTreeReestimation<Type>& process);

   /** Destructor */
   void remove();

   /** Compute FrequencyDistribution nb_element, mean, etc.
       in empirical generation distributions */
   // void frequency_distribution_computation() const;

public:

   /** Create a MarkovOutTreeReestimation using the number of states,
       of factors and parent states determining the number inb_generation
       of branching processes and the empirical generation distributions */
   MarkovOutTreeReestimation(int inb_state = 0,
                             unsigned int inb_children_branching = 0,
                             unsigned int inb_factors = 0, const index_array *ifactor_values = NULL,
                             unsigned int inb_generation = 0,
                             unsigned int unordered_rank = 0,
                             unsigned int nb_vomc = 0,
                             mvHistogram * const * pgeneration = NULL,
                             const sequence_analysis::MarkovianSequences * const * pseq = NULL);
   /** Copy constructor of MarkovOutTreeReestimation class */
   MarkovOutTreeReestimation(const MarkovOutTreeReestimation& process);
   /** Destructor for MarkovOutTreeReestimation class */
   ~MarkovOutTreeReestimation();

   /** Code a configuration of factors with an unsigned int */
   unsigned int code(const index_array& fact) const;
   /** Decode a configuration of factors coded with an unsigned int */
   index_array decode(unsigned int code) const;

   /** Code a configuration of factors with an unsigned int
       with control of errors */
   unsigned int code_error(StatError& error, const index_array& fact) const;

   /** Get number of states */
   int get_nb_state() const;

   /** Get rank (strictly) beyond which children are unordered */
   int get_unordered_rank() const;

   /** Get number of variable order Markov chains */
   int get_nb_vomc() const;

   /** Get number of ordered children involved in the branching processes
       as factors (due to dependencies regarding these states,
       excluding parent state) */
   unsigned int get_nb_children_branching() const;

   /** Get number of other factors involved in the branching process */
   unsigned int get_nb_factors() const;

   /** Get number of possible values for each factor */
   index_array get_factor_values() const;

   /** Get number of generation processes
      (combining all possible values of every factor) */
   unsigned int get_nb_generation_values() const;

   /** Get number of factors for generation processes
       (including parent state) */
   unsigned int get_nb_generation() const;

   /** Get every possible combination of factors
       for generation processes (including parent state).
       Parent state comes first (?) if included */
   std::vector<GenerationProcess::index_array> get_factor_combinations() const;

   /** Permutation of the states of \e self*/
   void state_permutation(StatError& error, int* perm) const;

   /** Get marginal (return a pointer; object should not be deallocated) */
   // const Reestimation< Type >* get_marginal_ptr(unsigned int code, int istate) const;

   /** Get marginal (a new instance is allocated) */
   uHistogram* get_marginal(unsigned int code, int istate) const;

   /** Get marginal FrequencyDistribution (a new instance is allocated) */
   // FrequencyDistribution* get_marginal_frequency(unsigned int code, int istate) const;

   /** Update frequency for a given configuration of factors and states */
   void update(const index_array& fact, const index_array& states, Type freq);

   /** Get frequency for a given configuration of factors and states */
   Type get_frequency(const index_array& fact, const index_array& states) const;

   /** Get reestimation for generation process
       (return a pointer; object should not be deallocated) */
   mvHistogram** get_generation_reestim_ptr() const;

   /** Get reestimation for generation process
       (return a copy; object should be deallocated) */
   mvHistogram** get_generation_reestim() const;

   /** Get reestimation for generation process
       (return a pointer; object should not be deallocated) */
    mvHistogram* get_generation_reestim_ptr(StatError& error,
                                            const index_array& fact) const;

   /** Get reestimation for generation process
       (return a copy; object should be deallocated) */
   mvHistogram* get_generation_reestim(StatError& error,
                                       const index_array& fact) const;

   /** Initialize empirical generation distributions */
   void init(const std::vector<unsigned int>& imin_offsets,
             const std::vector<unsigned int>& imax_nb_values);

   /** Compute frequencies for every configurations of states */
   void build_characteristics(const MarkovOutTreeData& data,
                              int index = I_DEFAULT,
                              const std::vector<unsigned int>* factor_indices = NULL);


};

#include "markov_out_tree.hpp"


}; // end namespace

#endif
