/*------------------------------------------------------------------------------
 *
 *        VPlants.Tree-Statistic : VPlants Tree-Statistic module
 *        HiddenMarkovTrees
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_hmt.cpp 9099 2010-06-08 09:03:00Z pradal $
 *
 *-----------------------------------------------------------------------------*/
// Includes ====================================================================
#include "tree/basic_visitors.h"
#include "tree/tree_traits.h"
#include "tree/tree_simple.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"


// required by "tree_statistic/multivariate_distribution.h"
#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_label.h"
#include "tree_statistic/tree_labels.h"
#include "statiskit/core/data/marginal/multivariate.h"
#include "tree_statistic/mv_histogram_tools.h"

#include "tree_statistic/multivariate_distribution.h"
#include "tree_statistic/tree_labels.h"
#include "tree_statistic/int_fl_containers.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include "tree_statistic/typed_edge_trees.h"
#include "tree_statistic/hidden_markov_tree.h"
#include "tree_statistic/markov_out_tree.h"

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::len
#include <boost/python/make_constructor.hpp>
// definition of boost::python::make_constructor

#include "stat_tool/wrapper_util.h"
#include "stat_tool/boost_python_aliases.h"

#include "../errors.h"

// Using =======================================================================
using namespace boost::python;
using namespace sequence_analysis;
using namespace Stat_trees;
using namespace tree_statistic;

// Declarations ================================================================
template<int num> struct UniqueInt { int v; enum { value=num };
UniqueInt(int _v) : v(_v) { } operator int() const { return v; } };

template <class T>
struct make_object {
    typedef T obj_type;

    const obj_type& __c_obj;

    make_object(const T& c_obj) : __c_obj(c_obj) {}
    boost::python::object operator()() const { return boost::python::object(__c_obj); }
};

struct dict_converter {
    typedef std::map<int, bool> dict_type;
    typedef std::map<int, bool>::const_iterator dict_const_iterator;
    typedef std::map<int, bool>::key_type dict_key_type;
#if (defined(USING_UNORDERED_MAP)) || defined(WIN32_STL_EXTENSION) || defined(__GNUC__)
    typedef make_object<std::map<int, bool>::mapped_type> ValueTranslator;
    typedef make_object<std::map<int, bool>::key_type> KeyTranslator;
    typedef std::map<int, bool>::mapped_type ValueType;
#else
    typedef make_object<std::map<int, bool>::data_type> ValueTranslator;
    typedef std::map<int, bool>::data_type ValueType;
#endif


    dict_const_iterator __c_dict_begin;
    dict_const_iterator __c_dict_end;

    dict_converter(const std::map<int, bool>& c_dict): __c_dict_begin(c_dict.begin()),__c_dict_end(c_dict.end()){}
    dict_converter(const dict_const_iterator& c_dict_begin,const dict_const_iterator& c_dict_end):
        __c_dict_begin(c_dict_begin),__c_dict_end(c_dict_end){}

    boost::python::dict convert() const {
        boost::python::dict d;
        for (dict_const_iterator it = __c_dict_begin; it != __c_dict_end; ++it)
            d[KeyTranslator(it->first)()] = ValueTranslator(it->second)();
        return d;
    }
    inline boost::python::object operator()() const { return convert(); }
    inline operator boost::python::object () const { return convert(); }
};

#define WRAP MarkovOutTreeWrap
class WRAP
{

public :

   static MarkovOutTreeData*
   Mot_wrapper_extract_data(const MarkovOutTree& markov)
   {
      MarkovOutTreeData *res = NULL;
      StatError error;
      ostringstream error_message;

      res = markov.extract_data(error);
      if (res == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static str Mot_wrapper_display(const MarkovOutTree& markov,
                                  bool exhaustive)
   {
      std::stringstream s;
      str res;

      markov.ascii_write(s, exhaustive);
      res = str(s.str());
      return res;
   }

   static bool Mot_wrapper_file_ascii_write2(const MarkovOutTree &markov,
                                             const char * path, bool exhaustive)
   {
      StatError error;
      bool res;
      ostringstream error_message;

      res = markov.ascii_write(error, path, exhaustive);
      if (!res)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static bool Mot_wrapper_file_ascii_write1(const MarkovOutTree &markov,
                                             const char * path)
   { return Mot_wrapper_file_ascii_write2(markov, path, false); }


   /*
   static int Mot_wrapper_nb_values(const MarkovOutTree &markov, int variable)
   {
      int res = 0;
      ostringstream error_message;

      if ((variable >= 0) && (variable < markov.get_nb_output_process()))
         res = markov.get_nb_values(variable);
      else
      {
         error_message << "Bad variable: " << variable << endl;
         throw_stat_tree_error(error_message);
      }
      return res;
   } */

   static MarkovOutTreeData*
   Mot_wrapper_simulate_depths(const MarkovOutTree &markov,
                               boost::python::list &factor_indices,
                               boost::python::list &depths)
   {
      StatError error;
      MarkovOutTreeData* markov_data= NULL;
      ostringstream error_message;
      const unsigned int nb_trees = boost::python::len(depths);
      const unsigned int nb_factors = boost::python::len(factor_indices);
      int i = 0;
      std::vector<unsigned int> dlist(nb_trees,0);
      std::vector<unsigned int> flist(nb_factors);

      if (nb_trees == 0)
         throw_stat_tree_error("Input list of depths cannot be empty");

      for(i = 0; i < nb_trees; i++)
         dlist[i] = boost::python::extract<unsigned int>(depths[i]);

      for(i = 0; i < nb_factors; i++)
         flist[i] = boost::python::extract<unsigned int>(factor_indices[i]);

      markov_data = markov.simulation(error, dlist, flist, false, false);


      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return markov_data;
   }

   // static MarkovOutTree*
   static boost::shared_ptr<MarkovOutTree>
   Mot_wrapper_ascii_read(const char * path, int depth = I_DEFAULT_TREE_SIZE,
                          double cumul_threshold = OCCUPANCY_THRESHOLD)
   {
      StatError error;
      MarkovOutTree *mot = NULL;
      ostringstream message;

      mot = Stat_trees::markov_out_tree_ascii_read(error, path, depth,
                                                   cumul_threshold);
      if (mot == NULL)
      {
         message << error;
         throw_stat_tree_error(message);
      }

      return boost::shared_ptr<MarkovOutTree>(mot);
   }

   // to be used with return_value_policy< manage_new_object >()
   // use DEF_RETURN_VALUE otherwise ?
   static MarkovOutTree*
   Mot_wrapper_ascii_read3(char *filename, int depth,
                           double cumul_threshold)
   {
      StatError error;
      ostringstream message;

      MarkovOutTree *mot = NULL;

      mot = markov_out_tree_ascii_read(error, filename, depth,
                                       cumul_threshold);

      if (mot == NULL)
      {
         message << error;
         throw_stat_tree_error(message);
      }

      return mot;
   }

   // to be used with return_value_policy< manage_new_object >()
   // use DEF_RETURN_VALUE otherwise ?
   static MarkovOutTree*
   Mot_wrapper_ascii_read2(char *filename, int depth)
   { return Mot_wrapper_ascii_read3(filename, depth, OCCUPANCY_THRESHOLD); }

   // to be used with return_value_policy< manage_new_object >()
   // use DEF_RETURN_VALUE otherwise ?
   static MarkovOutTree*
   Mot_wrapper_ascii_read1(char *filename)
   { return Mot_wrapper_ascii_read2(filename, I_DEFAULT_TREE_SIZE); }

   static unsigned int Mot_wrapper_nb_state(const MarkovOutTree &markov, int dummy_variable = 0)
   { return markov.get_nb_state(); }

   static str Mot_wrapper_ascii_write0(const MarkovOutTree& mot)
   {
      std::stringstream s;
      str res;

      mot.line_write(s);
      res = str(s.str());
      return res;
   }

   static MultiPlotSet*
   Mot_wrapper_get_plotable(const MarkovOutTree& markov)
   {
     StatError error;
     MultiPlotSet *ret = markov.get_plotable();
     if (ret == NULL)
       stat_tool::wrap_util::throw_error(error);
     return ret;
   }

}; // WRAP


// Module ======================================================================
void class_mot()
{
  class_< MarkovOutTree >
    ("_MarkovOutTree", "MarkovOutTree", init< const MarkovOutTree&, optional< bool> >())
        .def("__init__", make_constructor(WRAP::Mot_wrapper_ascii_read))
        .def("__init__", make_constructor(WRAP::Mot_wrapper_ascii_read1))

        /* .def("ExtractData", WRAP::Mot_wrapper_extract_data,
                            return_value_policy< manage_new_object >()) */
        DEF_RETURN_VALUE_NO_ARGS("ExtractData", WRAP::Mot_wrapper_extract_data,
                                 "ExtractData(self) -> _MarkovOutTreeData  \n\n")
        .def("Display", WRAP::Mot_wrapper_display,
                        "Display(self, bool) -> str \n\n"
                        "Print MarkovOutTree definition")
        .def("FileAsciiWrite", WRAP::Mot_wrapper_file_ascii_write1)
        .def("FileAsciiWrite", WRAP::Mot_wrapper_file_ascii_write2)
        /* .def("Simulate", WRAP::Mot_wrapper_simulate_trees,
                         return_value_policy< manage_new_object >()) */

        /* .def("_LineWrite", WRAP::Mot_wrapper_ascii_write0,
                           "_LineWrite(self) -> str \n\n"
                           "Print Markov out-tree on a single line\n") */

        DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::Mot_wrapper_get_plotable,
                                 "return plotable")

        .def("_NbDiscreteParametric", &MarkovOutTree::get_nb_discrete_parametric_output_process,
                        "_NbDiscreteParametric(self) -> int \n\n"
                        "return the number of discrete parametric output processes "
                        "of the Markov out-tree\n")

        .def("NbFloat", &MarkovOutTree::get_nb_continuous_output_process,
                        "NbFloat(self) -> int \n\n"
                        "return the number of continuous output processes "
                        "of the Markov out-tree\n")

        .def("_NbGenerationPlotSet", &MarkovOutTree::nb_generation_plot_set_computation,
                                     "_NbGenerationPlotSet(self) -> int \n\n"
                                     "return the number of plot views related to "
                                     "generation processes of the Markov out-tree\n")

        .def("_NbNonParametric", &MarkovOutTree::get_nb_nonparametric_output_process,
                        "_NbNonParametric(self) -> int \n\n"
                        "return the number of discrete nonparametric output processes "
                        "of the Markov out-tree\n")

        .def("NbStates", &MarkovOutTree::get_nb_state,
                         "NbStates(self) -> int \n\n"
                         "return the number of states "
                         "of the Markov out-tree\n")

        /*
        .def("NbValues", WRAP::Mot_wrapper_nb_values,
                         "NbValues(self) -> int \n\n"
                         "return the number of values for a given process "
                         "of the Markov out-tree\n")*/

        .def("NbVariable", &MarkovOutTree::get_nb_output_process,
                          "NbVariable(self) -> int \n\n"
                          "return the number of output processes "
                          "of the Markov out-tree\n")

        DEF_RETURN_VALUE("Simulate", WRAP::Mot_wrapper_simulate_depths,
                         args("NbTrees", "Factors", "Depths"),
                         "Simulate(self) -> _MarkovOutTreeData  \n\n")

        .def(self_ns::str(self)) //__str__
        .def("__str__", WRAP::Mot_wrapper_ascii_write0)
    ;

    def("MotAsciiRead", WRAP::Mot_wrapper_ascii_read1,
                        return_value_policy< manage_new_object >());

};

#undef WRAP

#define WRAP MarkovOutTreeDataWrap
class WRAP
{

public :

   static str Mot_data_wrapper_display(const MarkovOutTreeData& data,
                                       bool exhaustive)
   {
      std::stringstream s;
      str res;

      data.ascii_write(s, exhaustive);
      res = str(s.str());
      return res;
   }

   static MarkovOutTree* Mot_data_wrapper_markov_out_tree_estimation1(const MarkovOutTreeData& data,
                                                                      unsigned int inb_vomc,
                                                                      unsigned int inb_ordered_children,
                                                                      unsigned int inb_children_branching,
                                                                      boost::python::list &factor_variables,
                                                                      bool parent_dependent,
                                                                      boost::python::list &generation_types,
                                                                      unsigned int generation_min_inf_bound,
                                                                      bool generation_min_inf_bound_flag)
   {
      bool status = true;
      unsigned int i;
      StatError error;
      std::vector<unsigned int> vfactor_variables, error_indices;
      std::vector<int> vgeneration_types;
      ostringstream error_message;
      object o;
      MarkovOutTree *markov = NULL;

      vfactor_variables.resize(boost::python::len(factor_variables));
      vgeneration_types.resize(boost::python::len(generation_types));

      for(i = 0; i < boost::python::len(factor_variables); i++)
      {
         o = factor_variables[i];
         boost::python::extract<unsigned int> x(o);
         if (x.check())
            vfactor_variables[i] = x();
         else
         {
            status = false;
            error_indices.push_back(i);
         }
      }

      if (!status)
      {
         for(i = 0; i < error_indices.size(); i++)
         error_message << "incorrect type for element " << i
                       << " of argument list: expecting an int ";
         throw_stat_tree_error(error_message);
      }

      else
      {
         for(i = 0; i < boost::python::len(generation_types); i++)
         {
            o = generation_types[i];
            boost::python::extract<int> x(o);
            if (x.check())
               vgeneration_types[i] = x();
            else
            {
               status = false;
               error_indices.push_back(i);
            }
         }
         if (!status)
         {
            for(i = 0; i < error_indices.size(); i++)
            error_message << "incorrect type for element " << i
                          << " of argument list: expecting an int ";
            throw_stat_tree_error(error_message);
         }
      }
      if (status)
      {
         markov = data.markov_out_tree_estimation(error, inb_vomc, inb_ordered_children,
                                                  inb_children_branching, vfactor_variables,
                                                  parent_dependent, vgeneration_types,
                                                  generation_min_inf_bound,
                                                  generation_min_inf_bound_flag);

         if (error.get_nb_error() > 0)
         {
            error_message << error;
            throw_stat_tree_error(error_message);
         }
      }
      return markov;
   }

   // static unsigned int Mot_data_wrapper_nb_state(const MarkovOutTreeData &data)
   // { return data.get_nb_state(); }

   static str Mot_data_wrapper_ascii_write0(const MarkovOutTreeData& data)
   {
      std::stringstream s;
      str res;

      data.line_write(s);
      res = str(s.str());
      return res;
   }

   static MultiPlotSet*
   Mot_data_wrapper_get_plotable(const MarkovOutTreeData& data)
   {
      StatError error;
      MultiPlotSet *ret = data.get_plotable();
      if (ret == NULL)
        stat_tool::wrap_util::throw_error(error);
      return ret;
   }

   static boost::python::list
   Mot_wrapper_get_factor_combinations(const MarkovOutTreeData& data)
   {
     int i, j;
     std::vector<GenerationProcess::index_array> comb = data.get_factor_combinations();
     GenerationProcess::index_array *mfactor; // = new GenerationProcess::index_array();
     boost::python::list ret, *current_factor = NULL;

     for (i = 0; i < comb.size(); i++)
     {
        current_factor = new boost::python::list();
        mfactor = &comb[i];
        for (j = 0; j < mfactor->size(); j++)
           current_factor->append((int)(*mfactor)[j]);
        ret.append(*current_factor);
        delete current_factor;
        current_factor = NULL;
     }
     return ret;
   }

   static boost::python::list
   Mot_data_wrapper_get_distribution_data(const MarkovOutTreeData &data,
                                          boost::python::list &factors)
   {
      bool status = true;
      unsigned int i;
      StatError error;
      ostringstream error_message;
      std::vector<unsigned int> error_indices(0);
      std::vector<DiscreteDistributionData*> res_ptr(0);
      boost::python::list res;
      object o;
      MarkovOutTreeData::index_array fact(boost::python::len(factors), 0);


      for(i = 0; i < boost::python::len(factors); i++)
      {
         o = factors[i];
         boost::python::extract<unsigned int> x(o);
         if (x.check())
            fact[i] = x();
         else
         {
            status = false;
            error_indices.push_back(i);
         }
      }
      if (!status)
      {
         for(i = 0; i < error_indices.size(); i++)
         error_message << "incorrect type for element " << i
                       << " of argument list: expecting an int ";
         throw_stat_tree_error(error_message);
      }
      else
      {
         res_ptr = data.get_distribution_data(error, fact);
         if (error.get_nb_error() > 0)
         {
            error_message << error;
            throw_stat_tree_error(error_message);
         }
         else
            for(i = 0; i < res_ptr.size(); i++)
               res.append(res_ptr[i]);
      }
      return res;
   }

   static boost::python::list
   Mot_data_wrapper_get_joint_distribution_data_pair(const MarkovOutTreeData &data,
                                                     boost::python::list &factors)
   {
      bool status = true;
      unsigned int i, j;
      StatError error;
      ostringstream error_message;
      std::vector<unsigned int> error_indices(0);
      std::vector<std::pair<std::vector<unsigned int>, unsigned int> > res_ptr(0);
      boost::python::list pair;
      boost::python::list res;
      boost::python::list element;
      object o;
      MarkovOutTreeData::index_array fact(boost::python::len(factors), 0);

      for(i = 0; i < boost::python::len(factors); i++)
      {
         o = factors[i];
         boost::python::extract<unsigned int> x(o);
         if (x.check())
            fact[i] = x();
         else
         {
            status = false;
            error_indices.push_back(i);
         }
      }
      if (!status)
      {
         for(i = 0; i < error_indices.size(); i++)
         error_message << "incorrect type for element " << i
                       << " of argument list: expecting an int ";
         throw_stat_tree_error(error_message);
      }
      else
      {
         res_ptr = data.get_joint_distribution(error, fact);

         if (error.get_nb_error() > 0)
         {
            error_message << error;
            throw_stat_tree_error(error_message);
         }
         else
            for(i = 0; i < res_ptr.size(); i++)
            {
               pair = boost::python::list();
               element = boost::python::list();
               // add number of descendants for each state
               for (j = 0; j < res_ptr[i].first.size(); j++)
                  element.append(res_ptr[i].first[j]);
               pair.append(element);
               pair.append(res_ptr[i].second);
               res.append(pair);
            }
      }
      return res;
   }

   static boost::python::list
   Mot_data_wrapper_get_joint_distribution_data_list(const MarkovOutTreeData &data,
                                                     boost::python::list &factors)
   {
      bool status = true;
      unsigned int i, j;
      StatError error;
      ostringstream error_message;
      std::vector<unsigned int> error_indices(0);
      std::vector<std::pair<std::vector<unsigned int>, unsigned int> > res_ptr(0);
      // const boost::python::list empty_list();
      boost::python::list res;
      boost::python::list element;
      object o;
      MarkovOutTreeData::index_array fact(boost::python::len(factors), 0);

      for(i = 0; i < boost::python::len(factors); i++)
      {
         o = factors[i];
         boost::python::extract<unsigned int> x(o);
         if (x.check())
            fact[i] = x();
         else
         {
            status = false;
            error_indices.push_back(i);
         }
      }
      if (!status)
      {
         for(i = 0; i < error_indices.size(); i++)
         error_message << "incorrect type for element " << i
                       << " of argument list: expecting an int ";
         throw_stat_tree_error(error_message);
      }
      else
      {
         res_ptr = data.get_joint_distribution(error, fact);

         if (error.get_nb_error() > 0)
         {
            error_message << error;
            throw_stat_tree_error(error_message);
         }
         else
            for(i = 0; i < res_ptr.size(); i++)
            {
               element = boost::python::list();
               // add number of descendants for each state
               for (j = 0; j < res_ptr[i].first.size(); j++)
                  element.append(res_ptr[i].first[j]);
               // duplicate configuration of number of descendants
               // as much times as occurrences
               for(j = 0; j < res_ptr[i].second; j++)
                  res.append(element);
            }
      }
      return res;
   }

   static boost::python::dict
   Mot_data_wrapper_get_virtual_vertices(const MarkovOutTreeData &data, int itree)
   // validity of arguments to be checked in python
   {
      const MarkovOutTreeData::virtual_vdic * vdic
         = data.get_virtual_vertices_ptr(itree);

      if (vdic == NULL)
      {
         boost::python::dict pdic;
         return pdic;
      }
      else
      {
         // return make_dict<MarkovOutTreeData::virtual_vdic>(*vdic)();
         dict_converter conv(*vdic);
         boost::python::dict pdic = conv.convert();
         return pdic;
      }
   }

   static void
   Mot_data_wrapper_set_virtual_vertex2(const MarkovOutTreeData &data,
                                        int itree, int ivertex)
   // validity of arguments to be checked in python
   { data.set_virtual_vertex(itree, ivertex); }

   static void
   Mot_data_wrapper_set_virtual_vertex3(const MarkovOutTreeData &data,
                                        int itree, int ivertex, bool bvirtual)
   // validity of arguments to be checked in python
   { data.set_virtual_vertex(itree, ivertex, bvirtual); }

   static void
   Mot_data_wrapper_update_markov_reestimation(MarkovOutTreeData &data,
                                               unsigned int inb_vomc,
                                               unsigned int inb_ordered_children,
                                               unsigned int inb_children_branching,
                                               boost::python::list &factor_values,
                                               bool parent_dependent)
   {
      const unsigned int nb_factors = boost::python::len(factor_values);
      int i = 0;
      ostringstream error_message;
      StatError error;
      std::vector<unsigned int> factor_variables(nb_factors);

      for(i = 0; i < nb_factors; i++)
         factor_variables[i] = boost::python::extract<unsigned int>(factor_values[i]);

      data.update_markov_reestimation(error, inb_vomc, inb_ordered_children,
                                      inb_children_branching,
                                      factor_variables, parent_dependent);

      if (error.get_nb_error() > 0)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
   }

}; // WRAP


void class_mot_data()
{
  class_< MarkovOutTreeData, bases<Trees>  >
    ("_MarkovOutTreeData", "MarkovOutTreeData", init< const MarkovOutTreeData&,  bool, bool>())
        .def(init< const Trees&, optional< int > > ())

        .def("Display", WRAP::Mot_data_wrapper_display,
                        "Display(self, bool) -> str \n\n"
                        "Print MarkovOutTreeData definition")

        .def("EstimationMot",
             WRAP::Mot_data_wrapper_markov_out_tree_estimation1,
             return_value_policy< manage_new_object >(),
             "EstimationMot(self) -> _MarkovOutTree \n\n"
             "Estimate a _MarkovOutTreeData from self \n")

        .def("ExtractMarginalGenerationProcess",
              WRAP::Mot_data_wrapper_get_distribution_data,
              "ExtractMarginalGenerationProcess(self, list) -> list[DistributionData] \n\n"
              "Return marginal distributions of generation process")

        .def("ExtractJointGenerationProcessPairs",
              WRAP::Mot_data_wrapper_get_joint_distribution_data_pair,
              "ExtractJointGenerationProcessPairs(self, list) -> list[list, int] \n\n"
              "Return lists of number of descendants in each state "
              "with the number of occurrences")

        .def("ExtractJointGenerationProcessList",
              WRAP::Mot_data_wrapper_get_joint_distribution_data_list,
              "ExtractJointGenerationProcessList(self, list) -> list[list] \n\n"
              "Return lists of number of descendants in each state "
              "(as much copies as occurrences)")

        .def("GetVirtualVertices",
             WRAP::Mot_data_wrapper_get_virtual_vertices,
             "GetVirtualVertices(self, int) \n\n"
             "Get the dictionary of virtual or censored vertices"
             "in a given tree of self \n")

        DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::Mot_data_wrapper_get_plotable,
                                 "return plotable")

        .def("_GetFactorCombinations",
             WRAP::Mot_wrapper_get_factor_combinations,
             "_GetFactorCombinations(self) -> list \n\n"
             "return every possible combination of factors "
             "for the Markov out-tree\n")

        .def("MotSize",
             &MarkovOutTreeData::get_size,
             "MotSize(self, int) \n\n"
             "Return the number of vertices of a given tree"
             "in _MarkovOutTreeData \n")

        .def("NbStates", &MarkovOutTreeData::get_nb_states,
                         "NbStates(self) -> int \n\n"
                         "return the number of states "
                         "of a _MarkovOutTreeData\n")

        .def("SetVirtualVertex",
             WRAP::Mot_data_wrapper_set_virtual_vertex2,
             "SetVirtualVertex(self, int, int) \n\n"
             "Identify a given vertex as virtual or censored in _MarkovOutTreeData \n")

        .def("SetVirtualVertex",
             WRAP::Mot_data_wrapper_set_virtual_vertex3,
             "SetVirtualVertex(self, int, int, bool) \n\n"
             "Identify a given vertex as virtual / censored (or not)"
             "in _MarkovOutTreeData \n")

        .def("StateTrees",
             &MarkovOutTreeData::get_state_markov_out_tree_data,
             return_value_policy< manage_new_object >(),
             "StateTrees(self) -> _MarkovOutTreeData \n\n"
             "Return a _MarkovOutTreeData containing the states as a variable \n")

        .def("UpdateMarkovReestimation",
              WRAP::Mot_data_wrapper_update_markov_reestimation,
             "UpdateMarkovReestimation(self, int, int, int, list, bool) \n\n"
             "Set parameters for visualization of generation processes \n")

        .def("UpdateNbChildrenFrequencyDistribution",
              &MarkovOutTreeData::update_frequency_distributions,
             "UpdateNbChildrenFrequencyDistribution(self) \n\n"
             "Update frequency distribution for the number of children, "
             "Taking into account virtual vertices \n")

        .def(self_ns::str(self)) //__str__
        .def("__str__", WRAP::Mot_data_wrapper_ascii_write0)
    ;

};

#undef WRAP

