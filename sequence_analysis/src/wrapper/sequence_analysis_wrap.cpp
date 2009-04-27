/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: stat_tool_wrap.cpp 6080 2009-03-13 16:11:35Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/



/* WRAPPER Boost.python for sequences class */
#include "export_function.h"
#include "export_tops.h"
#include "export_sequences.h"
#include "export_correlation.h"
#include "export_markovian_sequences.h"
#include "export_nonhomogeneous_markov.h"
#include "export_nonparametric_sequence_process.h"
#include "export_renewal.h"
#include "export_time_events.h"

#include <boost/python.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 103400
#include <boost/python/docstring_options.hpp>
#endif

using namespace boost::python;


// Define python module "_sequence_analysis"
BOOST_PYTHON_MODULE(_sequence_analysis)
{
#if BOOST_VERSION >= 103400
  docstring_options doc_options(true, false);
#endif

  class_function();

  class_sequences();

  class_markovian_sequences();
  class_self_transition();

  class_correlation();

  class_tops();
  class_top_parameters();

  class_nonhomogeneous_markov();
  class_nonhomogeneous_markov_data();

  class_nonparametric_sequence_process();

  class_renewal();    
  class_renewal_data();    

  class_time_events();
}

