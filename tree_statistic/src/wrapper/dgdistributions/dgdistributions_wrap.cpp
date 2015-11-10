#ifndef _WRAP_CPP
#define _WRAP_CPP
#include "export_dgdistributions.h"
#include <Python.h>
#include <boost/python.hpp>

BOOST_PYTHON_MODULE(_dgdistributions)
{
	class_discrete_graphical_distribution();
	enum_dge_algorithms_ids();
	enum_dge_criterions_ids();
	class_discrete_graphical_parametric();
	class_discrete_graphical_data();
}

#endif
