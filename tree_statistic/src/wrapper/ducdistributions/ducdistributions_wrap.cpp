#ifndef _WRAP_CPP
#define _WRAP_CPP
#include "export_ducdistributions.h"
#include <Python.h>
#include <boost/python.hpp>

BOOST_PYTHON_MODULE(_ducdistributions)
{
	class_discrete_univariate_conditional_distribution();
	enum_ducp_link_ids();
	class_discrete_univariate_conditional_parametric();
	class_discrete_univariate_conditional_data();
}

#endif
