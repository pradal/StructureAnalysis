#ifndef _WRAP_CPP
#define _WRAP_CPP
#include "export_dmcdistributions.h"
#include <Python.h>
#include <boost/python.hpp>

BOOST_PYTHON_MODULE(_dmcdistributions)
{
	class_discrete_multivariate_conditional_distribution();
	enum_dmcp_distributions_ids();
	class_discrete_multivariate_conditional_parametric();
	class_discrete_multivariate_conditional_data();
}

#endif
