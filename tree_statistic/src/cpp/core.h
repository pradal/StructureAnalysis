#ifndef CORE_H
#define CORE_H

#include <exception>
#include <boost/thread/thread.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/tokenizer.hpp>
// #include <eigen2/Eigen/Dense>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>
#include <vector>
#include <map>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <set>

bool SILENT = false;
bool WAIT = false;
unsigned int GRAPH=0;
unsigned int CRITERION=2;

#ifndef CUMUL_THRESHOLD
#define CUMUL_THRESHOLD 0.9999
#endif

#ifndef DSTEP
#define DSTEP 25e-2
#endif

#ifndef MAXITS
#define MAXITS 1e4
#endif

#ifndef MIN
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x,y)  ((x) > (y) ? (x) : (y))
#endif

#ifndef I_INF
#define I_INF std::numeric_limits<unsigned int>::max()
#endif

#ifndef I_NAN
#define I_NAN std::numeric_limits<unsigned int>::quiet_NaN()
#endif

/*#ifndef D_INF
#define D_INF std::numeric_limits<double>::infinity()
#endif*/

#ifndef D_NAN
#define D_NAN std::numeric_limits<double>::quiet_NaN()
#endif

#ifndef EPSILON
#define EPSILON 1e-4
#endif

double generator();
void dbg();
void set_criterion(unsigned int ident);
template<typename T> std::vector< std::vector<T> > t(const std::vector< std::vector<T> >& M);
void set_wait();
void unset_wait();
template<typename T> std::string toString(T v);
template<typename T> std::string toString(const std::vector<T>& v);
template<typename T> std::vector<std::string> InString(const std::vector<T>& v);
template<typename T> std::vector< std::vector<std::string> > InString(const std::vector< std::vector<T> >& M);
std::vector<std::string> split_string(const std::string& string, const std::string& separator);
template<typename T> std::ostream& display(const std::vector<T>& V, const std::vector<std::string>& labels=std::vector<std::string>(), std::ostream& os=std::cout);
template<typename T> std::ostream& display(const std::vector< std::vector<T> >& M, std::ostream& os=std::cout);
template<typename T> std::ostream& display(const std::vector< std::vector<T> >& M, const std::vector<std::string>& clabels, const std::vector<std::string>& rlabels=std::vector<std::string>(), std::ostream& os=std::cout);
template<typename T> std::vector< std::set<T> > compute_categories(const std::vector< std::vector<T> >& values);
template<typename T> std::vector< std::vector<T> > compute_combinations(const std::vector< std::set<T> >& icategories);
std::vector< std::vector<unsigned int> > compute_powerset(const std::vector< std::vector<unsigned int> >& set, unsigned int invariables);
std::vector<double> Eigen2std(const Eigen::VectorXd& v);
std::vector< std::vector<double> > Eigen2std(const Eigen::MatrixXd& M);
template<typename T> Eigen::VectorXd std2Eigen(const std::vector<T>& v);
template<typename T> Eigen::MatrixXd std2Eigen(const std::vector< std::vector<T> >& M);
template<typename T> std::vector< std::vector<T> > v2M(const std::vector<T>& v);

#include "./core.hpp"

#include "./glm.h"
//#include "./univariate_regression.h"
#include "./contingency_table.h"
#endif
