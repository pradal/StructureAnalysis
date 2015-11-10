/****************************************************************
 *
 *  Test of the MultivariateDiscreteDistribution methods
 *  as defined in multivariate_distributions.h
 */

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"

// required by "tree_statistic/multivariate_distribution.h"
#include <boost/math/distributions/binomial.hpp>
#include "stat_tool/stat_label.h"

#include "tree_statistic/tree_labels.h"
#include "statiskit/core/data/marginal/multivariate.h"

#include "tree_statistic/multivariate_distribution.h"
#include <iostream>
#include <string>
#include <cstdlib>

#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>

#include "assert.h"

#include <boost/math/special_functions/gamma.hpp>

#include "tree_statistic/core.h"

using namespace stat_tool;
using namespace Stat_trees;

int main(void) {
    unsigned int i,j,nb_variable = 3, inf_bound = 0, sample_size = 1000;
    std::vector<unsigned int> msup_bound(nb_variable,2);
    std::vector< std::vector<unsigned int> > sample(sample_size);
    std::vector<double> mparameter(nb_variable+1), mprobability(nb_variable);
    mparameter[0] = 1;
    mparameter[1] = 0.1;
    mparameter[2] = 0.2;
    mprobability[0] = 0.3;
    mprobability[1] = 0.4;
    mprobability[2] = 0.2;
    DiscreteParametric* poisson = new DiscreteParametric(POISSON, 0, I_DEFAULT, 0.3,D_DEFAULT);
    DiscreteMultivariateParametric* mpoisson = new DiscreteMultivariateParametric(nb_variable, MMULTINOMIAL, inf_bound, msup_bound, mparameter, mprobability);
    for(i = 0; i < sample.size(); i++) {
        sample[i] = mpoisson->simulation();
        sample[i].insert(sample[i].end(),poisson->simulation());
        for(j = 0; j < sample[i].size(); j++)
            std::cout << sample[i][j] << " ";
        std::cout << std::endl;
    }
    ContingencyTable<unsigned int> *CT = new ContingencyTable<unsigned int>(sample);
    CT->Parsimoniest();


    return 0;
}
