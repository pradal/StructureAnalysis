# -*- coding: utf-8 -*-
# a test for the class mot.MarkovOutTree: constructor and basic methods
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.dmdistributions as dmd
from openalea.stat_tool.plot import gnuplot, set_plotter
from stat_tool.plot import gnuplot, mplotlib, PLOTTER, DISABLE_PLOT
PLOTTER = mplotlib()
#stat_tool.plot.PLOTTER = stat_tool.plot.mplotlib()
stat_tool.plot.DISABLE_PLOT = True

#############################
# Negative Multinomials
#############################

d = dmd.NegativeMultinomial(2, 1, 0.15, [0.3, 0.85])
# print d
# d.plot()
d.Save("tmpfile.dmd", Overwrite=True)
d2 = dmd.DiscreteMultivariateParametric(d)
# d = dmd.DiscreteMultivariateDistribution("tmpfile.dmd")
d3 = dmd.DiscreteMultivariateParametric("tmpfile.dmd")
print d3.Simulate()
# d3.plot()

try:
    import openalea.tree_statistic.trees as trees
    d4 = dmd.DiscreteMultivariateParametric("nmulti_error.dst")
except dmd.StatTreeError, msg:
    print "catch exception:" + str(msg)
else:
    raise dmd.StatError, "Failed to raise StatError"

try:
    db = dmd.NegativeMultinomial(2, 1, 0.1, [0.1, 0.2])
except Exception, e:
    print e
else:
    raise dmd.StatError, "Failed to raise StatError"

d = dmd.NegativeMultinomial(2, 1, 0.1, [0.5, 0.90])
d.plot()
v2 = d.Simulate(100000)
d4 = dmd.NegativeMultinomialEstimation(v2)
# d4.plot()
d4.Display()

vr = stat_tool.Vectors("nmvectors.vec")
d4 = dmd.NegativeMultinomialEstimation(vr)
d5 = dmd.MultivariateParametricEstimation(vr, "NegativeMultinomial")

#############################
# Multinomials
#############################

nb_variable = 4
inf_bound = 3
parameter = 20
bad_probability = [0.4, 0.9, 0.3, 0.2]
probability = [0.4, 0.1, 0.3, 0.2]
try:
    db = dmd.Multinomial(nb_variable, inf_bound,
                        parameter, bad_probability)
except dmd.StatTreeError, msg:
    print "catch exception:" + str(msg)
else:
    raise dmd.StatError, "Failed to raise StatError"

bad_probability[1] = -0.1
try:
    db = dmd.Multinomial(nb_variable, inf_bound,
                        parameter, bad_probability)
except dmd.StatTreeError, msg:
    print "catch exception:" + str(msg)
else:
    raise dmd.StatError, "Failed to raise StatError"

d = dmd.Multinomial(nb_variable, inf_bound,
                    parameter, probability)
d.Display()
v = d.Simulate(1000)
# v.file_ascii_data_write("mvectors.vec", True)

de = dmd.MultivariateParametricEstimation(v, "Multinomial")
de.Display()

dr = dmd.DiscreteMultivariateParametric("multi.dst")
# dr.Display()

try:
    de = dmd.MultivariateParametricEstimation(vr, "Multinomial")
except dmd.StatTreeError, msg:
    print "catch exception:" + str(msg)
else:
    raise dmd.StatError, "Failed to raise StatError"

#############################
# Multivariate Poisson
#############################

dp = dmd.DiscreteMultivariateParametric("mpoisson.dst")
d = dmd.MultivariatePoisson(4, 4, [5., 0.5, 0.4, 0.9, 1.3])
print d
v = d.Simulate(2000)
print v
de1 = dmd.MultivariateParametricEstimation(v, "MultivariatePoisson")
de1.Display()

try:
    dp = dmd.MultivariatePoisson(4, 2, [5., 0.5, 0.4, 0.9])
except dmd.StatTreeError, msg:
    print "catch exception:" + str(msg)
else:
    raise dmd.StatError, "Failed to raise StatError"

try:
    dp = dmd.MultivariatePoisson(4, 2, [-5., 0.5, 0.4, 0.9, 1.3])
except dmd.StatTreeError, msg:
    print "catch exception:" + str(msg)
else:
    raise dmd.StatError, "Failed to raise StatError"

d = dmd.MultivariatePoisson(3, 0, [3., 4., 2., 3.])
d.Display()
v = d.Simulate(3129)
print v
de1 = dmd.MultivariateParametricEstimation(v, "MultivariatePoisson")
de1.Display()

#############################
# Compound multinomials
#############################

d = dmd.CompoundMultinomial(4, [0.2, 0.3, 0.4, 0.1], stat_tool.Binomial(2,10,0.5))

d.Display()
v = d.Simulate(2000)
de = dmd.CompoundMultinomialEstimation(v, "Binomial")
# ED = dmd.MultivariateParametricEstimation(v, "CompoundMultinomial")
de.Display()

file_name = "dmp_write.dst"
d.Save(file_name, "ASCII", True)
d2 = dmd.MultinomialCompoundDiscreteParametric(d)
if d._display() != d2._display():
    equal = False
import os
os.remove(file_name)
