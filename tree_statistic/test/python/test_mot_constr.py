# -*- coding: utf-8 -*-
# a test for the class mot.MarkovOutTree: constructor and basic methods
import sys, os
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.trees as trees, openalea.tree_statistic.mot as mot
import openalea.tree_statistic.hmt as hmt
stat_tool.plot.PLOTTER = stat_tool.plot.mplotlib()
stat_tool.plot.PLOTTER = stat_tool.plot.gnuplot()
import cPickle as pickle


# use MTG
filename = ".." + os.sep + "data" + os.sep + "multivariate_mtg.pkl"
pkl_file = open(filename, 'rb')
g = pickle.load(pkl_file)
pkl_file.close()

from openalea.mtg import treestats
f1 = lambda v: g.node(v).label
f2 = lambda v: g.node(v).Date
f3 = lambda v: g.node(v).cov1
f4 = lambda v: g.node(v).cov2
lf = [f1, f2, f3, f4]
lp = ['State', 'CovI1', 'CovI2', 'CovD1']
T = treestats.extract_trees(g, 0, None, lf, lp)
T = T.SelectIndividual([0])
TInt = T.SelectVariable([0, 1, 2])
TMotInt = mot.MarkovOutTreeData(TInt, StateVariable = 0)

    

TMotInt2 = T.Shift(1, -2)
TMotInt2 = TMotInt2.SelectVariable([0, 1, 2])
TMotInt2 = mot.MarkovOutTreeData(TMotInt2, StateVariable = 0)

filename = "mot3s_2pop.mot"
# filename = "mbp3s.mot"
# filename = "mot3s_2pop_sim.mot"
# m = mot.MarkovOutTree(filename, 10, .9999)
m = mot.MarkovOutTree(filename)
print(m)
print "Number of states", m.NbStates()
m.Display()
m.Display(Detail=2)
# m.plot()

s = m.Simulate([3])
print "Simulated trees: "
print s
print "Display simulated trees: "
s.Display()
# s.plot()
sst = s.StateTrees()

# branching process
filename = "mbp3s.mot"
mbp = mot.MarkovOutTree(filename)
print(mbp)
mbp.Display()

s2 = mbp.Simulate([5, 5])
while (s2.Size() < 100):
    s2 = mbp.Simulate([5, 5])

print "Simulated trees: "
print s2
print "Display simulated trees: "
s2.Display()
sst2 = s2.StateTrees()
# s2.plot()

# s2.ExtractMarginalGenerationProcess([s2])

l = s2.ExtractMarginalGenerationProcess([0])
print l[0].display(Detail=2)

l = s2.ExtractVectorsGenerationProcessByFactor([0], True)
print l

import numpy
A = numpy.array(l)

l = s2.ExtractVectorsGenerationProcessByFactor([0], False)
print l

l = s2._get_factor_combinations()
print l

# branching process with multinomials
filename = "mbp3s_multi.mot"
mbp = mot.MarkovOutTree(filename)
print(mbp)
mbp.Display()

s2 = mbp.Simulate([5, 5])
while (s2.Size() < 100):
    s2 = mbp.Simulate([5, 5])

print "Simulated trees: "
print s2
print "Display simulated trees: "
s2.Display()
sst2 = s2.StateTrees()
# s2.plot()

# s2.ExtractMarginalGenerationProcess([s2])

l = s2.ExtractMarginalGenerationProcess([0])
print l[0].display(Detail=2)




