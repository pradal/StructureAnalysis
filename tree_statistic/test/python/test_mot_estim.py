# -*- coding: utf-8 -*-
import openalea.tree_statistic.mot as mot
import openalea.stat_tool as stat_tool
from stat_tool.plot import gnuplot
stat_tool.plot.PLOTTER = gnuplot()

ENABLE_PLOT = False

mbp_name1 = "mbp3s.mot"
MBP = mot.MarkovOutTree(mbp_name1)
MBPData = MBP.Simulate([5, 5])
while (MBPData.Size() < 100):
    MBPData = MBP.Simulate([5, 5])


SMBPData = MBPData.StateTrees()
# print SMBPData.Tree(1)
# SMBPData.Save("mbp3s.mtg", Overwrite=True)
# MBPData.Save("mbp3s2.mtg", Overwrite=True)

# Estimation
print "Data set:"
MBPData.Display()

print "True model:"
MBP.Display()

EMBP = MBPData.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["NegativeMultinomial"], 0, True)
EMBP = SMBPData.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["NegativeMultinomial"], 0, True)

print "Estimated model:"
EMBP.Display()

##############################################
# multinomial distributions
##############################################
mbp_name1 = "mbp3s_multi.mot"
MBP = mot.MarkovOutTree(mbp_name1)
MBPData = MBP.Simulate([5, 5])
while (MBPData.Size() < 10000):
    MBPData = MBP.Simulate([5, 5])


SMBPData = MBPData.StateTrees()

# Estimation
print "Data set:"
MBPData.Display()

print "True model:"
MBP.Display()

EMBP = MBPData.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["Multinomial"], 0, True)

print "Estimated model:"
EMBP.Display()

# Estimate multivariate Poisson distribution on multinomials
try:
    EMBP = MBPData.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["MultivariatePoisson"], 0, True)
except mot.StatError, msg:
    print msg


##############################################
# multivariate Poisson distributions
##############################################
mbp_name1 = "mbp3s_mpoisson.mot"
MBP = mot.MarkovOutTree(mbp_name1)
MBPData = MBP.Simulate([5, 5])
while (MBPData.Size() < 10000):
    MBPData = MBP.Simulate([5, 5])


SMBPData = MBPData.StateTrees()

# Estimation
print "Data set:"
MBPData.Display()

print "True model:"
MBP.Display()

EMBP = MBPData.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["MultivariatePoisson"], 0, True)

print "Estimated model:"
EMBP.Display()

# EMBP.plot()

##############################################
# Estimation when output processes are present
##############################################


mbp_name2 = "mbp3s_2pop.mot"
MBP2 = mot.MarkovOutTree(mbp_name2)
MBPData2 = MBP2.Simulate([10, 10, 10])
while (MBPData2.Size() < 1000):
    MBPData2 = MBP2.Simulate([10, 10, 10])


# Estimation
print "Data set:"
MBPData2.Display()

print "True model:"
MBP2.Display()

EMBP2 = MBPData2.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["NegativeMultinomial"], 0, True)

print "Estimated model:"
EMBP2.Display()

# EMBP2.Plot(ViewPoint="Observation", variable = 0)
MBP2.Plot(ViewPoint="Observation", variable = 0)
MBP2.Plot(ViewPoint="StateProcess")
