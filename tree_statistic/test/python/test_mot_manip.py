# -*- coding: utf-8 -*-
import sys, os
import openalea.tree_statistic.mot as mot
import openalea.tree_statistic.trees as trees
import openalea.stat_tool as stat_tool
stat_tool.plot.PLOTTER = stat_tool.plot.mplotlib()
stat_tool.plot.PLOTTER = stat_tool.plot.gnuplot()
import cPickle as pickle

DISABLE_PLOT = False

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
# functions to count number of descendants
def NbDescendantStates(g, v, NbStates, PropertyName):
    """Count the number of descendants in each state"""
    c = g.children(v)
    prop = g.property(PropertyName)
    res = []
    for s in range(NbStates):
        res += [len([u for u in c if not(prop[u] is None) and prop[u] == s])]
    return res

def NbDescendantCount(g, NbStates, PropertyName):
    """Count the number of configurations of descendants in each state"""
    d = {}
    prop = g.property(PropertyName)
    for i in range(NbStates):
        d[i] = {}
    for v in g.vertices():
        if not(prop[v] is None):
            s = prop[v]
            conf = str(NbDescendantStates(g, v, NbStates, PropertyName))
            if d[s].has_key(conf):
                d[s][conf] += 1
            else:
                d[s][conf] = 1
    return d
# dependence wrt parent only
TMotInt._UpdateMarkovReestimation(ParentDependent=True)
l = {}
for f in TMotInt._get_factor_combinations():
    l[str(f)] = TMotInt.ExtractVectorsGenerationProcessByFactor(f, False)

nbcount = 0
for i in l.items():
    for v in i[1]:
        nbcount += v[1]

assert(nbcount == TInt.Tree(0).Size())

# dependence wrt parent and one brother child
# vertices with less than 1 child are ignored
censored = {}
from openalea.tree_statistic.trees import etrees
for t in range(TMotInt.NbTrees()):
    Tr = etrees.Tree(TMotInt.Tree(t))
    for v in Tr.Preorder():
        censored[v] = (Tr.NbChildren(v) < 1)

def DiscreteDistributionDataToArray(h):
    import numpy
    l = 0
    try:
        while True:
            h[l]
            l += 1
    except Exception:
        pass
    
    res = []
    for i in range(l):
        res += [h[i]]

    res = numpy.array(res)
    return res

res = DiscreteDistributionDataToArray(TMotInt.ExtractHistogram("NbChildren"))
res[0]
TMotInt.SetCensoredDict([censored])

res = DiscreteDistributionDataToArray(TMotInt.ExtractHistogram("NbChildren"))
assert(res[0] == 0)

TMotInt._UpdateMarkovReestimation(NbVariableOrder=1, NbOrderedChildren=1, NbChildrenBranching=1, ParentDependent=True)
TMotInt._get_factor_combinations()

l = {}
for f in TMotInt._get_factor_combinations():
    l[str(f)] = TMotInt.ExtractVectorsGenerationProcessByFactor(f, False)

nbcount = 0
for i in l.items():
    for v in i[1]:
        nbcount += v[1]

assert(res.sum() == nbcount)

TMotInt2 = T.Shift(1, -2)
TMotInt2 = TMotInt2.SelectVariable([0, 1, 2])
TMotInt2 = mot.MarkovOutTreeData(TMotInt2, StateVariable = 0)
TMotInt2.SetCensoredDict([censored])
TMotInt2._UpdateMarkovReestimation(FactorVariables=[0,1], ParentDependent=True)
TMotInt2._get_factor_combinations()
l = {}
for f in TMotInt2._get_factor_combinations():
    l[str(f)] = TMotInt2.ExtractVectorsGenerationProcessByFactor(f, False)
    print f

nbcount = 0
for i in l.items():
    for v in i[1]:
        nbcount += v[1]

assert(res.sum() == nbcount)

mbp_name1 = "mbp3s.mot"
MBP = mot.MarkovOutTree(mbp_name1)
MBPData = MBP.Simulate([5, 5])
MBPData = MBP.Simulate([15, 15])
MBPData = MBP.Simulate([16, 16])

cd = MBPData.GetCensoredDict(0)

assert(MBPData)
assert(MBPData.NbTrees()==2)

if not(DISABLE_PLOT):
    MBP.plot()

if not(DISABLE_PLOT):
    MBPData.plot()

MBPData2 = mot.MarkovOutTreeData(MBPData)
assert(MBPData2)
assert(MBPData2._display() == MBPData._display())

SMBPData = MBPData.StateTrees()
# print SMBPData.Tree(1)
# SMBPData.Save("mbp3s.mtg")

assert(SMBPData)
assert(SMBPData.NbTrees() == MBPData.NbTrees())

l = MBPData.ExtractMarginalGenerationProcess([0])

msg = "Bad number of states: "
msg += "model: " + str(MBP.NbStates())
nbstate = (SMBPData._max(0) + 1)
msg += " / data: " + str(nbstate)
assert (MBP.NbStates() == nbstate), msg
assert(len(l) == SMBPData._max(0) + 1)
    
##############################################
# MBP with output processes
##############################################


mbp_name2 = "mbp3s_2pop.mot"
MBP2 = mot.MarkovOutTree(mbp_name2)
MBPData2 = MBP2.Simulate([5, 5])
assert(MBPData2.ExtractMarginalGenerationProcess([0]))
# MBPData2 = MBP2.Simulate([15, 15])
MBP2.plot()
MBPData2.plot()
# MBPData2.Plot(ViewPoint="StateProcess")
MBP2.Plot(ViewPoint="Observation", variable = 1)

# build MarkovOutTreeData from trees.Trees and variable
TF = trees.Trees("sample_mtg_forest.txt")
T = mot.MarkovOutTreeData(TF)
T1 = mot.MarkovOutTreeData(TF, StateVariable=3)

# Check on estimation

mbp_name2 = "mbp3s_2pop.mot"
MBP2 = mot.MarkovOutTree(mbp_name2)
MBPData2 = MBP2.Simulate([10, 10])
while (MBPData2.Size() < 300):
    MBPData2 = MBP2.Simulate([10, 10])

# add attribute names
MBPData2 = mot.MarkovOutTreeData(MBPData2, AttributeNames = ["MV1", "MV2"])

print "Data set:"
MBPData2.Display()
sst2 = MBPData2.StateTrees()
t2 = trees.Trees(sst2)
# build state trees
t3 = t2.SelectVariable([0])
MBP2S = mot.MarkovOutTreeData(t3)
assert(MBP2S.Attributes()[0] == "StateVariable")
MBP2SB = mot.MarkovOutTreeData(sst2, 0)
vd = MBP2SB.GetCensoredDict()
MBP2SC = mot.MarkovOutTreeData(t2, StateVariable=0)
assert(len(MBP2SC.Attributes()) == 2)
# build state MarkovOutTreeData
t3s = mot.MarkovOutTreeData(t3, StateVariable=0)
# command below fails
# st3s = t3s.StateTrees()

MBP2SC.SetCensoredDict(vd)
assert(MBP2SC.GetCensoredDict() == MBPData2.GetCensoredDict())
sst2.Tree(0)._display(0) == MBP2SB.Tree(0)._display(0)
sst2.Tree(0)._display(0) == MBPData2.Tree(0)._display(0)
EMBP1 = MBPData2.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["NegativeMultinomial"], 0, True)
EMBP2 = MBP2SB.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["NegativeMultinomial"], 0, True)
EMBP3 = MBP2SC.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, ["NegativeMultinomial"], 0, True)

assert(EMBP1._display() == EMBP3._display())

assert(MBPData2.get_plotable())
# _UpdateMarkovReestimation() should be called by plot
MBP2SC._UpdateMarkovReestimation()
assert(MBP2SC.get_plotable())

MBP2SC.Plot(ViewPoint="StateProcess")
MBP2SC.Plot(ViewPoint="Observation", variable = 0)
MBP2SC.Plot(ViewPoint="Observation", variable = 1)

t2s = t2.SelectVariable([0])
MBP2SD = mot.MarkovOutTreeData(t2s, StateVariable=0)
MBP2SD._UpdateMarkovReestimation()
assert(len(MBP2SD.ExtractMarginalGenerationProcess([0])) == MBP2SD.NbStates())
