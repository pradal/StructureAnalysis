## -*- coding: utf-8 -*-
"""Markov out-tree models"""
import string
import openalea.stat_tool.error as check_error
import openalea.tree_statistic._errors as _errors
import openalea.stat_tool as stat_tool, openalea.tree_statistic.trees as trees
import openalea.tree_statistic.hmt._hmt as _hmt
import openalea.tree_statistic.dmdistributions.dmdistributions as dmd

VariableType = stat_tool.VariableTypeBis
StatTreeError = _errors.StatTreeError
CharacteristicType = trees.CharacteristicType
VariableTypeDict = VariableType.values

from openalea.stat_tool import interface
interface.extend_class(_hmt._MarkovOutTree, interface.StatInterface)
interface.extend_class(_hmt._MarkovOutTreeData, interface.StatInterface)

class MarkovOutTree(interface.StatInterface):
    """An implementation of Markov out-tree models."""
    
    def __init__(self, arg, aliasing=False):
        """Initialize a MarkovOutTree by copy or by reading into a file.

        :Usage:
            M = MarkovOutTree("file_name")
            M = MarkovOutTree(MarkovOutTree)."""
        ## Aliasing is used to make an alias between argument and
        ## self -> useful for the connection between MarkovOutTree
        ## and MarkovOutTreeData
        if type(arg)==str:
            # arg is expectedly a file name...
            self.__chmt = _hmt.MotAsciiRead(arg)
            nbvariables = self.NbVariable()
            self._attributes = ["Variable " + str(i) for i in range(nbvariables)]
        elif issubclass(arg.__class__, MarkovOutTree):
            # ... or a MarkovOutTree object...
            if aliasing:
                self.__chmt = arg.__chmt
            else:
               self.__chmt = _hmt._MarkovOutTree(arg.__chmt)
            self._attributes = list(arg._attributes)
        elif issubclass(arg.__class__, _hmt._MarkovOutTree):
            # ... or a _MarkovOutTree object
            if aliasing:
                self.__chmt = arg
            else:
               self.__chmt = _hmt._MarkovOutTree(arg)
            nbvariables = self.NbVariable()
            self._attributes = ["Variable " + str(i) for i in range(nbvariables)]
        else:
            msg = "bad argument type: "+str(type(arg))
            raise TypeError, msg

    def Display(self, ViewPoint = "Data", Detail=1):
        """Display MarkovOutTree object.

        :Usage:
            Display()
            Display(ViewPoint="Data", Detail=2)"""
        print self._display(ViewPoint, Detail)

    def get_plotable(self):
        """Return plotable"""
        return self.__chmt.get_plotable()

    def NbStates(self):
        """Return the number of states of the MarkovOutTree."""
        return self.__chmt.NbStates()

    def NbVariable(self):
        """Return the number of variables of the MarkovOutTree."""
        return self.__chmt.NbVariable()

    def Plot(self, ViewPoint=None,Title=None, variable=None):
        """Graphical output using the Geom 3D viewer for Trees
           or stat_tool.plot for features.

        :Usage:

            Plot()
            Plot(ViewPoint="FirstOccurrenceRoot", variable=0)
            Plot(ViewPoint="Observation", variable=1)
            Plot(ViewPoint="StateProcess")

        :Parameters:

          `ViewPoint` (str) - type of graphical output
                Possible values for ViewPoint:
                    "FirstOccurrenceLeaves"
                    "SojournSize"
                    "NbZones"
                    "NbOccurrences"
                    "Observation"
                    "StateProcess"
          `Title` (str) - main title for figures
          `variable` (int) - variable to use in graph (0 for state variable)"""
        import stat_tool.plot as plot
        if (ViewPoint is None):
            multiplotset = self.get_plotable()
            try:
                plotter = plot.get_plotter()
            except:
                import warnings
                warnings.warn("Cannot use new plotter. Use old style plot.")
                multiplotset = None
            if plot.DISABLE_PLOT:
                return multiplotset
            if(multiplotset is not None):
                plotter.plot(multiplotset, Title)
        else:
            check_error.CheckType([ViewPoint], [str])
            if (Title is None):
                Title = ""
            if (type(variable) == int):
                self._check_variable(variable)
            check_error.CheckType([Title], [str])            
            multiplotset = self.get_plotable()
            viewpoints = [x for x in multiplotset.viewpoint]
            plotvariables = [x for x in multiplotset.variable]
            nb_generation_plot_set = self.__chmt._NbGenerationPlotSet()
            nodata_plot = (len(multiplotset) == nb_generation_plot_set + \
                                self.NbVariable())
            if (ViewPoint.upper() == "STATEPROCESS"):
                check_error.CheckType([variable], [type(None)])
                if (nodata_plot):
                    plotable = [multiplotset[i] for i in range(nb_generation_plot_set)]
                else:
                    plotable = [multiplotset[i] for i in range(len(multiplotset))
                                    if plotvariables[i] == 0]
            elif (ViewPoint.upper() == "OBSERVATION"):
                ftype = CharacteristicType.OBSERVATION
                check_error.CheckType([variable], [int])
                if (nodata_plot):
                    plotable = [multiplotset[nb_generation_plot_set+variable]]
                else:
                    plotable = [multiplotset[i] for i in range(len(multiplotset))
                                    if plotvariables[i] == variable + 1]
            try:
                plotter = plot.get_plotter()
            except:
                import warnings
                warnings.warn("Cannot use new plotter. Use old style plot.")
                plotable = None
            if plot.DISABLE_PLOT:
                return plotable
            if (multiplotset is not None):
                plotter.plot(plotable, Title)

    def Save(self, FileName, Format="ASCII", Overwrite=False):
        """Save MarkovOutTree object into a file.

        :Usage:

            Save(FileName, Format, Overwrite)

        :Parameters:
          `FileName` (str) - name and path of output file
          `Format` (str) - output format: "ASCII" or "SpreadSheet".
          `Overwrite` (bool) - indicates that existing files should not
                be overwritten.
        """
        check_error.CheckType([FileName, Format, Overwrite], [str, str, bool])
        if not Overwrite:
            try:
                f = file(FileName, 'r')
            except IOError:
                f = file(FileName, 'w+')
            else:
                msg = "File " + FileName + " already exists"
                raise IOError, msg
            f.close()
        import string
        if not (string.upper(Format)=="ASCII"
                or string.upper(Format)=="SPREADSHEET"):
            msg = "unknown file format: " + str(Format)
            raise ValueError, msg
        elif (string.upper(Format)=="ASCII"):
            _hmt._MarkovOutTree.FileAsciiWrite(self.__chmt, FileName)
        else:
            _hmt._MarkovOutTree.SpreadsheetWrite(self.__chmt, FileName)

    def Simulate(self, Depths, Factors=None):
        """Generate a sample of trees from self.

        :Usage:

            Simulate(Depths)

        :Parameters:
          `Depths` (list) - List of tree maximal depths
          `Factors` (list) - List of indices of variable considered 
            as factors in the MarkovOutTree model

        :Returns:
            A :ref:`MarkovOutTreeData` instance is returned, with maxmimal depths given
            by `Depths` argument. Maximal depths is used as stopping criterion.

        """
        if (Factors is None):
            cfactors = []
        else:
            cfactors = list(Factors)
        for variable in cfactors:
            self._check_variable(self, variable)
        chmt_data = _hmt._MarkovOutTree.Simulate(self.__chmt, cfactors, Depths)
        return MarkovOutTreeData(chmt_data)

    def _check_variable(self, variable):
        """Check that a given int corresponds to a model variable"""
        check_error.CheckType([variable], [int])
        if ((variable < 0) or (variable >= self.NbVariable())):
            msg = "variable index out of range: "+str(variable)
            raise ValueError(msg)

    def _chmt(self):
        return self.__chmt

    def _display(self, ViewPoint = "Data", Detail=1):
        """Return str for displaying MarkovOutTree object.

        Usage: _display()
               _display(ViewPoint="Data", Detail=2)"""
        check_error.CheckType([ViewPoint], [str])
        check_error.CheckType([Detail], [int])
        if (ViewPoint.upper() != "DATA"):
            msg = 'bad value for argument "ViewPoint": ' + ViewPoint
            raise ValueError, msg
        if ((Detail != 1) and (Detail != 2)):
            msg = 'bad value for argument "Detail": ' + str(Detail)
            raise ValueError, msg
        return _hmt._MarkovOutTree.Display(self.__chmt, Detail==2)

    def __str__(self):
        classstr = str(self.__class__)
        res = classstr + ": " + _hmt._MarkovOutTree.__str__(self.__chmt)
        return res


class MarkovOutTreeData(trees.Trees, interface.StatInterface):
    """A set of trees associated with a Markov out-tree model."""

    def __init__(self, TreesObject, Aliasing = False,
                 AttributeNames = None, StateVariable = None):
        """Initialize a MarkovOutTreeData by copy.

        :Usage:

            M1 = MarkovOutTreeData(MarkovOutTreeData)
            M2 = MarkovOutTreeData(trees.Trees, StateVariable = 0)"""
        ## Aliasing is used to make an alias between TreesObject and
        ## self._Trees__ctrees -> useful for the connection between MarkovOutTree
        ## and MarkovOutTreeData
        ## The attributes names can be provided
        # and the attribute types should have the possibility to be specified
        if issubclass(TreesObject.__class__, MarkovOutTreeData):
            # TreesObject is supposed to be a MarkovOutTreeData object...
            if not(StateVariable is None):
                msg = "Bad type for StateVariable argument: " + str(type(StateVariable))
                msg += ": should be " + str(type(None))
                raise TypeError, msg
            super(MarkovOutTreeData, self).__init__(TreesObject, attribute_names=AttributeNames)
            self._Trees__ctrees = _hmt._MarkovOutTreeData(TreesObject._ctrees(), True, True)
        elif issubclass(TreesObject.__class__, _hmt._MarkovOutTreeData):
            # ... or a _hmt._MarkovOutTreeData object...
            if not(StateVariable is None):
                msg = "Bad type for StateVariable argument: " + str(type(StateVariable))
                msg += ": should be " + str(type(None))
                raise TypeError, msg
            super(MarkovOutTreeData, self).__init__(TreesObject, attribute_names=AttributeNames)
#            if issubclass(markov.__class__, MarkovOutTree):
#                self._Trees__ctrees = _hmt._MarkovOutTreeData(self._ctrees(), markov._chmt())
#            else:
            if Aliasing:
                self._Trees__ctrees = TreesObject
            else:
                self._Trees__ctrees = _hmt._MarkovOutTreeData(TreesObject, True, True)
        elif issubclass(TreesObject.__class__, trees.Trees):
            super(MarkovOutTreeData, self).__init__(TreesObject, attribute_names=AttributeNames)
            if (StateVariable is None):
                self._Trees__ctrees = _hmt._MarkovOutTreeData(TreesObject._ctrees())
                self._Trees__attributes = TreesObject.Attributes()
            else:
                check_error.CheckType([StateVariable], [int])
                TreesObject._valid_cvariable(StateVariable)
                self._Trees__ctrees = _hmt._MarkovOutTreeData(TreesObject._ctrees(), StateVariable)
                attributes = TreesObject.Attributes()
                self._Trees__attributes = [attributes[i] for i in range(len(attributes)) if i != StateVariable]
                types = TreesObject.Types()
                self._Trees__types = [types[i] for i in range(len(types)) if i != StateVariable]
        else:
            msg = "bad type for first argument: "+str(type(TreesObject))
            raise TypeError, msg
        if not(AttributeNames is None):
            self._Trees__attributes = list(AttributeNames)
            if len(self._Trees__attributes) != len(self._Trees__types):
                msg = "Number of variables ("+str(len(self._Trees__types))+\
                      ") and number of attribute names ("+str(len(self._Trees__attributes))+\
                       ") do not match"
                raise Warning, msg
        self.__factor_combinations = None

    def Display(self, ViewPoint = "Data", Detail=1):
        """Display MarkovOutTreeData object.

        :Usage:

            Display()
            Display(ViewPoint="Data", Detail=2)"""
        print self._display(ViewPoint, Detail)

    def Estimate(self, ModelName, NbVariableOrder=0, NbOrderedChildren=0, NbChildrenBranching=0,
                 FactorVariables=[], ParentDependent=True, GenerationTypes=[], GenerationMinInfBound=0,
                 GenerationMinInfBoundEstimation=True):
        """Estimate a Markov out-tree.

        :Usage:

            Estimate("MARKOV_OUT_TREE", NbVariableOrder, NbOrderedChildren, NbChildrenBranching,
                 FactorVariables=[1], ParentDependent=True, GenerationTypes=["NEGATIVE_MULTINOMIAL"], GenerationMinInfBound=0,
                 GenerationMinInfBoundEstimation=True)

        :Parameters:
          * `NbVariableOrder` (int) - number of variable order Markov chains for ordered children
          * `NbOrderedChildren` (int) - number of ordered children among the set of children
          * `NbChildrenBranching` (int) - number of ordered children on which depend the generation processes
          * `FactorVariables` (list of int) - list of discrete variables on which depend the generation processes
          * `ParentDependent` (bool) - True iif the children states depend on the parent states
          * `GenerationTypes` (list of str) - list of distribution families for the generation processes. If empty,
                the family will be estimated. If the length is 1, the same family will be used for every distribution
          * `GenerationMinInfBound` (int) - mininmal inferior bound of the supports of the generation processes
          * `GenerationMinInfBoundEstimation` (bool) - True iif the mininmal inferior bound of the supports
                of the generation processes have to be estimated

        :Returns:
            An object of class :ref:`openalea.tree_statistic.mot.MarkovOutTree` is returned.

        :Examples:

        .. doctest::

            >>> M = T.Estimate("MARKOV_OUT_TREE", 1, 1, 1, [], True, ["NEGATIVE_MULTINOMIAL"], 0, True)
        """
        check_error.CheckType([ModelName, NbVariableOrder, NbOrderedChildren, NbChildrenBranching], [str, int, int, int])
        check_error.CheckType([ParentDependent, GenerationMinInfBound, GenerationMinInfBoundEstimation], [bool, int, bool])
        if (ModelName.upper() != "MARKOV_OUT_TREE"):
            msg = 'bad value for argument "ModelName: "' + str(ModelName)
            msg += '; should be "MARKOV_OUT_TREE"'
            raise ValueError, msg
        for i in range(len(FactorVariables)):
            check_error.CheckType([FactorVariables[i]], [int])
        VGenerationTypes = []
        for i in range(len(GenerationTypes)):
            try:
                gtype = GenerationTypes[i]
            except TypeError:
                raise TypeError, 'Argument "GenerationTypes" must be a list'
            check_error.CheckType([gtype], [str])
            VGenerationTypes += [dmd._DiscreteMultivariateParametricId(gtype)]
        chmt = _hmt._MarkovOutTreeData.EstimationMot(self._Trees__ctrees,
                                                    NbVariableOrder,
                                                    NbOrderedChildren,
                                                    NbChildrenBranching,
                                                    FactorVariables,
                                                    ParentDependent,
                                                    VGenerationTypes,
                                                    GenerationMinInfBound,
                                                    GenerationMinInfBoundEstimation)
        estimated_hmt = MarkovOutTree(chmt)
        estimated_hmt._attributes = self.Attributes()
        return estimated_hmt

    def ExtractMarginalGenerationProcess(self, Factors):
        """Return marginal distributions associated with generation process
        for a given configuration of factors"""
        try:
            l = self._Trees__ctrees.ExtractMarginalGenerationProcess(Factors)
        except StatTreeError, e:
            msg = str(e)
            if (len(msg) > 0):
                msg += "; "
            msg += "Please retry to run ExtractMarginalGenerationProcess"
            self._generation_process_warning(msg)
        else:
            return l

    def ExtractVectorsGenerationProcess(self, Duplicate=True):
        """Return the list of number of descendants in each state
        for every possible configuration of factors.

        :Parameters:
          `Duplicate` (bool) - True if each configuration of descendant must be
                duplicated according to the number of occurrences. False is each
                configuration must be added once, together with the number of occurrences.

        :Returns:
            If Duplicate is True, a list of list is returned. If Duplicate is False,
            a list of list of list is returned ([configuration of vectors, number of occurrences])
        """
        check_error.CheckType([Duplicate], [bool])
        r = {}
        try:
            factors = self._get_factor_combinations()
            for Factor in factors:
                if Duplicate:
                    l = self._Trees__ctrees.ExtractJointGenerationProcessList(Factor)
                else:
                    l = self._Trees__ctrees.ExtractJointGenerationProcessPairs(Factor)
                r[str(Factor)] = l
        except StatTreeError, e:
            msg = str(e)
            if (len(msg) > 0):
                msg += "; "
            msg += "Please retry to run ExtractMarginalGenerationProcess"
            self._generation_process_warning(msg)
        else:
            return r

    def ExtractVectorsGenerationProcessByFactor(self, Factors, Duplicate=True):
        """Return the list of number of descendants in each state
        for a given configuration of factors.


        :Parameters:
          `Factors` (list) - List of values for variables
                on which the number of descendants depends
          `Duplicate` (bool) - True if each configuration of descendant must be
                duplicated according to the number of occurrences. False is each
                configuration must be added once, together with the number of occurrences.

        :Returns:
            If Duplicate is True, a list of list is returned. If Duplicate is False,
            a list of list of list is returned ([configuration of vectors, number of occurrences])
        """
        try:
            factors = self._get_factor_combinations()
            if not(Factors in factors):
                msg = "Invalid factor value: " + str(Factors)
                raise ValueError, msg
            if Duplicate:
                l = self._Trees__ctrees.ExtractJointGenerationProcessList(Factors)
            else:
                l = self._Trees__ctrees.ExtractJointGenerationProcessPairs(Factors)
        except StatTreeError, e:
            msg = str(e)
            if (len(msg) > 0):
                msg += "; "
            msg += "Please retry to run ExtractMarginalGenerationProcess"
            self._generation_process_warning(msg)
        else:
            return l

    def get_plotable(self):
        """Return plotable"""
        multiplotset = self._Trees__ctrees.get_plotable()
        return multiplotset


    def GetCensoredDict(self, TreeId=None):
        """Get dictionary of censored vertices in self.

        :Usage:

            GetCensoredDict(0)

        :Parameters:

          `TreeId` (int) - Identifier of the tree

        :Returns:
             If TreeId is given, the dictionary of censored vertices in tree "TreeId" is returned.
             Otherwise, the list of dictionaries for every tree in self is returned.
             The censored vertices are the leaf vertices whose number of descendants
             is unknown, or the ordered children whose state is unknown,
        """
        if (TreeId is None):
            l = []
            for t in range(self.NbTrees()):
                if trees.Trees._valid_tree(self, t):
                    l.append(self._ctrees().GetVirtualVertices(t))
            return l
        else:
            if trees.Trees._valid_tree(self, TreeId):
                return self._ctrees().GetVirtualVertices(TreeId)

    def NbStates(self):
        """Return the number of states of the state tree process.

        :Usage:

            NbStates()
        """
        return self._ctrees().NbStates()

    def Plot(self, ViewPoint=None, Length=None, BottomDiameter=None,
             DressingFile=None, Title=None, variable=None):
        """Graphical output using the Geom 3D viewer for Trees
           or stat_tool.plot for features.

        :Usage:

            Plot(ViewPoint="Data")
            Plot()
            Plot(ViewPoint="FirstOccurrenceRoot", variable=0)
            Plot(ViewPoint="Observation", variable=1)
            Plot(ViewPoint="StateProcess")

        :Parameters:

          `ViewPoint` (str) - type of graphical output
                Other possible values for ViewPoint:
                    "FirstOccurrenceLeaves"
                    "SojournSize"
                    "NbZones"
                    "NbOccurrences"
                    "Observation"
                    "StateProcess"
          `Length` (function) - function defining vertex length
          `BottomDiameter` (function) - function defining vertex diameter
          `DressingFile` (str) - name of dressing file
          `Title` (str) - main title for figures
          `variable` (int) - variable to use in graph (0 for state variable)"""
        import stat_tool.plot as plot
        if (ViewPoint is None):
            multiplotset = self.get_plotable()
            try:
                plotter = plot.get_plotter()
            except:
                import warnings
                warnings.warn("Cannot use new plotter. Use old style plot.")
                multiplotset = None
            if plot.DISABLE_PLOT:
                return multiplotset
            if(multiplotset is not None):
                plotter.plot(multiplotset, Title)
        else:
            check_error.CheckType([ViewPoint], [str])
            if (ViewPoint.upper() == "Data"):
                # use Plot for Trees
                Color = lambda x: 0
                trees.Trees.Plot(self, ViewPoint, Length, BottomDiameter,
                                Color, DressingFile, Title="", variable=0)
            else:
                if (Title is None):
                    Title = ""
                check_error.CheckType([Length, BottomDiameter, DressingFile],
                                    [type(None), type(None), type(None)])
                check_error.CheckType([Title], [str])
                multiplotset = self.get_plotable()
                viewpoints = [x for x in multiplotset.viewpoint]
                plotvariables = [x for x in multiplotset.variable]
                if (ViewPoint.upper() == "STATEPROCESS"):
                    check_error.CheckType([variable], [type(None)])
                    plotable = [multiplotset[i] for i in range(len(multiplotset))
                                    if plotvariables[i] == 0]
                    if len(plotable) == 0:
                        msg = 'Please retry to run self.Plot(ViewPoint="StateProcess")'
                        self._generation_process_warning(msg)
                elif (ViewPoint.upper() == "OBSERVATION"):
                    ftype = CharacteristicType.OBSERVATION
                    check_error.CheckType([variable], [int])
                    if ((self.Type(variable) == VariableType.INT_VALUE) or
                        (self.Type(variable) == VariableType.REAL_VALUE) or
                        (self.Type(variable) == VariableType.STATE)):
                        cvariable = self._valid_cvariable(variable)+1
                        plotable = [multiplotset[i] for i in range(len(multiplotset))
                                        if plotvariables[i] == cvariable]
                    else:
                        msg = "Bad type for variable " + str(variable)
                        msg += ": " + str(self.Type(variable))
                        raise TypeError, msg
                else:
                    msg = "Bad value for argument ViewPoint: " + str(ViewPoint)
                    raise ValueError, msg
                try:
                    plotter = plot.get_plotter()
                except:
                    import warnings
                    warnings.warn("Cannot use new plotter. Use old style plot.")
                    plotable = None
                if plot.DISABLE_PLOT:
                    return plotable
                if (multiplotset is not None):
                    plotter.plot(plotable, Title)


    def Save(self, FileName, Overwrite=False, VariableNames=None):
        """Save tree into a file as a MTG.

        :Usage:

            Save(FileName, Overwrite, VariableNames)

        :Parameters:
          `FileName` (str) - name and path of output file
          `Overwrite` (bool) - indicates that existing files should not
            be overwritten.
          `VariableNames` (list of str) - variable names 
        """
        if (self.NbVariables() > 0):
            hmt_data = self
        else:
            hmt_data = self.StateTrees()
        super(MarkovOutTreeData, hmt_data).Save(FileName, Overwrite, VariableNames)

    def SetCensoredDict(self, Dict, TreeId=None):
        """Set dictionary of censored vertices in self.
            The censored vertices are the leaf vertices whose number of descendants
            is unknown, or the ordered children whose state is unknown,

        :Usage:

            SetCensoredDict(Dict)
            SetCensoredDict(Dict, TreeId)

        :Parameters:

          `Dict` (dict or list) - List of dictionaries in the case where TreeId = None
            a single dictionary otherwise. Keys must be vertex identifiers, and values
            must be booleans
          `TreeId` (int) - Identifier of the tree whose dictionary is to be set
            (all trees if None)"""
        if not(TreeId is None):
            if trees.Trees._valid_tree(self, TreeId):
                d = []
                for t in range(self.NbTrees()):
                    if (t == TreeId):
                        d.append(Dict)
                    else:
                        d.append({})
        else:
            if (len(Dict) != self.NbTrees()):
                msg = "bad number of dictionaries: " + str(len(Dict))
                msg += " - should be " + str(self.NbTrees())
                raise ValueError, msg
            d = []
            for t in range(self.NbTrees()):
                try:
                    d.append(dict(Dict[t]))
                except ValueError:
                    msg = "Bad type for element " + str(t)
                    msg += " in list of dictionaries: "
                    msg += str(type(Dict[t]))
                    msg += " - dict expected"
                    raise TypeError, msg
        for t in range(self.NbTrees()):
            cd = d[t]
            for k in cd.keys():
                if (self._valid_vid(t, k)):
                    check_error.CheckType([cd[k]], [bool])
                    self._Trees__ctrees.SetVirtualVertex(t, k, cd[k])
        self._Trees__ctrees.UpdateNbChildrenFrequencyDistribution()
        self.__factor_combinations = None

    def Size(self, TreeId=None):
        """Return the number of vertices of a given tree."""
        if (TreeId is None):
            return trees.Trees.Size(self)
        elif trees.Trees._valid_tree(self, TreeId):
            return self._Trees__ctrees.MotSize(TreeId)

    def StateTrees(self):
        """Add State variable to the list of output variables
        
        :Returns:
             A MarkovOutTreeData is returned. If self contained to state tree,
             a copy of self is returned. Otherwise, a copy of self with a new
             variable corresponding to the state variable is returned."""
        cp = self._Trees__ctrees.StateTrees()
        attributes = self.Attributes()
        if (cp.NbInt() == self.NbInt()+1):
            # self contained a state tree
            attributes = ["StateVariable"] + attributes
        # otherwise self contained no state tree and a mere copy is returned
        return MarkovOutTreeData(cp, AttributeNames = attributes, Aliasing = True)

    def Tree(self, TreeId):
        """Return the tree corresponding to the given identifier."""
        if (self.NbVariables() == 0):
            res = self.StateTrees().Tree(TreeId)
        else:
            res = trees.Trees.Tree(self, TreeId)
        return res

    def _display(self, ViewPoint = "Data", Detail=1):
        """Return str for displaying MarkovOutTreeData object.

        Usage: _display()
               _display(ViewPoint="Data", Detail=2)"""
        check_error.CheckType([ViewPoint], [str])
        check_error.CheckType([Detail], [int])
        if (ViewPoint.upper() != "DATA"):
            msg = 'bad value for argument "ViewPoint": ' + ViewPoint
            raise ValueError, msg
        if ((Detail != 1) and (Detail != 2)):
            msg = 'bad value for argument "Detail": ' + str(Detail)
            raise ValueError, msg
        return _hmt._MarkovOutTreeData.Display(self._Trees__ctrees, Detail==2)

    def _ctrees(self):
        return self._Trees__ctrees

    def _generation_process_warning(self, post_message):
        """Return a warning in the case where mot_reestimation has not been computed.
        Warning is updated with post_message
        Compute default values"""
        msg = "Characteristics for generation process have not bee computed.\n"
        msg = "Most likely, this dataset is not connected to any model.\n"
        msg += "Default values will be used to define the number of "
        msg += "ordered and unordered children.\n"
        msg += "To change default parameters, please use "
        msg += "the _UpdateMarkovReestimation method.\n"
        if (len(post_message) > 0):
            msg += str(post_message) + "\n"
        self._UpdateMarkovReestimation()
        raise Warning, msg
    
    def _get_factor_combinations(self):
        """Return every possible combination of factors"""
        if self.__factor_combinations is None:
            self.__factor_combinations = self._Trees__ctrees._GetFactorCombinations()
        return self.__factor_combinations

    def _UpdateMarkovReestimation(self, NbVariableOrder=0, NbOrderedChildren=0, NbChildrenBranching=0,
                                  FactorVariables=[], ParentDependent=True):
        """Set parameters for visualization of generation processes

        :Usage:

            _UpdateMarkovReestimation(NbVariableOrder, NbOrderedChildren, NbChildrenBranching,
                 FactorVariables=[1], ParentDependent=True)

        :Parameters:

          `NbVariableOrder` (int) - number of variable order Markov chains for ordered children
          `NbOrderedChildren` (int) - number of ordered children among the set of children
          `NbChildrenBranching` (int) - number of ordered children on which depend the generation processes
          `FactorVariables` (list of int) - list of discrete variables on which depend the generation processes
          `ParentDependent` (bool) - True iif the children states depend on the parent states
        """
        check_error.CheckType([NbVariableOrder, NbOrderedChildren, NbChildrenBranching], [int, int, int])
        check_error.CheckType([ParentDependent], [bool])
        for i in range(len(FactorVariables)):
            check_error.CheckType([FactorVariables[i]], [int])
            msg = 'bad value for factor variable ' + str(i) + ': ' + str(FactorVariables[i])
            assert(FactorVariables[i] < self.NbVariables()), msg
        _hmt._MarkovOutTreeData.UpdateMarkovReestimation(self._Trees__ctrees, NbVariableOrder, NbOrderedChildren,
                                                         NbChildrenBranching, FactorVariables, ParentDependent)
        self.__factor_combinations = None

    def __str__(self):
        classstr = str(self.__class__)
        res = classstr + ": " + _hmt._MarkovOutTreeData.__str__(self._Trees__ctrees)
        return res
