# -*- coding: utf-8 -*-
"""Markov out-tree models"""
import string
import openalea.stat_tool.error as check_error
import openalea.stat_tool as stat_tool
import _dmdistributions
# StatError = _dmdistributions.StatErrorClass
import openalea.tree_statistic._errors as _errors
# StatError = _errors.StatError
StatTreeError = _errors.StatTreeError

PMDistributionIds = _dmdistributions.PMDistributionIds

from openalea.stat_tool import interface
interface.extend_class(_dmdistributions._DiscreteMultivariateParametric, interface.StatInterface)

class DiscreteMultivariateDistribution(_dmdistributions._DiscreteMultivariateDistribution):
    """Discrete Multivariate Distributions."""

    def __init__(self, arg):
        """Initialize a DiscreteMultivariateDistribution by copy, 
        or by reading into a file.

        :Usage:

            M = DiscreteMultivariateDistribution("Name")
            M = DiscreteMultivariateDistribution(DiscreteMultivariateDistribution).
        """
        _dmdistributions._DiscreteMultivariateDistribution.__init__(self, arg)

class DiscreteMultivariateParametric(_dmdistributions._DiscreteMultivariateParametric):
    """Parametric Discrete Multivariate Distributions."""

    def __init__(self, utype, *args):
        """Initialize a DiscreteMultivariateParametric by copy, using its name and parameters,
        or by reading into a file.

        :Usage:

            M = DiscreteMultivariateParametric("file_name")
            M = DiscreteMultivariateParametric(DiscreteMultivariateParametric).
        """
        if ((type(utype)==str)
            or (issubclass(utype.__class__, _dmdistributions._DiscreteMultivariateParametric))):
            # utype is expectedly a file name or a _DiscreteMultivariateParametric
            _dmdistributions._DiscreteMultivariateParametric.__init__(self, utype)

    def Display(self, Detail=1):
        """Display DiscreteParametricDistribution  object.

        :Usage:

            Display()
            Display(Detail=2)
        """
        print self._display(Detail)

    def Save(self, FileName, Format="ASCII", Overwrite=False):
        """Save DiscreteMultivariateParametric object into a file.

        :Usage:

            Save(FileName, Format, Overwrite)

        :Parameters:
          `Filename` (str) - name and path of output file
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
            _dmdistributions._DiscreteMultivariateParametric.FileAsciiWrite(self, FileName)
        else:
            _dmdistributions._DiscreteMultivariateParametric.SpreadsheetWrite(self, FileName)

    def Simulate(self, SampleSize=1):
        """Generate a vector from self.

        :Usage:  Simulate(SampleSize)

        :Parameters:
          `SampleSize` (int) - Number of vectors to simulate

        :Returns:
            A list of lists with size SampleSize is simulated if SampleSize > 1
            a random list of int is returned otherwise

        """
        check_error.CheckType([SampleSize], [int])
        if (SampleSize == 1):
            return _dmdistributions._DiscreteMultivariateParametric.Simulate(self)
        elif (SampleSize > 1):
            # return _dmdistributions._DiscreteMultivariateParametric.Simulate(self, SampleSize)
            return _dmdistributions._DiscreteMultivariateParametric.Simulate(self, SampleSize)
        else:
            return []

    def _display(self, Detail=1):
        """Return str for displaying DiscreteParametricDistribution  object.

        Usage: _display()
               _display(Detail=2)"""
        check_error.CheckType([Detail], [int])
        if ((Detail != 1) and (Detail != 2)):
            msg = 'bad value for argument "Detail": ' + str(Detail)
            raise ValueError, msg
        return _dmdistributions._DiscreteMultivariateParametric.Display(self, Detail==2)

    def __str__(self):
        classstr = str(self.__class__)
        res = classstr + ": " +  _dmdistributions._DiscreteMultivariateParametric.__str__(self)
        return res

class MultinomialCompoundDiscreteParametric(_dmdistributions._MultinomialCompoundDiscreteParametric):
    """Compound multinomial Distributions."""

    def __init__(self, utype, *args, **kargs):
        """Initialize a MultinomialCompoundDiscreteParametric by copy, using its name and parameters,
        or by reading into a file.

        :Usage:

            M = MultinomialCompoundDiscreteParametric("file_name")
            M = MultinomialCompoundDiscreteParametric(MultinomialCompoundDiscreteParametric)
            M = MultinomialCompoundDiscreteParametric(NbVariable, InfBound, Probability, CompoundDistribution).
        """
        if ((type(utype)==str)
            or (issubclass(utype.__class__, _dmdistributions._MultinomialCompoundDiscreteParametric))):
            _dmdistributions._MultinomialCompoundDiscreteParametric.__init__(self, utype)
        else:
            # check names, types and values of arguments
            _dmdistributions._MultinomialCompoundDiscreteParametric.__init__(NbVariable, InfBound, Probability, CompoundDistribution)

    def Display(self, Detail=1):
        """Display DiscreteParametricDistribution  object.

        :Usage:

            Display()
            Display(Detail=2)
        """
        print self._display(Detail)

    def Save(self, FileName, Format="ASCII", Overwrite=False):
        """Save MultinomialCompoundDiscreteParametric object into a file.

        :Usage:

            Save(FileName, Format, Overwrite)

        :Parameters:
            `Filename` (str) - name and path of output file
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
        if not (string.upper(Format)=="ASCII" or string.upper(Format)=="SPREADSHEET"):
            msg = "unknown file format: " + str(Format)
            raise ValueError, msg
        elif (string.upper(Format)=="ASCII"):
            _dmdistributions._MultinomialCompoundDiscreteParametric.FileAsciiWrite(self, FileName)
        else:
            _dmdistributions._MultinomialCompoundDiscreteParametric.SpreadsheetWrite(self, FileName)


    def Simulate(self, SampleSize=1):
        """Generate a vector from self.

        :Usage:

            Simulate(SampleSize)

        :Parameters:
            `SampleSize` (int) - Number of vectors to simulate

        :Returns:
            A list of lists with size SampleSize is simulated if SampleSize > 1
            a random list of int is returned otherwise

        """
        check_error.CheckType([SampleSize], [int])
        if (SampleSize == 1):
            return _dmdistributions._MultinomialCompoundDiscreteParametric.Simulate(self)
        elif (SampleSize > 1):
            return _dmdistributions._MultinomialCompoundDiscreteParametric.Simulate(self, SampleSize)
        else:
            return []

    def _display(self, Detail=1):
        """Return str for displaying DiscreteParametricDistribution  object.

        :Usage:

            _display()
            _display(Detail=2)
        """
        check_error.CheckType([Detail], [int])
        if ((Detail != 1) and (Detail != 2)):
            msg = 'bad value for argument "Detail": ' + str(Detail)
            raise ValueError, msg
        return _dmdistributions._MultinomialCompoundDiscreteParametric.Display(self, Detail==2)

    def __str__(self):
        classstr = str(self.__class__)
        res = classstr + ": " +  _dmdistributions._MultinomialCompoundDiscreteParametric.__str__(self)
        return res

def NegativeMultinomial(NbVariable, InfBound, Parameter, Probability):
    """Construction of a Negative Multinomial distribution

    :Parameters:
        `NbVariable` (int) - Number of variables
        `InfBound` (int) - minimal possible value
        `Parameter` (float) - distribution parameter
        `Probability` (list of float) - probabilities of 'success'

    :Returns:
        A DiscreteMultivariateParametric object with NegativeMultinomial Id
        is returned
    """
    check_error.CheckType([NbVariable, InfBound, Parameter], [int, int, float])
    for i in Probability:
        check_error.CheckType([i], [float])
    try:
        res = _dmdistributions._DiscreteMultivariateParametric(NbVariable,
                                                            PMDistributionIds.MNEGATIVE_MULTINOMIAL,
                                                            InfBound, [], [Parameter], Probability)
    except _errors.StatTreeError, error:
        raise _errors.StatTreeError(error)
    else:
        return DiscreteMultivariateParametric(res)

def NegativeMultinomialEstimation(Vectors, InfBound=0, InfBoundFlag=True):
    """Estimate a Negative Multinomial distribution from vectors"""
    res = _dmdistributions.NegativeMultinomialEstimation(Vectors, InfBound, InfBoundFlag)
    return DiscreteMultivariateParametric(res)

def MultivariatePoisson(NbVariable, InfBound, Parameter):
    """Construction of a Multivariate Poisson distribution

    :Parameters:
        `NbVariable` (int) - Number of variables
        `InfBound` (int) - minimal possible value
        `Parameter` (list of float) - distribution parameters

    :Returns:
        A DiscreteMultivariateParametric object with MultivariatePoisson Id is returned
    """
    check_error.CheckType([NbVariable, InfBound], [int, int])
    for i in Parameter:
        check_error.CheckType([i],[float])
    res = _dmdistributions._DiscreteMultivariateParametric(NbVariable,
                                                           PMDistributionIds.MPOISSON,
                                                           InfBound, [], Parameter, [])
    return DiscreteMultivariateParametric(res)


def Multinomial(NbVariable, InfBound, Parameter, Probability):
    """Construction of a Multinomial distribution

    :Parameters:
        `NbVariable` (int) - Number of variables
        `InfBound` (int) - minimal possible value
        `Parameter` (int) - distribution parameter
        `Probability` (list of float) - distribution probabilities

    :Returns:
        A DiscreteMultivariateParametric object with Multinomial Id is returned
    """
    check_error.CheckType([NbVariable, InfBound, Parameter], [int, int, int])
    for i in Probability:
        check_error.CheckType([i],[float])
    res = _dmdistributions._DiscreteMultivariateParametric(NbVariable,
                                                           PMDistributionIds.MMULTINOMIAL,
                                                           InfBound, [], [Parameter], Probability)
    return DiscreteMultivariateParametric(res)

def CompoundMultinomial(NbVariable, Probability, CompoundDistribution):
    """Construction of a Compound Multinomial distribution

    :Usage:

        CompoundMultinomial(4, 2, [0.2, 0.3, 0.4, 0.1], stat_tool.Binomial(0,10,0.5))

    :Parameters:
       `Identifier` (int) - Code for the compound distribution
       `Probability` (list of float) - multinomial distribution probabilities
       `CompoundDistribution` (stat_tool._stat_tool._DiscreteParametric) - distributions
           of the sum of random vector

    :Returns:
        A MultinomialCompoundDiscreteParametric instance is returned
    """
    check_error.CheckType([NbVariable], [int])
    for p in Probability:
        check_error.CheckType([p], [float])
    if not(issubclass(CompoundDistribution.__class__,
            stat_tool._stat_tool._DiscreteParametricModel)):
        msg = "Bad type for argument 'CompoundDistribution':"
        msg += str(type(CompoundDistribution))
        msg += " - should be: openalea.stat_tool._stat_tool._DiscreteParametricModel"
        raise TypeError, msg
    res = _dmdistributions._MultinomialCompoundDiscreteParametric(NbVariable, Probability, CompoundDistribution)
    return MultinomialCompoundDiscreteParametric(res)
    

def MultivariateParametricEstimation(Vectors, Identifier, InfBound=None):
    """Estimate a parametric discrete multivariate distribution from vectors

    :Parameters:
        `Vectors` (stat_tool.Vectors) - Vectors
        `Identifier` (str) - name of the family of parametric discrete multivariate distributions
        `InfBound` (int) - inferior bound of the support of the distribution (if known)

    :Returns:
        A DiscreteMultivariateParametric object is returned

    Usage: MultivariateParametricEstimation(Vectors, "NegativeMultinomial")
           MultivariateParametricEstimation(Vectors, "MultivariatePoisson")
           MultivariateParametricEstimation(Vectors, "Multinomial")
    """
    if (InfBound is None):
        iInfBound = 0
        InfBoundFlag = True
    else:
        iInfBound = InfBound
        InfBoundFlag = False
    check_error.CheckType([Identifier, iInfBound], [str, int])
    iIdent = _DiscreteMultivariateParametricId(Identifier)
    if not(issubclass(Vectors.__class__, stat_tool._stat_tool._Vectors)):
        msg = "Bad type for Vectors: " + str(Vectors.__class__)
        raise TypeError, msg
    res = _dmdistributions.MultivariateParametricEstimation(Vectors, iIdent, iInfBound, InfBoundFlag)
    return DiscreteMultivariateParametric(res)

def OrdinalAssociationModel(Vectors):
    res = _dmdistributions.OrdinalAssociationModel(Vectors)
    return res

# should be included in MultivariateParametricEstimation
def CompoundMultinomialEstimation(Vectors, Identifier, InfBound=None):
    """Estimate a parametric discrete compound multinomial distribution from vectors

        :Parameters:
            `Vectors` (stat_tool.Vectors) - Vectors
            `Identifier` (str) - name of the compound distribution family
            `InfBound` (int) - inferior bound of the support
                of the compounding distribution (if known)

        :Returns:
            A MultinomialCompoundDiscreteParametric object is returned

        Usage:
            CompoundMultinomialEstimation(Vectors, "BINOMIAL")
						CompoundMultinomialEstimation(Vectors, "DEFAULT")
  """
    if (InfBound is None):
        iInfBound = 0
        InfBoundFlag = True
    else:
        iInfBound = InfBound
        InfBoundFlag = False
    check_error.CheckType([iInfBound], [int])
    if not(issubclass(Vectors.__class__, stat_tool._stat_tool._Vectors)):
        msg = "Bad type for Vectors: " + str(Vectors.__class__)
        raise TypeError, msg
    check_error.CheckType([Identifier], [str])
    # string identifying the stat_tool._stat_tool.DistributionIdentifierType
    ciIdent = "stat_tool._stat_tool.DistributionIdentifierType."
    ciIdent += Identifier.upper()
    try:
        iIdent =  eval(ciIdent)
    except AttributeError:
        msg = "Bad name for family of univariate distributions: "
        msg += Identifier
        raise ValueError, msg
    res = _dmdistributions.CompoundMultinomialEstimation(Vectors, iIdent, iInfBound, InfBoundFlag)
    return MultinomialCompoundDiscreteParametric(res)

def _DiscreteMultivariateParametricId(Identifier):
    """Return the identifier of a  parametric discrete multivariate distribution
       identified by a string

    :Parameters:
        Identifier` (str)  - name of the family of parametric discrete multivariate distributions
                among NegativeMultinomial, MultivariatePoisson, Multinomial, Iid and CompoundMultinomial
       """
    if (Identifier.upper() == "NEGATIVEMULTINOMIAL"):
        iIdent = PMDistributionIds.MNEGATIVE_MULTINOMIAL
    elif (Identifier.upper() == "MULTIVARIATEPOISSON"):
        iIdent = PMDistributionIds.MPOISSON
    elif (Identifier.upper() == "MULTINOMIAL"):
        iIdent = PMDistributionIds.MMULTINOMIAL
    elif (Identifier.upper() == "IID"):
        iIdent = PMDistributionIds.MIID
    elif (Identifier.upper() == "COMPOUNDMULTINOMIAL"):
        iIdent = PMDistributionIds.MCOMPOUND_MULTINOMIAL
    else:
        msg = "Unknown family of parametric discrete multivariate distributions: "
        msg += str(Identifier)
        raise stat_tool.error.FormatError(msg)
    return iIdent
