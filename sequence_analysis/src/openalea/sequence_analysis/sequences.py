"""Sequences"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
import _sequence_analysis


from _sequence_analysis import _Sequences
#_Sequences = csequence._Sequences

__all__ = ['Sequences',
           '_Sequences']


# Extend dynamically class
interface.extend_class( _Sequences, interface.StatInterface)

# Add methods to _Vectors


def Sequences(*args, **kargs):
    """Construction of a set of sequences from multidimensional arrays of integers, from data generated by a renewal process or from an ASCII file.

    The data structure of type array(array(array(int))) should be constituted at the most internal level of arrays of constant size. If the optional argument IndexParameter is set at "Position" or "Time", the data structure of type array(array(array(int))) is constituted at the most internal level of arrays of size 1+n (index parameter, n variables attached to the explicit index parameter). If the optional argument IndexParameter is set at "Position", only the index parameter of the last array of size 1+n is considered and the first component of successive elementary arrays (representing the index parameter) should be increasing. If the optional argument IndexParameter is set at "Time", the first component of successive elementary arrays should be strictly increasing.
  
    :Parameters:
        array1 (array(array(int))): input data for univariate sequences
        arrayn (array(array(array(int)))): input data for multivariate sequences,
        timev (renewal_data),
        file_name (string).
    
    :Optional Parameters:
    Identifiers (array(int)): explicit identifiers of sequences. This optional argument can only be used if the first argument is of type array(array(int/array(int))).
    IndexParameter (string): type of the explicit index parameter: "Position" or "Time" (the default: implicit discrete index parameter starting at 0). This optional argument can only be used if the first argument is of type array(array(int/array(int))).
    
    :Returns:
        If the construction succeeds, an object of type sequences or discrete_sequences is returned, otherwise no object is returned. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is <= 15.

    :Examples:
        >>> Sequences(array1, Identifiers->[1, 8, 12])
        >>> Sequences(arrayn, Identifiers->[1, 8, 12],
        >>> IndexParameter->"Position")
        >>> Sequences(timev)
        >>> Sequences(file_name)    
    
    .. seealso::    
        :class:`~openalea.stat_tool.output.Save`,
        `ExtractHistogram`, 
        `ExtractVectors`, 
        `AddAbsorbingRun`,
        `Cluster`, 
        `Cumulate`, 
        `Difference`, 
        `Indexextract`, 
        `Lengthselect`, 
        `Merge`, 
        `MergeVariable`, 
        `MovingAverage`, 
        `RecurrenceTimeSequences`, 
        `RemoveRun`, 
        `Reverse`, 
        `SegmentationExtract`, 
        `SelectIndividual`, 
        `SelectVariable`, 
        `Shift`, 
        `Transcode`, 
        `ValueSelect`,
        `VariableScaling`, 
        `ComputeCorrelation`, 
        `ComputePartialAutoCorrelation`, 
        `ComputeSelfTransition`, 
        `Compare` (sequences), 
        `Compare` (Markovian models of seuqences), 
        `Compare` (Markovian models), 
        `Estimate` (Markovian models), 
        `ComputeStateSequences`, 
        `Simulate` (Markovian models).
    """ 
    
    type_map = {
        "INT": _sequence_analysis.INT_VALUE, 
        "REAL" : _sequence_analysis.REAL_VALUE,
        "STATE": _sequence_analysis.STATE,
        "NB_INTERNODE":  _sequence_analysis.NB_INTERNODE,
        "AUXILIARY":  _sequence_analysis.AUXILIARY,
        }
    
    index_parameter_type_map = {
        "IMPLICIT_TYPE": _sequence_analysis.IMPLICIT_TYPE,
        "TIME": _sequence_analysis.TIME,
        "TIME_INTERVAL": _sequence_analysis.TIME_INTERVAL,
        "POSITION": _sequence_analysis.POSITION,
        "POSITION_INTERVAL": _sequence_analysis.POSITION_INTERVAL
        }
 
    
    # by default, we use implicit
    index_parameter_type = kargs.get("IndexParameterType", "IMPLICIT_TYPE")
    #wrapper is looking for int or real 
    #type = kargs.get("Type","REAL_INT")    
        
    try:
        index_parameter_type = index_parameter_type_map[index_parameter_type]
    except KeyError:
        raise KeyError("Possible types are : " + 
                       str(index_parameter_type_map.keys()))
    
    
    # First case: a filename constructor. So, let us check we have only one 
    # argument, which is a string
    if len(args)==1 and isinstance(args[0], str):
        filename = args[0]
        if os.path.isfile(filename):
            return _Sequences(filename)
        else:
            raise IOError("bad file name")
    # otherwise, we switch to a list constructor that requires a list of sequences
    # and a list of identifiers. The latter being optional
    elif len(args)==1 and isinstance(args[0], list):
        return _Sequences(args[0], range(0,len(args[0])), index_parameter_type)
    # or may be provided by the user.
    elif len(args)==2 and isinstance(args[0], list) and isinstance(args[1], list):
        #if len(args[0])!=len(args[1]):
        #    raise TypeError("Expect the list of sequences and list of identifiers to have the same length")
        return _Sequences(args[0], args[1], index_parameter_type)
    else:
        raise TypeError("Expected a valid filename or a list of lists (e.g., [[1,0],[0,1]])")
    
    




