"""comparison test

histo done
vector done

.. todo:: sequence and markov
"""
__revision__ = "$Id$"

from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.comparison import Compare, ComparisonTest

import os

class TestHisto:
    """a simple unittest class"""

    def __init__(self):
        self.meri1 = Histogram("meri1.his")
        self.meri2 = Histogram("meri2.his")
        self.meri3 = Histogram("meri3.his")
        
    def test_comparisontest(self):
        meri1 = self.meri1
        meri2 = self.meri2

        assert ComparisonTest("F", meri1, meri2) == meri1.f_comparison(meri2)
        assert ComparisonTest("T", meri1, meri2) == meri1.t_comparison(meri2)
        assert ComparisonTest("W", meri1, meri2) == meri1.wmw_comparison(meri2)

    def test_comparison_histo(self):
        # check both the long and short argument (O, S, N)
        meri1 = self.meri1
        meri2 = self.meri2
        meri3 = self.meri3

        c1 = Compare(meri1, meri2, meri3, 'N')
        c2 = Compare(meri1, meri2, meri3, 'O')
        c3 = Compare(meri1, meri2, meri3, 'S')

        c1_long = Compare(meri1, meri2, meri3, 'NUMERIC')
        c2_long = Compare(meri1, meri2, meri3, 'ORDINAL')
        c3_long = Compare(meri1, meri2, meri3, 'SYMBOLIC')

        assert c1 == c1_long
        assert c2 == c2_long
        assert c3 == c3_long
        
        assert meri1.compare_histo(meri2, meri3, "N") == c1
        assert meri1.compare_histo(meri2, meri3, "S") == c3 
        assert meri1.compare_histo(meri2, meri3, "O") == c2
        
        assert meri1.compare_histo(meri2, meri3, "N") == c1_long
        assert meri1.compare_histo(meri2, meri3, "S") == c3_long
        assert meri1.compare_histo(meri2, meri3, "O") == c2_long
        
    def test_comparison_histo_filename(self):
        meri1 = self.meri1
        meri2 = self.meri2
        meri3 = self.meri3
        
        c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat')
        os.remove('result.dat')
        c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat', \
                     Format='ASCII')
        os.remove('result.dat')
        c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat', \
                     Format='SpreadSheet')
        os.remove('result.dat')
        try:
            c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat', \
                         Format='badname')
            assert False
        except TypeError:
            assert True
    
class TestVectors:
    
    def __init__(self):
        pass
        
    def test_compare_vectors(self):
        
        vec10 = Vectors("chene_sessile.vec")
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        assert vec15

        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        assert matrix10
        assert str(vec15.compare(VectorDistance("N", "N", "N"))) == str(matrix10)


class TestSequence:
    def test_compare_sequence(self):
        pass
        
class TestMarkov:
    def test_compare_markov(self):
        pass



