""" Cluster tests"""
__revision__ = "$Id: test_cluster.py 6258 2009-04-22 15:27:12Z cokelaer $"


from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.cluster import Transcode, Clustering, \
    ToDistanceMatrix, Cluster

from tools import interface

class Test(interface):
    
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/distribution1.dist",
                           ToDistanceMatrix)
    
    def build_data(self):
        vec10 = Vectors("data/chene_sessile.vec")
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        c1 = Clustering(matrix10, "Partition", 3, Prototypes=[1, 3, 12],
                        Algorithm="Divisive")
       # c1_bis = Clustering(matrix10, "Partition", 3, Prototypes=[1, 3, 12],
#                            Algorithm="Ordering")

#        c2 = Clustering(matrix10, "Hierarchy", Algorithm="Agglomerative")
#        c3 = Clustering(matrix10, "Hierarchy", Algorithm="Divisive")
#        c4 = Clustering(matrix10, "Hierarchy", Algorithm="Ordering")
        return ToDistanceMatrix(c1)

    def test_len(self):
        pass

    def test_empty(self):
        pass
        #self.empty()

    def test_constructor_from_file(self):
        pass
        #self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        pass
        #self.constructor_from_file_failure()
        
    def test_print(self):
        self.print_data()
    
    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()
    
    def test_plot(self):        
        self.plot()

    def test_save(self):
        pass
        #self.save()
                
    def test_extract(self):
        pass
    
    def test_extract_data(self):
        pass
    
    def test_symmetrize(self):
        data = self.data
        assert data.symmetrize()

    def test_unnormalize(self):
        data = self.data
        assert data.symmetrize().unnormalize()

    def test_select_individual(self):
        data = self.data
        keep_false = data.select_individual([1], keep=False)
        keep_true = data.select_individual([1], keep=True)

        assert keep_false.get_nb_row() == 2
        assert keep_false.get_nb_column() == 2
        
        assert keep_true.get_nb_row() == 1
        assert keep_true.get_nb_column() == 1
        
    def test_test_symmetry(self):
        data = self.data
        assert data.test_symmetry()
        
    def test_get_distance(self):
        data = self.data
        data.get_distance(0,0)

    def test_get_length(self):
        data = self.data
        data.get_length(0,0)

