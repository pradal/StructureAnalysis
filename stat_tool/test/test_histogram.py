"""histogram tests"""
__revision__ = "$Id: $"

from openalea.stat_tool.distribution import ToHistogram
from openalea.stat_tool.distribution import ToDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Distribution


class Test:
    """a simple unittest class"""

    def test_empty(self):
        """Test that empty constructor fails"""
        try:
            _h = Histogram()
            assert False
        except Exception:
            assert True
    
    def test_constructor_fromfile(self):
        """run constructor with filename argument"""
        h = Histogram("peup1.his")
        assert h

        return h
    
    def test_constructor_fromfile_failure(self):
        """run constructor with wrong filename argument"""
        try:
            _h = Histogram("peup1.hi")
            assert False
        except Exception:
            assert True
        
    def test_constructor_fromlist(self):
        """constructor from list """
        h = Histogram([0, 1, 2, 3])
        assert h
        
    def test_fromfile2(self):
        """From file"""
        import os
        h = Histogram("meri1.his")
        assert h
        # File
        h.save("test.his")

        h1 = Histogram("meri1.his")
        h2 = Histogram("test.his")
        assert len(h) == len(h2)
        assert list(h) == list(h2)
        os.remove("test.his")

    def test_ascii_and_display(self):
        """ test display"""
        h = self.test_constructor_fromfile()
        
        s = str(h)
        assert h.display() == s
        
        h.ascii_write(True)
        h.survival_ascii_write()
        
    def test_fromlist(self):
        """test constructor with list argument"""
        h = Histogram([1,2,3])
        assert h

        try:
            h = Histogram([0, 'a', 'b', 3])
            assert False
        except TypeError:
            assert True

    def test_plot(self):
        """ plot"""
        h = self.test_constructor_fromfile()
        h.plot()

    def test_len(self):
        h = Histogram(range(10))
        assert len(h) == 10
        
    def test_container(self):
        """ container / iterator"""
        h = Histogram("meri1.his")

        assert h[0] == 0
        assert h[10] == 1

    def test_to_histogram(self):
        """Test the ToHistogram function"""
        
        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        h = d.simulate(1000)
        d2 = ToDistribution(h)
        assert h and d2

        h2 = ToHistogram(d2)
        assert h2
        assert h == h2
        
    

if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    test = Test()
    for method in dir(test):
        if method.startswith('_'):
            continue
        if callable(getattr(test, method)):
            getattr(test, method)()
        else:
            print 'skipping'
    # and functional tests.    


