from unittest import TestCase

from gridded.gridded import Gridded, TestDimensionAdapter

class TestGridded(TestCase):
    def setUp(self):
        Gridded._load_grid_objs()

    def test_test_dimension_adapter_loaded(self):
        assert TestDimensionAdapter in Gridded._grid_obj_classes

    def test_test_dimension_adapter_load(self):
        lda = Gridded.load(True)
        assert isinstance(lda, TestDimensionAdapter)

    def test_gridded_no_load(self):
        lda = Gridded.load(False)
        assert lda == None

