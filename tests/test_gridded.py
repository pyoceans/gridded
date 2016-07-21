from unittest import TestCase

from gridded.gridded import Gridded


class TestGridAdapter(object):
    @classmethod
    def is_mine(cls, prime, *args, **kwargs):
        if isinstance(prime, bool) and prime:
            return True

        return False

    def __init__(self, *args, **kwargs):
        pass


class TestGridded(TestCase):
    def setUp(self):
        Gridded._grid_obj_classes = [ TestGridAdapter ]

    def test_grid_adapter_loaded(self):
        assert TestGridAdapter in Gridded._grid_obj_classes

    def test_dimension_adapter_load(self):
        lda = Gridded.load(True)
        assert isinstance(lda, TestGridAdapter)

    def test_gridded_no_load(self):
        lda = Gridded.load(False)
        assert lda is None
