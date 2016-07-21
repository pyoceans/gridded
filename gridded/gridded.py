class Gridded:

    _grid_obj_classes = []
    _grids_loaded = False

    @classmethod
    def _load_grid_objs(cls):
        from pkg_resources import working_set
        for ep in working_set.iter_entry_points('gridded.grid_objects'):
            cls._grid_obj_classes.append(ep.load())

    @classmethod
    def load(cls, *args, **kwargs):
        for go in cls._grid_obj_classes:
            if hasattr(go, 'is_mine') and go.is_mine(*args, **kwargs):
                return go(*args, **kwargs)

class TestGridObject:
    @classmethod
    def is_mine(cls, nc, *args, **kwargs):
        return True

