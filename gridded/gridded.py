class Gridded:

    _grid_obj_classes = []
    _grids_loaded = False

    @classmethod
    def _load_grid_objs(cls):
        from pkg_resources import working_set
        for ep in working_set.iter_entry_points('gridded.grid_objects'):
            cls._grid_obj_classes.append(ep.load())

    @classmethod
    def load(cls, nc, *args, **kwargs):
        for go in self._grid_obj_classes:
            if hasattr(go, 'is_mine') and go.is_mine(nc):
                return go(nc, *args, **kwargs)

