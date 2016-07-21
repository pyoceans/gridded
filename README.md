[![Build Status](https://travis-ci.org/pyoceans/gridded.svg?branch=master)](https://travis-ci.org/pyoceans/gridded)

# gridded
A high-level general-use API for working with regular, curvilinear, and unstructured grids.

This project is under **heavy development** and is likely to change drastically until the 1.0 release.

## Development (conda)

### Setup

```
conda create -n gridded -c conda-forge --file requirements.txt
source activate gridded
```

NOTE: If running Python less than 3.5, you should install `cyordereddict`:

```
conda install -c conda-forge cyordereddict
```

#### Tests

Be sure to run `python setup.py develop`. Then you can run the suite with `py.test`

### Creating Your Own gridded-Compatible object


`gridded` uses Python's setuptools `entry_points` mechanism to discover gridded-Compatible GridAdapter objects. To provide a GridAdapter from your package, create an object similar to:

```python
class SomeGridAdapter(object):
    ...
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """
        Check passed arguments to see if this class should interpret the arguments
        as being compatible with itself. For example, a `CF` adapter object
        may check for the global attribute `Conventions` in this method if the
        first argument is a netCDF4-python Dataset object or an xarray Dataset
        object. Return True if so. The possible argument are yet to be defined
        but will most likely be one of the following:
            * xarray Dataset
            * netCDF4-python Dataset
        """
        if (...):
            return True
    ...
```

and in your package's `setup.py`, include the following in your `setup()` object:

```python
setup(
    ...
    entry_points = {
        'gridded.grid_adapters': [
            'someGridAdapter = python.package.module:SomeGridAdapter'
        ]
    }
    ...
```

The API which a GridAdapter must conform to to be compatible with `gridded` is
still in development. Check back soon!
