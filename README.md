[![Build Status](https://travis-ci.org/pyoceans/gridded.svg?branch=master)](https://travis-ci.org/pyoceans/gridded)

# gridded
API for interpolation on regular grid, curvilinear orthoganal grid, unstructured grid

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

### Creating Your Own gridded-Compatible Dimension Adapter object

gridded uses Python's setuptools `entry_points` mechanism to discover gridded Dimension Adapter objects.
To provide a Dimension Adapter from your package, create an object similar to:

```python
class SomeGridGriddedDimensionAdapter(TBDBaseClass):
    ...
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """
        Check passed arguments to see if our class is compatible with the passed in objects, 
        such as a NetCDF Dataset object. Return True if so.
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
        'gridded.dimension_adapters': [
            'somegridDA = python.package.module:SomeGridGriddedDimensionAdapter'
        ]
    }
    ...
```


