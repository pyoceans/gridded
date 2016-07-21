from setuptools import setup, find_packages
from gridded import __version__

reqs = [line.strip() for line in open('requirements.txt')]
def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name             = "gridded",
    version          = __version__,
    description      = "API for interpolation on regular grid, curvilinear orthogonal grid, unstructured grid",
    long_description = readme(),
    license          = "MIT License",
    author           = "Rob Hetland, Rich Signell, Kyle Wilcox",
    author_email     = "hetland@tamu.edu, rsignell@usgs.gov, kyle@axiomdatascience.com",
    url              = "https://github.com/pyoceans/gridded",
    packages         = find_packages(),
    install_requires = reqs,
    tests_require    = ['pytest'],
    classifiers      = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
)

