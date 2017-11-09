from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import unittest

def my_test_suite():
    #For test purposes, currently there is no test script
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('test', pattern='test_*.py')
    return test_suite

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'Clusterkit',
    packages = ['Clusterkit'],
    version = '0.1',
    description = 'The code is mainly designed for applying cluster algorithms for description of energetics, measurement of thermodynamic properties and diffusion simulations.',
    long_description=long_description,
    author = 'Bin Ouyang',
    author_email = 'jeffrey.oakley.ouyang@gmail.com',
    license='MIT',  # LICENSE.txt
    url = 'https://github.com/Jeff-oakley/Clusterkit', # use the URL to the github repo
    download_url = 'https://github.com/DallasTrinkle/onsager/tarball/v0.1.tar.gz', # for when we upload
    keywords = ['Cluster algorithm', 'Solid solution', 'Disordered lattice',
                'Monte Carlo simulation', 'Atomic diffusion'],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7.11',
        'Programming Language :: Python :: 3.4.3', #Does not test rigorously
    ],
    install_requires=['numpy', 'pyyaml'],
    test_suite = 'setup.my_test_suite'
)
