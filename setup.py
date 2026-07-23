#!/usr/bin/env python

from setuptools import setup
from glob import glob

__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Robyn Wright"

long_description = ("Genome coverage checker")

setup(name='Genome Coverage Checker',
      version=__version__,
      description=('Genome Coverage Checker'),
      maintainer=__maintainer__,
      url='https://github.com/R-Wright-1/genome_coverage_checker/wiki',
      packages=['genome_coverage_checker'],
      scripts=glob('scripts/*py'),
      install_requires=[],
      package_data={'genome_coverage_checker':
                    ['numpy',
			'h5py',
                        'joblib']},
      long_description=long_description)
