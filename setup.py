#!/usr/bin/env python

import os
from setuptools import setup

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

version_fname = os.path.join(THIS_DIR, 'flowws_unit_cell', 'version.py')
with open(version_fname) as version_file:
    exec(version_file.read())

readme_fname = os.path.join(THIS_DIR, 'README.md')
with open(readme_fname) as readme_file:
    long_description = readme_file.read()

module_names = [
    'BasisSelection',
    'CenterSpaceGroup',
    'Projection',
]

flowws_modules = []
for name in module_names:
    flowws_modules.append('{0} = flowws_unit_cell.{0}:{0}'.format(name))
    flowws_modules.append(
        'flowws_unit_cell.{0} = flowws_unit_cell.{0}:{0}'.format(name))

setup(name='flowws-unit-cell',
      author='Matthew Spellings',
      author_email='mspells@umich.edu',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 3',
      ],
      description='Stage-based scientific workflows for crystal unit cell identification',
      entry_points={
          'flowws_modules': flowws_modules,
      },
      extras_require={},
      install_requires=[
          'flowws',
          'flowws-analysis',
          'flowws-freud',
          'freud-analysis',
          'plato-draw',
          'rowan',
          'scikit-learn',
          'spglib',
      ],
      license='BSD-3-Clause',
      long_description=long_description,
      long_description_content_type='text/markdown',
      packages=[
          'flowws_unit_cell',
      ],
      python_requires='>=3',
      version=__version__
      )
