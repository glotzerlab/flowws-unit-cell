[![ReadTheDocs](https://img.shields.io/readthedocs/flowws-unit-cell.svg?style=flat)](https://flowws-unit-cell.readthedocs.io/en/latest/)

## Introduction

`flowws-unit-cell` is a set of modules to identify crystalline unit
cells. At a high level, this analysis proceeds in 4 steps:

1. Manually select a clean, single grain of well-ordered particles
2. Select three vectors specifying the periodic directions of the crystal
3. Project the observations into the unit cell and cluster the resulting coordinates
4. Detect the space group and center the system accordingly

`flowws-unit-cell` implements this workflow interactively in the
desktop or jupyter notebook as a set of modules using
[flowws-analysis](https://flowws-analysis.readthedocs.io).

## Installation

Install `flowws-unit-cell` from source:

```
pip install git+https://github.com/glotzerlab/flowws-unit-cell.git#egg=flowws-unit-cell
```

## API Documentation

Browse more detailed documentation
[online](https://flowws-unit-cell.readthedocs.io) or build the sphinx
documentation from source:

```
git clone https://github.com/glotzerlab/flowws-unit-cell
cd flowws-unit-cell/doc
pip install -r requirements.txt
make html
```
