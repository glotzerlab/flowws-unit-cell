
.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. mdinclude:: ../../README.md

Usage
=====

The modules provided in this project are designed to be used as part
of a live :std:doc:`flowws-analysis<flowws-analysis:index>` workflow,
sandwiched between data sources and interactive backends and/or data
sinks. Consult the :std:doc:`flowws<flowws:index>` documentation or
the `flowws-examples
project<https://github.com/klarh/flowws-examples>`_ for information on
how to create and execute workflows.

A suggested workflow:

- Read in data (:py:class:`flowws_analysis.Pyriodic` for simple test data, otherwise read in a file with :py:class:`flowws_analysis.GTAR` or :py:class:`flowws_analysis.Garnett`)
- Select a clean crystalline grain (:py:class:`flowws_analysis.Selection` and :py:class:`flowws_analysis.Plato`)
- Compute the RDF (:py:class:`flowws_freud.RDF`)
- Select directions and distances to form the basis of the unit cell (:py:class:`flowws_unit_cell.BasisSelection`)
- Project all particles into the smaller unit cell and cluster (:py:class:`flowws_unit_cell.Projection`)
- Center the system according to the automatically-detected space group (:py:class:`flowws_unit_cell.CenterSpaceGroup`)
- Visualize the resulting unit cell (:py:class:`flowws_analysis.Plato`)
- Display plots and add interactivity (:py:class:`flowws_analysis.ViewNotebook` or :py:class:`flowws_analysis.ViewQt`)
- Optionally, save the final unit cell (:py:class:`flowws_analysis.SaveGarnett`)

Modules
=======

.. automodule:: flowws_unit_cell
   :members: BasisSelection, CenterSpaceGroup, Projection

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
