FEFF85 EXAFS C++ Documentation
===============================

.. warning::

   This C++ port is under active development and testing. It has been validated
   for **K-edge** and **L3-edge EXAFS** calculations, including **polarized**
   calculations. Testing against the original Fortran baseline is ongoing.
   **Use at your own risk** and always verify results against known references
   for your system.

FEFF85 EXAFS C++ is a complete C++ port of the FEFF85 code for *ab initio*
calculations of X-ray Absorption Fine Structure (EXAFS). It reproduces the
results of the original Fortran FEFF85 code and includes integration tests
against the Fortran baseline for 8 materials in both SCF and non-SCF modes.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   usage
   testing

.. toctree::
   :maxdepth: 2
   :caption: Theory & Architecture

   theory
   constants

.. toctree::
   :maxdepth: 2
   :caption: Module Reference

   modules
