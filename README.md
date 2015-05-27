MB-pol plugin for OpenMM
=======================

[![Build Status](https://travis-ci.org/paesanilab/mbpol_openmm_plugin.svg?branch=master)](https://travis-ci.org/paesanilab/mbpol_openmm_plugin)

`mbpol` is a plugin for the `OpenMM` toolkit for molecular simulations using the `MB-pol` potential energy surface for water.

Previously, we carried out a detailed analysis of the two- and three-body water interactions evaluated at the CCSD(T) level to quantitatively assess the accuracy of several force fields, DFT models, and ab initio based interaction potentials that are commonly used in molecular simulations [1]. On the basis of this analysis, we have developed `MB-pol`, a "first principles" potential energy surface built upon a many-body expansion of the interaction energy of water molecules [2,3]. In particular, `MB-pol` explicitly treats the one-body (intramolecular distortion energy) term and the short-ranged two- and three-body terms through permutationally invariant polynomials. Long-range two-body interactions are described through electrostatic interactions and by dispersion energies. Higher-order terms in the expansion are modeled through many-body polarizable electrostatic interactions.

Without containing any empirical parameters, `MB-pol` has been found to accurately describe the dimer vibration-rotation tunneling spectrum to within a few wavenumbers of experiment [2], capture the second and third virial coefficients [2,3], and obtain good agreement with benchmark calculations of cluster structures and energies [3]. Finally, after including nuclear quantum effects into molecular simulations, we have demonstrated that `MB-pol` provides a highly accurate description of the liquid phase of water at ambient conditions in comparison with experiment for several structural, thermodynamic, and dynamical properties [4].

References:

1. [G.R. Medders, V. Babin, and F. Paesani, J. Chem. Theory Comput. 9, 1103–1114 (2013)](http://pubs.acs.org/doi/abs/10.1021/ct300913g).
2. [V. Babin, C. Leforestier, and F. Paesani, J. Chem. Theory Comput. 9, 5395–5403 (2013)](http://pubs.acs.org/doi/abs/10.1021/ct400863t).
3. [V. Babin, G.R. Medders, and F. Paesani, J. Chem. Theory Comput. 10, 1599–1607 (2014)](http://pubs.acs.org/doi/abs/10.1021/ct500079y).
4. [G.R. Medders, V. Babin, and F. Paesani, J. Chem. Theory Comput. 10, 2906–2910 (2014)](http://pubs.acs.org/doi/abs/10.1021/ct5004115).


available components are:

* `MBPolReferenceElectrostaticsForce`
* `MBPolReferenceOneBodyForce`
* `MBPolReferenceTwoBodyForce`
* `MBPolReferenceThreeBodyForce`
* `MBPolReferenceDispersionForce`

The parameters of each component are defined in [`python/mbpol.xml`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/mbpol.xml).
As of version `1.0`, only the `Reference` platform, i.e. single threaded C++ on CPU, is supported. It currently can simulate clusters of water molecules and bulk water with Particle Mesh Ewald (PME).

## How to install

See the [`INSTALL.md`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/INSTALL.md) file in the package.

## Automated script

`OpenMM` uses a framework approach, users need to create a Python script for setting up and running simulations, the easiest way is to create a script with the online configuration tool <http://builder.openmm.org>, however the builder does not support `MBPol`.
In order to simplify running `MBPol` for new `OpenMM` users we created a script that loads a configuration file, runs `OpenMM` and writes outputs to disk.

The script covers these use-cases:

* Cluster simulations: energy computation, geometry optimization and constant energy simulations
* Bulk simulations: cubic box Periodic Mesh Ewald energy computation, constant energy, constant temperature and constant pressure simulations

How to use the `run_mbpol` script:

* Create a new folder for the simulation
* Copy a `pdb` format initial positions file in the folder, this plugin provides a sample `water14_cluster.pdb` and `water256_bulk.pdb` in the `python/` folder
* Copy a sample configuration file in the folder, the plugin provides `mbpol_config.ini` and `mbpol_example_bulk.ini`
* Modify the configuration file parameters, see the comments
* run `run_mbpol yourconfiguration.ini`

## Example simulation

Simulation of a cluster of 14 water molecules:

* [IPython Notebook example on `nbviewer`](http://nbviewer.ipython.org/gist/zonca/54c7040c1cf3f583930f)
* IPython Notebook example source: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.ipynb>
* Python Script: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.py>
* C++ example:
  <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/examples/Simulate14WaterCluster/Simulate14WaterCluster.cpp>, see the documentation in the `examples/` folder.

## Support

* For installation issues or software bugs open an issue on Github: <https://github.com/paesanilab/mbpol_openmm_plugin/issues>
* For other questions email the main developer Andrea Zonca (username: zonca, domain: sdsc.edu)
