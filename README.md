MB-pol plugin for OpenMM
=======================

`mbpol` is a plugin for the `OpenMM` toolkit for molecular simulations.

This plugin implements the `MB-pol` potential energy surface for water. MB-pol is built upon a many-body expansion of the interaction energy of water molecules and is fit exclusively to ``first principles'' calculations [1,2]. In particular, MB-pol explicitly treats the one-body (intramolecular distortion energy) term and the short-ranged two- and three-body terms through permutationally invariant polynomials. Long-range two-body interactions are described through electrostatic interactions and by dispersion energies. Finally, higher-order terms in the expansion are modeled through many-body polarizable electrostatic interactions.

Without containing any empirical parameters, MB-pol has been found to accurately describe the dimer vibration-rotation tunneling spectrum to within a few wavenumbers of experiment [1], capture the second and third virial coefficients [1,2], and obtain good agreement with benchmark calculations of cluster structures and energies [2]. Comparisons with the available experimental data for several structural, thermodynamic, and dynamical properties indicate that MB-pol provides a highly accurate description of the liquid phase of water at ambient conditions [3].

References:

1. V. Babin, C. Leforestier, and F. Paesani, J. Chem. Theory Comput. 9, 5395–5403 (2013).
2. V. Babin, G.R. Medders, and F. Paesani, J. Chem. Theory Comput. 10, 1599–1607 (2014).
3. G.R. Medders, V. Babin, and F. Paesani, J. Chem. Theory Comput. 10, 2906–2910 (2014).


available components are:

* `MBPolReferenceElectrostaticsForce`
* `MBPolReferenceOneBodyForce`
* `MBPolReferenceTwoBodyForce`
* `MBPolReferenceThreeBodyForce`
* `MBPolReferenceDispersionForce`

The parameters of each component are defined in [`python/mbpol.xml`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/mbpol.xml).
As of version `0.4.0`, only the `Reference` platform, i.e. single threaded C++ on CPU, is supported.

## How to install

See the [`INSTALL.md`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/INSTALL.md) file in the package.

## Example simulation

Simulation of a cluster of 14 water molecules:

* [IPython Notebook example on `nbviewer`](http://nbviewer.ipython.org/gist/zonca/54c7040c1cf3f583930f)
* IPython Notebook example source: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.ipynb>
* Python Script: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.py>
* C++ example:
  <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/examples/Simulate14WaterCluster/Simulate14WaterCluster.cpp>, see the documentation in the `examples/` folder.
