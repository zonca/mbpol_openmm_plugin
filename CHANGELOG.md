# v0.4.0 (18-09-2014)

Complete code refactor based on the [OpenMM Example Plugin](https://github.com/peastman/openmmexampleplugin)  ([#21](https://github.com/paesanilab/mbpol_openmm_plugin/pull/21)), allows to:

* install the plugin on a pre-existing OpenMM install
* create Python wrappers

Other notes:

* Created IPython notebook example
* Prepared README.md and license for release

# v0.3.0 (10-06-2014)

Complete MBPol force field implementation for clusters of water molecules

* Implementation of 2B interaction [#12](https://github.com/sdsc/mbpol_openmm_plugin/pull/12)
* Implementation of ThreeBody Neighbor list that scales in linear time [#14](https://github.com/sdsc/mbpol_openmm_plugin/pull/14)
* Implementation of 3B interaction [#13](https://github.com/sdsc/mbpol_openmm_plugin/pull/13)
* Implementation of Dispersion force [#15](https://github.com/sdsc/mbpol_openmm_plugin/pull/15)
* Large code cleanup and consistent renaming [#17](https://github.com/sdsc/mbpol_openmm_plugin/pull/17)
* Implementation of a full integration test with comparison to the legacy `mbpol` code

# v0.2.0 (28-04-2014)

Complete MBPol electrostatics and 1B interaction

* Implementation of multiple Thole factors for charge-charge/charge-dipole/dipole-dipole [#7](https://github.com/sdsc/mbpol_openmm_plugin/pull/7)
* Virtual Site and charge redistribution [#3](https://github.com/sdsc/mbpol_openmm_plugin/pull/3)
* Implementation of 1B interaction based on `AmoebaReferenceStretchBendForce` [#9](https://github.com/sdsc/mbpol_openmm_plugin/pull/9)

# v0.1.0 (24-02-2014) 

* Electrostatics with no m-site, no charge redistribution and no thole factors [#2](https://github.com/sdsc/mbpol_openmm_plugin/pull/2)
