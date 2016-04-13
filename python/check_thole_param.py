from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

pdb = app.PDBFile("water14_cluster.pdb")
forcefield = app.ForceField("mbpol.xml")
#<MBPolElectrostaticsForce thole-charge-charge="0.4" thole-charge-dipole="0.4" thole-dipole-dipole="0.055" thole-dipole-dipole-singlebond="0.626">
expected = [0.4, 0.4, 0.055, 0.626, 0.055]
#<MBPolElectrostaticsForce thole-charge-charge="1" thole-charge-dipole="2" thole-dipole-dipole="3" thole-dipole-dipole-singlebond="4">
#expected = [1, 2, 3, 4, 3]
for param, expected_param in zip(forcefield._forces[0].thole, expected):
    assert ( param == expected_param)
print ("Test Passed")
