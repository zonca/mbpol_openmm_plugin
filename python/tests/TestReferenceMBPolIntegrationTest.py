from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestReferenceMBPolIntegration(unittest.TestCase):
    
    def testMBPolIntegrationTest(self):
        expected_energy = -2270.88890
        tolerance = 1.0e-2
        pdb = app.PDBFile("../water256_integration_test.pdb")
        forcefield = app.ForceField("../mbpol.xml")
        nonbondedCutoff = .09
        nonbondedMethod = app.CutoffPeriodic
        #if (nonbondedMethod == app.CutoffPeriodic):
        boxsize = [19.3996888399961804/10., 19.3996888399961804/10., 19.3996888399961804/10.]
        pdb.topology.setUnitCellDimensions( boxsize )
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod,
                                         nonBondedCutoff=nonbondedCutoff)
            
        integrator = mm.LangevinIntegrator(0.02*unit.femtoseconds)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        simulation.context.computeVirtualSites()
        state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
        potential_energy = state.getPotentialEnergy()
        potential_energy.in_units_of(unit.kilocalorie_per_mole)
        print(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value)
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value -expected_energy) < tolerance)
        
        
if __name__ == '__main__':
    unittest.main()