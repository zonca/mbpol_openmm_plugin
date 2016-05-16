from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestCustomForce(unittest.TestCase):
    """Test the functionality of Custom Dispersion Force xml file."""

    def test_three_water(self, nonbondedMethod=app.CutoffNonPeriodic):
        expected_energy = -6.84471477
        pdb = app.PDBFile("../water3.pdb")
        forcefield = app.ForceField("../mbpol.xml")
        nonbondedCutoff = 1e3*unit.nanometer
        
        if (nonbondedMethod == app.CutoffPeriodic):
            boxsize = [50, 50, 50]
            pdb.topology.setUnitCellDimensions( boxsize )
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonBondedCutoff=nonbondedCutoff)
        
        #forces must be removed from the system to test single force
        #order of forces:
        # elec
        # one
        # two
        # three
        # CMMotionRemover
        # CustomDispersion
        system.removeForce(0) #remove elec
        system.removeForce(0) #remove one 
        system.removeForce(0) #remove two
        system.removeForce(0) #remove three
        system.removeForce(0) #remove CMMotionRemover
           
        integrator = mm.VerletIntegrator(0.02*unit.femtoseconds)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        simulation.context.computeVirtualSites()
        state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
        potential_energy = state.getPotentialEnergy()
        potential_energy.in_units_of(unit.kilocalorie_per_mole)
        
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .01)
    
    def test_water_and_ion(self):
        expected_energy = -1.306598
        pdb = app.PDBFile("pdb_files/water_and_ion.pdb")
        forcefield = app.ForceField("../mbpol.xml")
        nonbondedMethod = app.CutoffNonPeriodic
        nonbondedCutoff = 1e3*unit.nanometer
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, 
                                         nonBondedCutoff=nonbondedCutoff) 
        
        #forces must be removed from the system to test single force
        #order of forces:
        # elec
        # one
        # two
        # three
        # CMMotionRemover
        # CustomDispersion
        system.removeForce(0) #remove elec
        system.removeForce(0) #remove one 
        system.removeForce(0) #remove two
        system.removeForce(0) #remove three
        system.removeForce(0) #remove CMMotionRemover
         
        integrator = mm.VerletIntegrator(0.02*unit.femtoseconds)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        simulation.context.computeVirtualSites()
        state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
        potential_energy = state.getPotentialEnergy()
        print(potential_energy.in_units_of(unit.kilocalorie_per_mole))
        
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .01)
       
        
        
if __name__ == '__main__':
    unittest.main()
