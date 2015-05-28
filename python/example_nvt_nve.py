from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import mbpol

pdb = app.PDBFile('water256_bulk.pdb')
boxDim = 19.3996888399961804/10.
temperature = 300
boxSize = (boxDim, boxDim, boxDim) * unit.nanometer
pdb.topology.setUnitCellDimensions(boxSize)
forcefield = app.ForceField(mbpol.__file__.replace('mbpol.py', 'mbpol.xml'))

ewaldErrorTolerance = 1e-8
timestep = 2*unit.femtoseconds
production_steps = 5000

################################## NVT

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=0.9*unit.nanometers, constraints=None, rigidWater=True,
    ewaldErrorTolerance=ewaldErrorTolerance)
integrator = mm.VerletIntegrator(timestep)

system.addForce(mm.AndersenThermostat(temperature, 1./unit.picoseconds))

platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.computeVirtualSites()

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(10)

simulation.reporters.append(app.StateDataReporter("mbpol_nvt.log", 50, step=True, time=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
    progress=True, remainingTime=True, speed=True, totalSteps=production_steps,
    separator='\t'))

print('Running Production...')
simulation.step(production_steps)
print('Done!')

## Save state to be used as starting configuration for NVE

final_nvt_state = simulation.context.getState(getVelocities=True, getPositions=True)
positions = final_nvt_state.getPositions()
velocities = final_nvt_state.getVelocities()

################################## NVE

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=0.9*unit.nanometers, constraints=None, rigidWater=True,
    ewaldErrorTolerance=ewaldErrorTolerance)
integrator = mm.VerletIntegrator(timestep)

simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(positions)
simulation.context.computeVirtualSites()
simulation.context.setVelocities(velocities)

simulation.reporters.append(app.StateDataReporter("mbpol_nve.log", 50, step=True, time=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
    progress=True, remainingTime=True, speed=True, totalSteps=production_steps,
    separator='\t'))

print('Running Production...')
simulation.step(production_steps)
print('Done!')
