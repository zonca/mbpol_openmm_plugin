from __future__ import print_function
import configparser

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

if len(sys.argv) < 2:
    print("run_mbpol configurationfile.ini")
    sys.exit(1)

## Read configuration
config_filename = sys.argv[1]
config = configparser.ConfigParser()
config.read(config_filename)

assert "mbpol" in config.sections(), "The configuration file needs a [mbpol] section"

pdb = app.PDBFile(config["mbpol"]["pdb_filename"])

## Load the mbpol force field

forcefield = app.ForceField("mbpol.xml")

## Set nonbonded interaction

nonbonded = getattr(app, config["mbpol"]["nonbonded"])
if config["mbpol"]["nonbonded"] == "PME":
    boxDim = float(config["mbpol"]["pme_box_size_nm"])
    boxSize = (boxDim, boxDim, boxDim) * unit.nanometer
    pdb.topology.setUnitCellDimensions(boxSize)

## Create the System, define an integrator, define the Simulation

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, nonBondedCutoff=0.9*unit.angstrom)
integrator = mm.VerletIntegrator(float(config["mbpol"]["timestep_fs"])*unit.femtoseconds)

platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.computeVirtualSites()

## Local energy minimization

if config.getboolean("mbpol", "local_minimization"):
    from simtk.openmm import LocalEnergyMinimizer
    LocalEnergyMinimizer.minimize(simulation.context, 1e-1)

## Run a constant energy simulation (Verlet integrator)

simulation.context.setVelocitiesToTemperature(float(config["mbpol"]["temperature_k"])*unit.kelvin)
# Equilibrate
simulation.step(10)

# Add a `reporter` that prints out the simulation status every 10 steps

simulation.reporters.append(app.StateDataReporter(sys.stdout, int(config["mbpol"]["print_stdout_every"]), step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=110, separator='\t'))

# Add a `PDBReporter` that writes molecules positions every 20 steps in a pdb file.

simulation.reporters.append(app.PDBReporter('trajectory.pdb', int(config["mbpol"]["print_pdb_every"])))

# Run 100 steps

simulation.step(int(config["mbpol"]["steps"]))
