from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import os.path
import pandas as pd
import datetime
import mbpol

configurations = (("amoeba2013", "CUDA"), ("amoeba2013", "Reference"),)#, ("mbpol", "Reference"),)
configurations = (("mbpol", "Reference"),)
steps = 10000

for forcefield_configuration, platform_name in configurations:

    pdb_filename = "pdb/TIP5P_PIMC_0C_LIQ"
    pdb_filename = "pdb/H2O_1"
    if forcefield_configuration.startswith("amoeba"):
        pdb_filename += "_AMOEBA"
    pdb_filename += ".pdb"
    pdb = app.PDBFile(pdb_filename)

    tag = os.path.basename(pdb_filename).replace(".pdb", "")

    # PME box size
    if tag.startswith("H2O"):
        boxSize = 24.8343697144598003
    elif tag.startswith("FM_25C") or tag.startswith("TIP5P_25C"):
        boxSize = 19.7316565863235596
    else:
        boxSize = 19.3996888399961804
    boxSize /= 10.
    pdb.topology.setUnitCellDimensions((boxSize,boxSize,boxSize))


    forcefield = app.ForceField(forcefield_configuration + ".xml")

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=.9*unit.nanometer)
    timestep = 0.02*unit.femtoseconds
    integrator = mm.VerletIntegrator(timestep)

    platform = mm.Platform.getPlatformByName(platform_name)
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.computeVirtualSites()

    state = simulation.context.getState(getForces=True, getEnergy=True)
    potential_energy = state.getPotentialEnergy()
    energy_kcal_per_mol = potential_energy.value_in_unit(unit.kilocalorie_per_mole)

    simulation.reporters.append(app.StateDataReporter("simulation_{tag}_{timestep}.log".format(**locals()).replace(" ",""), 1, step=True, totalEnergy=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=steps, separator=','))

    output = pd.DataFrame(dict(
                    box_size = boxSize,
                    openmm_energy_kcal_per_mol = energy_kcal_per_mol
                    ), index=[tag])

    start = datetime.datetime.now()
    simulation.step(steps)
    end = datetime.datetime.now()
    print("%s: %s %s time elapsed (%d steps): %f s" % (pdb_filename, forcefield_configuration, platform_name, steps, (end-start).total_seconds()))
