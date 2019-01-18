# Tests writing trajectory to .txt and .H5 file format
import msmrd2
import numpy as np
import h5py
from msmrd2.integrators import overdampedLangevin as odLangevin
from msmrd2.potentials import patchyParticleAngular
#from joblib import Parallel, delayed
import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
from functools import partial

# Particle list definition
D = 1
Drot = 1
position = np.array([0.5,0,0])
position2 = np.array([-0.5,0,0])
position3 = np.array([0.0,0,0.5])
orientation = np.array([1,0,0,0])
orientation3 = np.array([-1,0,0,0])
part1 = msmrd2.particle(D, Drot, position, orientation)
part2 = msmrd2.particle(D, Drot, position2, orientation)
part3 = msmrd2.particle(D, Drot, position3, orientation3)
partlist = msmrd2.integrators.particleList([part1, part2, part3])
Nparticles = np.size(partlist)



# Simulation wrapper for parallel runs
def runParallelSims(simnumber):
    # Integrator definition
    dt = 0.001
    seed = -simnumber # uses random seed, good for parallel simulations
    bodytype = "rigidbody"
    integrator = odLangevin(dt, seed, bodytype)

    # Define boundary (choose either spherical or box)
    radius = 3
    sphereBoundary = msmrd2.sphere(radius,'reflective')
    integrator.setBoundary(sphereBoundary)

    # Patchy Particle potential with angular dependence definition
    sigma = 1.0
    strength = 200.0
    patch1 = np.array([1.,0.,0.])
    patch2 = np.array([np.cos(3*np.pi/5.0),np.sin(3*np.pi/5.0),0.])
    patchesCoordinates = [patch1, patch2]
    potentialPatchyParticleAngular = patchyParticleAngular(sigma, strength, patchesCoordinates)
    integrator.setPairPotential(potentialPatchyParticleAngular)

    # Creates simulation
    sim = msmrd2.simulation(integrator)
    # Simulation parameters
    timesteps = 5000
    bufferSize = 1024
    stride = 1
    outTxt = False
    outH5 = True
    outChunked = True
    filename = "../data/trimer/simTrimer" + "{:04d}".format(simnumber)
    # Run simulation
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
Nsimulations = 10
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
iterator = [i for i in range(Nsimulations)]
pool.map(partial(runParallelSims), iterator)