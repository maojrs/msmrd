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

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):
    # Particle list definition
    D = 1
    Drot = 1
    position1 =  np.random.uniform(-1,1,3)
    position2 =  np.random.uniform(-1,1,3)
    position3 =  np.random.uniform(-1,1,3)
    orientation1 = np.random.uniform(-1,1,4)
    orientation2 = np.random.uniform(-1,1,4)
    orientation3 = np.random.uniform(-1,1,4)
    orientation1 = orientation1/np.linalg.norm(orientation1)
    orientation2 = orientation2/np.linalg.norm(orientation2)
    orientation3 = orientation3/np.linalg.norm(orientation3)
    part1 = msmrd2.particle(D, Drot, position1, orientation1)
    part2 = msmrd2.particle(D, Drot, position2, orientation2)
    part3 = msmrd2.particle(D, Drot, position3, orientation3)
    partlist = msmrd2.integrators.particleList([part1, part2, part3])
    Nparticles = np.size(partlist)

    # Integrator definition
    dt = 0.00001 #0.001
    seed = -simnumber # random seed (negative), it is also different for every simulation, good for parallel simulation
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
    timesteps = 500000 #5000
    bufferSize = 1024
    stride = 2500 #1
    outTxt = False
    outH5 = True
    outChunked = True
    filename = "../data/trimer/simTrimer" + "{:04d}".format(simnumber)

    # Runs simulation
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
Nsimulations = 10
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
iterator = [i for i in range(Nsimulations)]
pool.map(partial(runParallelSims), iterator)