import msmrd2
import numpy as np
from msmrd2.integrators import overdampedLangevin as odLangevin
from msmrd2.potentials import patchyParticleAngular
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
    orientation1 = np.random.uniform(-1,1,4)
    orientation2 = np.random.uniform(-1,1,4)
    orientation1 = orientation1/np.linalg.norm(orientation1)
    orientation2 = orientation2/np.linalg.norm(orientation2)
    part1 = msmrd2.particle(D, Drot, position1, orientation1)
    part2 = msmrd2.particle(D, Drot, position2, orientation2)
    partlist = msmrd2.integrators.particleList([part1, part2])
    Nparticles = np.size(partlist)

    # Integrator definition
    dt = 0.00001 #0.001
    seed = -simnumber # random seed (negative), it is also different for every simulation, good for parallel simulation
    bodytype = "rigidbody"
    integrator = odLangevin(dt, seed, bodytype)

    # Define boundary (choose either spherical or box)
    edge = 5
    boxBoundary = msmrd2.box(edge, edge, edge,'periodic')
    integrator.setBoundary(boxBoundary)

    # Patchy Particle potential with angular dependence definition
    sigma = 1.0
    strength = 160.0 # Deprecated does nothing
    angularStrength = 20.0
    angleDiff = 3*np.pi/5.0
    patch1 = np.array([np.cos(angleDiff/2),np.sin(angleDiff/2),0.])
    patch2 = np.array([np.cos(-angleDiff/2),np.sin(-angleDiff/2),0.])
    patchesCoordinates = [patch1, patch2]
    potentialPatchyParticleAngular = patchyParticleAngular(sigma, strength, angularStrength, patchesCoordinates)
    integrator.setPairPotential(potentialPatchyParticleAngular)

    # Creates simulation
    sim = msmrd2.simulation(integrator)
    # Simulation parameters
    timesteps = 50000000 #20000000 # 200000000 #5000
    bufferSize = 1024
    stride = 250 #1
    outTxt = False
    outH5 = True
    outH5 = True
    outChunked = True
    # Folder included in filename must exist (and preferably empty), otherwise H5
    # fails (might fail with parallel python) to write the data.
    filename = "../data/dimer/simDimer_t" + "{:.2E}".format(timesteps)  \
               + "_s{:d}".format(stride) + "_{:04d}".format(simnumber)
    #filename = "/group/ag_cmb/scratch/maojrs/msmrd2_data/dimer_discrete/simDimer" + "{:04d}".format(simnumber)
    trajtype = "patchyDimer"

    # Runs simulation
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked, trajtype)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
Nsimulations = 20 #200
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
iterator = [i for i in range(Nsimulations)]
pool.map(partial(runParallelSims), iterator)
