import msmrd2
import numpy as np
from msmrd2.integrators import overdampedLangevin as odLangevin
from msmrd2.potentials import patchyParticleAngular
import msmrd2.tools.particleTools as particleTools
import multiprocessing
from multiprocessing import Pool
from functools import partial

# Main parameters for particle and integrator
Nparticles = 2
dt = 0.00001
bodytype = "rigidbody"
D = 1
Drot = 1
separationDistance = 1.5 # minimum separation distance for initial condition
numSimulations = 10 #500 #1000 #20 #200


# Simulation parameters
timesteps = 2000000 #5000000 #400000 #20000000
bufferSize = 1024
stride = 25
outTxt = False
outH5 = True
outH5 = True
outChunked = True
trajtype = "patchyDimer"


# Define patchy Particle potential with angular dependence
sigma = 1.0
strength = 160.0 #50.0
angularStrength = 20.0 #10.0
angleDiff = 3*np.pi/5.0
patch1 = np.array([np.cos(angleDiff/2),np.sin(angleDiff/2),0.])
patch2 = np.array([np.cos(-angleDiff/2),np.sin(-angleDiff/2),0.])
patchesCoordinates = [patch1, patch2]
potentialPatchyParticleAngular = patchyParticleAngular(sigma, strength, angularStrength, patchesCoordinates)

# Define simulaion boundaries (choose either spherical or box)
boxsize = 4
boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, 'periodic')

# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = "../data/dimer/simDimer_t" + "{:.2E}".format(timesteps) + \
           "_s{:d}".format(stride)
#basefilename = "/group/ag_cmb/scratch/maojrs/msmrd2_data/dimer/simDimer_t" + \
#           "{:.2E}".format(timesteps)  + "_s{:d}".format(stride)

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list
    seed = int(simnumber)
    partlist = particleTools.randomParticleList(Nparticles, boxsize, separationDistance, D, Drot, seed)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative), it is also different for every simulation, good for parallel simulation
    integrator = odLangevin(dt, seed, bodytype)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyParticleAngular)

    # Creates simulation
    sim = msmrd2.simulation(integrator)

    # Output filename definition
    filename = basefilename + "_{:04d}".format(simnumber)

    # Runs simulation
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked, trajtype)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)
