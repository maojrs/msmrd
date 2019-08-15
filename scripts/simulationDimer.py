import msmrd2
import numpy as np
from msmrd2.integrators import overdampedLangevin as odLangevin
from msmrd2.potentials import patchyParticleAngular
import msmrd2.tools.particleTools as particleTools
import msmrd2.tools.analysis as analysisTools
import multiprocessing
from multiprocessing import Pool
from functools import partial
import sys
import os

# Main parameters for particle and integrator
Nparticles = 2
dt = 0.00001
bodytype = "rigidbody"
D = 1
Drot = 1
separationDistance = 1.25 # minimum separation distance for initial condition
numSimulations = 500 #500 #1000 #20 #200


# Simulation parameters
timesteps = 3000000 #2000000 #5000000 #400000 #20000000
bufferSize = 1024
stride = 25
outTxt = False
outH5 = True
outChunked = True
trajtype = "patchyDimer"

# Define patchy Particle potential with angular dependence
sigma = 1.0
strength = 75.0 #50.0
angularStrength = 10.0 #10.0
angleDiff = 3*np.pi/5.0
patch1 = np.array([np.cos(angleDiff/2),np.sin(angleDiff/2),0.])
patch2 = np.array([np.cos(-angleDiff/2),np.sin(-angleDiff/2),0.])
patchesCoordinates = [patch1, patch2]
potentialPatchyParticleAngular = patchyParticleAngular(sigma, strength, angularStrength, patchesCoordinates)

# Define simulaion boundaries (choose either spherical or box)
boxsize = 5
boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, 'periodic')

# Parent directory location
parentDirectory = "../data/dimer/"
# parentDirectory = "/group/ag_cmb/scratch/maojrs/msmrd2_data/dimer/"

# Create folder for data
foldername = "75strength"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder already exists. Previous data files might be overwritten. Continue, y/n?")
    proceed = input()
    if proceed != 'y':
        sys.exit()

# Create parameter dictionary to write to parameters reference file
parameterfilename = os.path.join(filedirectory, "parameters")
parameterDictionary = {'Nparticles' : Nparticles, 'dt' : dt, 'bodytype' : bodytype, 'D' : D, 'Drot' : Drot,
                       'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype, 'boxsize' : boxsize,
                       'potentialStrength' : strength, 'potentialAngularStrength' : angularStrength,
                       'potentialPatchesAngleDiff' : angleDiff}
analysisTools.writeParameters(parameterfilename, parameterDictionary, pickleDictionary = True)


# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simDimer")

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
