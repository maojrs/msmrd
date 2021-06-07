import numpy as np
import msmrd2
import msmrd2.visualization as msmrdvis
from msmrd2.integrators import langevin
from msmrd2.potentials import WCA
import msmrd2.tools.particleTools as particleTools
import msmrd2.tools.analysis as analysisTools

import multiprocessing
from multiprocessing import Pool
from functools import partial
import random
import sys
import os

# Main parameters
numBathParticles = 500
numparticles = 1 + numBathParticles #Added distinguished particle (index 0)
D = 0.1
particlemass = 1.0
separationDistance = 2 # minimum separation distance for initial condition
numSimulations = 250 #500

# Main parameters for integrator
dt = 0.005 #0.005
seed = -1 # Seed = -1 used random device as seed
bodytype = 'point'

# Define simulation boundaries (choose either spherical or box)
boxsize = 25
boundaryType = 'periodic'

# Parameters for WCA potential (rm=2^(1/6)sigma)
epsilon = 1.0
rm = 1.0
sigma = rm * 2**(-1/6)

# Simulation parameters
timesteps = 20000 #10000000 #3000000 #3000000
bufferSize = 1024
stride = 25
outTxt = False
outH5 = True
outChunked = True
trajtype = "positionType"


# Parent directory location
parentDirectory = "../../data/MoriZwanzig/"

# Create folder for data
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("Folder MoriZwanzig already exists.")
    proceed = True

# Create folder for benchmark data        
foldername = "benchmark"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder MoriZwanzig/benchmark already exists. Previous data files might be overwritten. Continue, y/n?")
    proceed = input()
    if proceed != 'y':
        sys.exit()

# Create parameter dictionary to write to parameters reference file
parameterfilename = os.path.join(filedirectory, "parameters")
parameterDictionary = {'numFiles' : numSimulations, 'numParticles' : numparticles, 'dt' : dt, 'bodytype' : bodytype,
                       'D' : D, 'sigma' : sigma, 'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype,
                       'boxsize' : boxsize, 'boundaryType' : boundaryType}
analysisTools.writeParameters(parameterfilename, parameterDictionary)


# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simMoriZwanzig")

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list
    seed = int(simnumber)
    random.seed(seed)
    partlist = particleTools.randomLangevinParticleList(numparticles, boxsize, separationDistance, D,
                                                        particlemass, seed)
    # Set distinguished particle (default type is zero)
    partlist[0].setType(1)

    # Define boundary
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)

    # Define potential
    potentialWCA = WCA(epsilon, sigma)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
    integrator = langevin(dt, seed, bodytype)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialWCA)

    # Creates simulation
    sim = msmrd2.simulation(integrator)

    # Output filename definitionconj
    filename = basefilename + "_{:04d}".format(simnumber)

    # Runs simulation
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked, trajtype)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
#num_cores = multiprocessing.cpu_count()
#pool = Pool(processes=num_cores)
#iterator = [i for i in range(numSimulations)]
#pool.map(partial(runParallelSims), iterator)

# Run serial for debugging
for i in range(numSimulations):
    runParallelSims(i)
