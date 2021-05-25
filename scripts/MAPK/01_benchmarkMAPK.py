import numpy as np
import msmrd2
from msmrd2.potentials import patchyProteinMAPK
from msmrd2.integrators import integratorMAPK
import msmrd2.tools.particleTools as particleTools
import msmrd2.tools.analysis as analysisTools

import multiprocessing
from multiprocessing import Pool
from functools import partial
import random
import sys
import os

# Main parameters for particles and integrator
numMAPKs = 1
numKinases = 1
numPhosphatases = 0
numparticles = numMAPKs + numKinases + numPhosphatases
D = 1.0
Drot = 1.0
particleTypes = [0, 1, 2]
# Main parameters for integrator
dt = 0.00001
bodytype = 'rigidbody'
anglePatches = np.pi/2
reactivationRateK = 1.0
reactivationRateP = 1.0
minimumUnboundRadius = 1.25
numSimulations = 4 #500

# Simulation parameters
timesteps = 3000000
bufferSize = 1024
stride = 25
outTxt = False
outH5 = True
outChunked = True
trajtype = "MAPK" # "trajectoryPositionOrientationState"

# Define Patchy Protein MAPK potential parameters (This values are fixed and should match
# those used to determine metastable states in potential and trajectory.)
sigma = 1.0
strength = 100 #65
angularStrength = 10 #2
patchesCoordinates1 = [np.array([np.cos(anglePatches/2), np.sin(anglePatches/2), 0.]),
                       np.array([np.cos(-anglePatches/2), np.sin(-anglePatches/2), 0.])]
patchesCoordinates2 = [ np.array([np.cos(-anglePatches/2), np.sin(-anglePatches/2), 0.]) ]
potentialPatchyProteinMAPK = patchyProteinMAPK(sigma, strength, patchesCoordinates1, patchesCoordinates2)

# Define simulation boundaries (choose either spherical or box)
boxsize = 3
boundaryType = 'periodic'
boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)

# Parent directory location
parentDirectory = "../../data/MAPK/"

# Create folder for data
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("Folder MAPK already exists.")
    proceed = True

# Create folder for benchmark data        
foldername = "benchmark"
if numMAPKs == 1 and numKinases == 1 and numPhosphatases == 0:
    foldername = "benchmark_kinase"
elif numMAPKs == 1 and numKinases == 0 and numPhosphatases == 1:
    foldername = "benchmark_phosphatase"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder MAPK/benchmark already exists. Previous data files might be overwritten. Continue, y/n?")
    proceed = input()
    if proceed != 'y':
        sys.exit()

# Create parameter dictionary to write to parameters reference file
parameterfilename = os.path.join(filedirectory, "parameters")
parameterDictionary = {'numFiles' : numSimulations, 'numParticles' : numparticles, 'dt' : dt, 'bodytype' : bodytype,
                       'D' : D, 'Drot' : Drot, 'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype,
                       'boxsize' : boxsize, 'boundaryType' : boundaryType, 'potentialStrength' : strength,
                       'potentialAngularStrength' : angularStrength, 'anglePatches' : anglePatches}
analysisTools.writeParameters(parameterfilename, parameterDictionary)


# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simMAPK")

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list
    seed = int(simnumber)
    random.seed(seed)
    defaultOrientvector = patchesCoordinates2[0]
    partlist, mapkIndex, kinaseIndex, phosIndex = particleTools.randomMAPKparticleList(numMAPKs, numKinases,
            numPhosphatases, boxsize, minimumUnboundRadius, D, Drot, seed, defaultOrientvector)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
    integrator = integratorMAPK(dt, seed, bodytype, anglePatches, reactivationRateK, reactivationRateP,
                                mapkIndex, kinaseIndex, phosIndex)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyProteinMAPK)
    integrator.disableDeactivation()

    # Creates simulation
    sim = msmrd2.simulation(integrator)

    # Output filename definitionconj
    filename = basefilename + "_{:04d}".format(simnumber)

    # Runs simulation
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked, trajtype)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)
