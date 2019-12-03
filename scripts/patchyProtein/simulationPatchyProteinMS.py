import msmrd2
import numpy as np
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.potentials import patchyProteinMarkovSwitch
from msmrd2.integrators import overdampedLangevinMarkovSwitch as odLangevinMS
import msmrd2.tools.particleTools as particleTools
import msmrd2.tools.analysis as analysisTools
import multiprocessing
from multiprocessing import Pool
from functools import partial
import random
import sys
import os

# Main parameters for particle and integrator
Nparticles = 2
dt = 0.00001
bodytype = "rigidbody"
particleTypes = [0, 1]
separationDistance = 1.25 # minimum separation distance for initial condition
numSimulations = 600 #500 #1000 #20 #200


# Simulation parameters
timesteps = 6000000 #2000000 #5000000 #400000 #20000000
bufferSize = 1024
stride = 50
outTxt = False
outH5 = True
outChunked = True
trajtype = "patchyProtein" #"positionOrienatationState"


# Define Patchy Protein potential parameters (This values are fixed and should match
# those used to determine metastable states in potential and trajectory.)
sigma = 1.0
strength = 65
angularStrength = 2
patchesCoordinates1 = [np.array([1.,0.,0.]), np.array([0.,1.,0.]), np.array([0.,0.,1.]),
                       np.array([-1.,0.,0.]), np.array([0.,-1.,0.]), np.array([0.,0.,-1.])]
patchesCoordinates2 = [np.array([1.,0.,0.])]


# Define simulation boundaries parameters (choose either spherical or box)
boxsize = 5
boundaryType = 'periodic'


# Parent directory location
parentDirectory = "../../data/patchyProtein/"
# parentDirectory = "/group/ag_cmb/scratch/maojrs/msmrd2_data/patchyProtein/"

# Create folder for data
foldername = "benchmark"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder already exists. Previous data files might be overwritten. Continue, y/n?")
    proceed = input()
    if proceed != 'y':
        sys.exit()

# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simPatchyProtein")

# Create parameter dictionary and writes parameters to reference file
parameterfilename = os.path.join(filedirectory, "parameters")
parameterDictionary = {'numFiles' : numSimulations, 'numParticles' : Nparticles, 'dt' : dt, 'bodytype' : bodytype,
                       'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype, 'boxsize' : boxsize,
                       'boundaryType' : boundaryType, 'potentialStrength' : strength,
                       'potentialAngularStrength' : angularStrength}
analysisTools.writeParameters(parameterfilename, parameterDictionary)


# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define seed
    seed = int(simnumber)
    random.seed(seed)

    # Create unbound MSMlist (CTMSMs)
    # MSM for particle type0
    MSMtype = 0
    markovModel0 = ctmsm(MSMtype, -1*seed) # random variable seed
    D0list = np.array([1.0])
    Drot0list = np.array([1.0])
    markovModel0.setD(D0list)
    markovModel0.setDrot(Drot0list)
    # MSM for particle type1
    MSMtype = 1
    ratematrix = np.array([[-0.5,0.5],[6.0,-6.0]])
    markovModel1 = ctmsm(MSMtype, ratematrix, -2*seed)
    D1list = np.array([0.5, 1.0])
    Drot1list = np.array([0.5, 1.0])
    markovModel1.setD(D1list)
    markovModel1.setDrot(Drot1list)
    unboundMSMlist = [markovModel0, markovModel1]

    # Define boundary
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)


    # Define particle list
    particleList = particleTools.randomParticleMSList(Nparticles, boxsize, separationDistance, particleTypes,
                                                 unboundMSMlist, seed)

    # Defines patchy protein potential
    potentialPatchyProteinMS = patchyProteinMarkovSwitch(sigma, strength, angularStrength,
                                                         patchesCoordinates1, patchesCoordinates2)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
    integrator = odLangevinMS(unboundMSMlist, dt, seed, bodytype)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyProteinMS)

    # Creates simulation
    sim = msmrd2.simulation(integrator)

    # Output filename definition
    filename = basefilename + "_{:04d}".format(simnumber)

    # Runs simulation
    sim.run(particleList, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked, trajtype)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)
