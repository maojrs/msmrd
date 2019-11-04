import msmrd2
import numpy as np
from msmrd2.integrators import overdampedLangevin as odLangevin
from msmrd2.potentials import patchyParticleAngular
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
D = 1
Drot = 1
separationDistance = 1.25 # minimum separation distance for initial condition
numSimulations = 600 #500 #1000 #20 #200


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
strength = 60.0 #50.0
angularStrength = 10.0 #10.0
angleDiff = 3*np.pi/5.0
patch1 = np.array([np.cos(angleDiff/2),np.sin(angleDiff/2),0.])
patch2 = np.array([np.cos(-angleDiff/2),np.sin(-angleDiff/2),0.])
patchesCoordinates = [patch1, patch2]
potentialPatchyParticleAngular = patchyParticleAngular(sigma, strength, angularStrength, patchesCoordinates)

# Define simulaion boundaries (choose either spherical or box)
boxsize = 5
boundaryType = 'periodic'
boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)

# Parent directory location
parentDirectory = "../../data/dimer/"
# parentDirectory = "/group/ag_cmb/scratch/maojrs/msmrd2_data/dimer/"

# Create folder for data
foldername = "60strength"
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
parameterDictionary = {'numFiles' : numSimulations, 'numParticles' : Nparticles, 'dt' : dt, 'bodytype' : bodytype,
                       'D' : D, 'Drot' : Drot, 'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype,
                       'boxsize' : boxsize, 'boundaryType' : boundaryType, 'potentialStrength' : strength,
                       'potentialAngularStrength' : angularStrength, 'potentialPatchesAngleDiff' : angleDiff}
analysisTools.writeParameters(parameterfilename, parameterDictionary)


# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simDimer")



def generateParticleList(state, boxsize, D, Drot, randomSeed = -1):
    '''
    Generate a random particle list of two particles corresponding to either state A or state B. As each
    state has several different configurations, this functions picks one randomly. The definition of the
    states is taken form patchyDimer.cpp
    :param state: bound state in which the two particle list should be in, A or B
    :param D and Drot: diffusion coefficients of particles.
    :param randomSeed: seed for python random generator. Important to specify in parallel runs. Default value of -1
    will use the default seed.
    :return: random particle list corresponding to either state A or state B.
    '''

    # Transform boxsize to vector if neccesarry.
    boxsize = boxsize - 1 # to ensure the whole compund is initially inside box (1 is the norm of relpos below)
    if np.isscalar(boxsize):
        boxsize = np.array([boxsize, boxsize, boxsize])

    position1 = np.array([0.0, 0.0, 0.0])
    #position1 = np.array([boxsize[0]*random.random()-0.5*boxsize[0],
    #                      boxsize[1]*random.random()-0.5*boxsize[1],
    #                      boxsize[2]*random.random()-0.5*boxsize[2]])
    orientation1 = np.array([1.0, 0.0, 0.0, 0.0])
    # Define relative position
    relpos1 = np.array([np.cos(angleDiff / 2.0),np.sin(angleDiff / 2.0), 0])
    relpos2 = np.array([np.cos(angleDiff / 2.0),np.sin(-1*angleDiff / 2.0), 0])
    relPos1orthogonal = np.array([-1.0 * np.sin(angleDiff / 2.0), np.cos(angleDiff / 2.0), 0.0])
    relPos2orthogonal = np.array([np.sin(angleDiff / 2.0), np.cos(angleDiff / 2.0), 0.0])
    # Define relative rotations
    rotations = [None]*8
    rotations[0] = np.pi * relPos1orthogonal
    rotations[1] = np.array([0.0, 0.0, -2 * np.pi / 5.0])
    rotations[2] = np.array([0.0, 0.0, np.pi])
    rotations[3] = np.array([0.0, np.pi, 0.0])
    # --first 4 rotations correspond to binding on top patch of particle 1, next 4 rotations to bottom patch
    rotations[4] = np.pi * relPos2orthogonal
    rotations[5] = np.array([0.0, 0.0, 2 * np.pi / 5.0])
    rotations[6] = np.array([0.0, 0.0, np.pi])
    rotations[7] = np.array([0.0, np.pi, 0.0])
    # Convert axis-angle rotations to quaternions
    quatRotations = [None]*8
    for i in range(8):
        quatRotations[i] = quaternions.angle2quat(rotations[i])
    # Assign correct position2 depending on the states
    pos2 = [None]*8
    pos2[0:4] = [position1 + relpos1]*4
    pos2[4:8] = [position1 + relpos2]*4

    # Assign position and orientation depending on initial state A or B
    if state == 'A':
        substate = random.choice([1,2,5,6])
        position2 = pos2[substate - 1]
        orientation2 = quatRotations[substate - 1]
    if state == 'B':
        substate = random.choice([3,4,7,8])
        position2 = pos2[substate - 1]
        orientation2 = quatRotations[substate - 1]
    part1 = msmrd2.particle(D, Drot, position1, orientation1)
    part2 = msmrd2.particle(D, Drot, position2, orientation2)
    partlist = msmrd2.integrators.particleList([part1, part2])
    return partlist

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list, sometimes begin in unbound state (or transition region) and sometimes in state A or B.
    seed = int(simnumber)
    random.seed(seed)
    choice = random.randint(1,3)
    if choice == 1:
        partlist = particleTools.randomParticleList(Nparticles, boxsize, separationDistance, D, Drot, seed)
    if choice == 2:
        partlist = generateParticleList('A', boxsize, D, Drot, seed)
    if choice == 3:
        partlist = generateParticleList('B', boxsize, D, Drot, seed)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
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
