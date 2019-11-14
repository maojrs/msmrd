import numpy as np
import msmrd2
import random
from msmrd2.potentials import patchyProteinMarkovSwitch
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.integrators import overdampedLangevinMarkovSwitch as odLangevinMS
import msmrd2.tools.quaternions as quaternions
import multiprocessing
from multiprocessing import Pool
import os

'''
Creates an MD simulation of two particle and calculates their first passage times (FPTs) from a given 
bound state to any unbound state. The data is written to '../data/patchyProtein/first_passage_times/filename_here.
'''

# Main parameters for particle and integrator
numparticles = 2
dt = 0.0001 #0.00001 #0.000005
bodytype = 'rigidbody'
particleTypes = [0, 1]
minimumUnboundRadius = 2.5
numTrajectories = 20 #10000
numBoundStates = 6
initialState = 1

# Define Patchy Protein potential parameters (This values are fixed and should match
# those used to determine metastable states in potential and trajectory.)
sigma = 1.0
strength = 65
patchesCoordinates1 = [np.array([1.,0.,0.]), np.array([0.,1.,0.]), np.array([0.,0.,1.]),
                       np.array([-1.,0.,0.]), np.array([0.,-1.,0.]), np.array([0.,0.,-1.])]
patchesCoordinates2 = [np.array([1.,0.,0.])]

# Define simulation boundaries parameters (choose either spherical or box)
boxsize = 6
boundaryType = 'periodic'

# Define bound states
#boundStates = [1,2,3,4,5,6]
boundStates = [1] # Main bound state


# Chooses parent directory
parentDirectory = "../../data/patchyProtein/first_passage_times/"
# parentDirectory = "/group/ag_cmb/scratch/maojrs/msmrd2_data/patchyProtein/first_passage_times/"

# Creates parent directory if it doesn't exist already
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("First passage times directory already exists. Simulation continues.")

# Create empty files to save the data in parallel algorithm
filename = parentDirectory + 'patchyProteinFPTs_' + initialState + '_2unbound_trajs' \
           + str(numTrajectories) +'_boxsize' + str(boxsize) + '.xyz'


def generateParticleList(state, boxsize, D, Drot, randomSeed = -1):
    '''
    Generate a random particle list of two particles corresponding to either state A or state B. As each
    state has several different configurations, this functions picks one randomly. The definition of the
    states is taken form patchyProtein.cpp (trajectory)
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
    # Define relative positions
    relPos = [None]*numBoundStates
    relPos[0] = np.array([1., 0., 0.])
    relPos[1] = np.array([0., 1., 0.])
    relPos[2] = np.array([0., 0., 1.])
    relPos[3] = np.array([-1., 0., 0.])
    relPos[4] = np.array([0., -1., 0.])
    relPos[5] = np.array([0., 0., -1.])
    # Define relative rotations
    rotations = [None]*numBoundStates
    rotations[0] = np.array([0.0, 0.0, np.pi])
    rotations[1] = np.array([0.0, 0.0, -np.pi / 2.0])
    rotations[2] = np.array([0.0, np.pi / 2.0, 0.0])
    rotations[3] = np.array([0.0, 0.0, 0.0])
    rotations[4] = np.array([0.0, 0.0, np.pi / 2.0])
    rotations[5] = np.array([0.0, -np.pi / 2.0, 0.0])
    # Convert axis-angle rotations to quaternions
    quatRotations = [None]*numBoundStates
    for i in range(numBoundStates):
        quatRotations[i] = quaternions.angle2quat(rotations[i])

    # Assign position and orientation depending on initial state A or B
    substate = initialState # Could be random among all of them, e.g. random.choice([1,2,3,4,5,6])
    position2 = relPos[substate - 1]
    orientation2 = quatRotations[substate - 1]
    part1 = msmrd2.particle(D, Drot, position1, orientation1)
    part2 = msmrd2.particle(D, Drot, position2, orientation2)
    partlist = msmrd2.integrators.particleList([part1, part2])
    return partlist




def simulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state a bound state
    :param trajectorynum: number of trajectories on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define dummy trajectory to extract bound states from python (needed to use getState function)
    dummyTraj = msmrd2.trajectories.patchyProtein(numparticles,1024)

    # Define base seed
    seed = int(-1*trajectorynum) # Negative seed, uses random device as seed

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

    # Define simulation boundaries (choose either spherical or box)
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, 'periodic')

    # Define potential
    potentialPatchyProtein = patchyProteinMarkovSwitch(sigma, strength, patchesCoordinates1, patchesCoordinates2)

    # Define integrator and boundary (over-damped Langevin)
    integrator = odLangevinMS(unboundMSMlist, dt, seed, bodytype)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyProtein)

    # Generate random position and orientation particle list with two particles
    seed = int(trajectorynum)
    partlist = generateParticleList(initialState, boxsize, D, Drot, seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    bound = True
    while(bound):
        integrator.integrate(partlist)
        boundState = dummyTraj.getState(partlist[0], partlist[1])
        if boundState == 0:
            bound = False
            return initialState, integrator.clock
        elif integrator.clock >= 10000.0:
            bound = False
            return 'Failed at:', integrator.clock



def multiprocessingHandler():
    '''
    Handles parallel processing of simulationFPT and writes to same file in parallel
    '''
    num_cores = multiprocessing.cpu_count()
    pool = Pool(processes=num_cores)
    trajNumList = list(range(numTrajectories))
    with open(filename, 'w') as file:
        for index, result in enumerate(pool.imap(simulationFPT, trajNumList)):
            state, time = result
            if state == 'A' or state == 'B':
                file.write(state + ' ' + str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
            else:
                print("Simulation " + str(index) + ", done. Failed :(")


# Run parallel code
multiprocessingHandler()




