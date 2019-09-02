import numpy as np
import msmrd2
import random
from msmrd2.potentials import patchyParticleAngular
from msmrd2.integrators import overdampedLangevin as odLangevin
import msmrd2.tools.quaternions as quaternions
import multiprocessing
from multiprocessing import Pool
import os

'''
Creates an MD simulation of two particle and calculates their first passage times (FPTs) from a random 
bound state to any unbound state. The data is written to '../data/dimer/first_passage_times/filename_here.
'''

# Main parameters for particle and integrator
numparticles = 2
dt = 0.0001 #0.00001 #0.000005
bodytype = 'rigidbody'
D = 1.0
Drot = 1.0
initialState = 'B'
targetState = 'A'
radialBounds = [1.25, 2.25] # must match patchyDimer discretization
minimumUnboundRadius = 2.5
numTrajectories = 10000
# Other important parameters
boxsize = 6

# Parameters of patchy Particle potential with angular dependence (make sure it is consistent with msmrd data)
sigma = 1.0
strength = 60.0
angularStrength = 10.0
angleDiff = 3*np.pi/5.0
patch1 = np.array([np.cos(angleDiff/2),np.sin(angleDiff/2),0.])
patch2 = np.array([np.cos(-angleDiff/2),np.sin(-angleDiff/2),0.])
patchesCoordinates = [patch1, patch2]

# Bound states, needed to calculate boundstate of patchydimer
boundStatesA = [1, 2, 5, 6] # U-shaped bound dimer, corresponds to A state
boundStatesB = [3, 4, 7, 8] # Zigzag-shaped bound dimer, corresponds to B state

# Chooses parent directory
parentDirectory = "../../data/dimer/first_passage_times/"
# parentDirectory = "/group/ag_cmb/scratch/maojrs/msmrd2_data/dimer/first_passage_times/"

# Creates parent directory if it doesn't exist already
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("First passage times directory already exists. Simulation continues.")

# Create empty files to save the data in parallel algorithm
filename = parentDirectory + 'patchyDimerFPTs_' + initialState + '2' + targetState + '_trajs' \
           + str(numTrajectories) +'_boxsize' + str(boxsize) + '.xyz'


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




def simulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state a bound state
    :param trajectorynum: number of trajectories on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define dummy trajectory to extract bound states from python (needed to use getState function)
    dummyTraj = msmrd2.trajectories.patchyDimer(2, 1, radialBounds[0], minimumUnboundRadius)

    # Define simulation boundaries (choose either spherical or box)
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, 'periodic')

    # Define potential
    potentialPatchyParticleAngular = patchyParticleAngular(sigma, strength, angularStrength, patchesCoordinates)

    # Define integrator and boundary (over-damped Langevin)
    seed = int(-1*trajectorynum) # Negative seed, uses random device as seed
    integrator = odLangevin(dt, seed, bodytype)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyParticleAngular)

    # Generate random position and orientation particle list with two particles
    seed = int(trajectorynum)
    partlist = generateParticleList(initialState, boxsize, D, Drot, seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # the target bound state is reached. The output in the files is the elapsed time.
    intransition = True
    while(intransition):
        integrator.integrate(partlist)
        currentState = dummyTraj.getState(partlist[0], partlist[1])
        if targetState == 'A' and currentState in boundStatesA:
            intransition = False
            return 'A', integrator.clock
        elif targetState == 'B' and currentState in boundStatesB:
            intransition = False
            return 'B', integrator.clock
        elif integrator.clock >= 2000.0:
            intransition = False
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




