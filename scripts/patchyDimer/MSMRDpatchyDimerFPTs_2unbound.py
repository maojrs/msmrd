import numpy as np
import pickle
import multiprocessing
from multiprocessing import Pool
import msmrd2
import random
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.markovModels import msmrdMarkovModel as msmrdMSM
from msmrd2.integrators import msmrdIntegrator
import msmrd2.tools.quaternions as quaternions
import os

'''
Creates an MSM/RD simulation of two particle that calculates first passage times (FPTs) from a random 
bound state (A or B) to any unbound configuration. The data is written 
to '../../data/patchyDimer/first_passage_times/MSMRDfilename_here.
'''

# Main parameters for particle and integrator
numParticles = 2
partTypes = 0 # All particles are type 0
dt = 0.0001 #0.002 # should be smaller than Gillespie inverse transition rates
bodytype = 'rigidbody'
initialState = 'A'
numBoundStates = 8
maxNumBoundStates = 10
radialBounds = [1.25, 2.25] # must match patchyDimer discretization
minimumUnboundRadius = 2.5
numParticleTypes = 1 # num. of particle types (not states) in unbound state
numTrajectories = 5000

# Other important parameters
lagtime = 150 #150 #300
boxsize = 6
angleDiff = 3*np.pi/5.0
dtMDsimulation = 0.00001
stride = 25
realLagtime = lagtime*dtMDsimulation*stride

# Discretization parameters (need to be consistent with the on used to generate the rate dictionary
numSphericalSectionsPos = 7 #7 #7
numRadialSectionsQuat = 5 #3 #5
numSphericalSectionsQuat = 7 #6 #7
totalnumSecsQuat = numSphericalSectionsQuat*(numRadialSectionsQuat -1) + 1
numTransitionsStates = numSphericalSectionsPos * totalnumSecsQuat #203

# Parameters to define continuous-time MSM for unbound dynamics: unboundMSM (assumed same for all particles)
MSMtype = 0
ratematrix = np.array([[0]]) # no unbound dynamics
Dlist = np.array([1.0])
Drotlist = np.array([1.0])

# Parameters to define coupling Markov model for bound dynamics: couplingMSM
Dbound = 0.5*np.ones(numBoundStates)
DboundRot = np.ones(numBoundStates)

# Bound states definition, needed to calculate boundstate
boundStatesA = [1, 2, 5, 6] # U-shaped bound dimer, corresponds to A state
boundStatesB = [3, 4, 7, 8] # Zigzag-shaped bound dimer, corresponds to B state

# Chooses parent directory
parentDirectory = "../../data/patchyDimer/first_passage_times/"
MSMdirectory = '../../data/patchyDimer/MSMs/'

# Creates parent directory if it doesn't exist already
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("First passage times directory already exists. Simulation continues.")

# Chooses filename for output file with the results of the parallel simulation
filename = parentDirectory + 'MSMRDpatchyDimerFPTs_' + initialState + '2unbound_trajs' \
           + str(numTrajectories) + '_lagt' + str(lagtime) + '_boxsize' + str(boxsize) + '.xyz'

def generateParticleList(state, boxsize, types, couplingMSM, randomSeed = -1):
    '''
    Generate a random particle list of two particles corresponding to either state A or state B. As each
    state has several different configurations, this functions picks one randomly. The definition of the
    states is taken form patchyDimer.cpp
    :param state: bound state in which the two particle list should be in, A or B
    :param boxsize: size of simulation box, if scalar it assumes the three box edges are the same in all dimensions
    :param types: array containing the particle types. The index should correspond to that of of the particle list.
    :param couplingMSM: MSM for MSM/RD, needed to extract diffusion coefficients of particles.
    :param randomSeed: seed for python random generator. Important to specify in parallel runs. Default value of -1
    will use the default seed.
    :return: random particle list corresponding to either state A or state B.
    '''

    # Transform boxsize and type to vector if neccesary.
    boxsize = boxsize - 1 # to ensure the whole compund is initially inside box (1 is the norm of relpos below)
    if np.isscalar(boxsize):
        boxsize = np.array([boxsize, boxsize, boxsize])

    if np.isscalar(types):
        newtypes = np.ones(2, dtype = 'int')
        types = (types*newtypes).astype(int)

    position1 = np.array([0.0, 0.0, 0.0])
    #position1 = np.array([boxsize[0]*random.random()-0.5*boxsize[0],
    #                    boxsize[1]*random.random()-0.5*boxsize[1],
    #                     boxsize[2]*random.random()-0.5*boxsize[2]])
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
    substate = 0
    if state == 'A':
        substate = random.choice([1,2,5,6])
        position2 = pos2[substate - 1]
        orientation2 = quatRotations[substate - 1]
    if state == 'B':
        substate = random.choice([3,4,7,8])
        position2 = pos2[substate - 1]
        orientation2 = quatRotations[substate - 1]
    # Obtain diffusion coefficients from unboundMSM
    Dbound = couplingMSM.D[substate - 1]
    DrotBound = couplingMSM.Drot[substate - 1]
    # Define particles
    part1 = msmrd2.particle(types[0], -1, Dbound, DrotBound, position1, orientation1)
    part2 = msmrd2.particle(types[1], -1, 0., 0., position2, orientation2)
    # Set up bound particles as if done by the code (deactivate one of them)
    part1.setBoundState(substate)
    part2.setBoundState(substate)
    part1.setBoundTo(1)
    part2.setBoundTo(0)
    part1.deactivateMSM()
    part2.deactivateMSM()
    part2.deactivate()
    part2.setPosition([10000000.0, 10000000.0, 10000000.0])
    partlist = msmrd2.integrators.particleList([part1, part2])
    return partlist

def MSMRDsimulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state a bound state
    :param trajectorynum: trajectory number (index on parallel computation) on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define dummy trajectory to extract bound states from python (needed to use getState function)
    dummyTraj = msmrd2.trajectories.patchyDimer(2, 1, radialBounds[0], minimumUnboundRadius)

    # Define discretization
    discretization = msmrd2.discretizations.positionOrientationPartition(radialBounds[1],
                                            numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat)

    # Define boundary
    boxBoundary = msmrd2.box(boxsize,boxsize,boxsize,'periodic')

    # Load rate dicitionary
    pickle_in = open(MSMdirectory + "MSM_patchyDimer_t3.00E+06_s25_lagt" + str(lagtime)
                     +  ".pickle","rb")
    mainMSM = pickle.load(pickle_in)
    tmatrix = mainMSM['transition_matrix']
    activeSet = mainMSM['active_set']

    # Set unbound MSM
    seed = int(-2*trajectorynum) # Negative seed, uses random device as seed
    unboundMSM = ctmsm(MSMtype, ratematrix, seed)
    unboundMSM.setD(Dlist)
    unboundMSM.setDrot(Drotlist)

    # Set coupling MSM
    seed = int(-3*trajectorynum) # Negative seed, uses random device as seed
    couplingMSM = msmrdMSM(numBoundStates, maxNumBoundStates,  tmatrix, activeSet, realLagtime, seed)
    couplingMSM.setDbound(Dbound, DboundRot)

    # Define integrator, boundary and discretization
    seed = -int(1*trajectorynum) # Negative seed, uses random device as seed
    integrator = msmrdIntegrator(dt, seed, bodytype, numParticleTypes, radialBounds, unboundMSM, couplingMSM)
    integrator.setBoundary(boxBoundary)
    integrator.setDiscretization(discretization)

    # Generate random position and orientation particle list with two particles
    seed = int(trajectorynum)
    partlist = generateParticleList(initialState, boxsize, partTypes, couplingMSM, seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # an unbound state is reached. The output in the files is the elapsed time.
    bound = True
    while(bound):
        integrator.integrate(partlist)
        #currentState = partlist[0].state
        currentState = dummyTraj.getState(partlist[0], partlist[1])
        if (partlist[0].boundTo == -1) and (currentState == 0):
            bound = False
            return initialState, integrator.clock
        elif integrator.clock >= 10000.0:
            unbound = False
            return 'Failed at:', integrator.clock



def multiprocessingHandler():
    '''
    Handles parallel processing of simulationFPT and writes to same file in parallel
    '''
    num_cores = multiprocessing.cpu_count()
    pool = Pool(processes=num_cores)
    trajNumList = list(range(numTrajectories))
    with open(filename, 'w') as file:
        for index, result in enumerate(pool.imap(MSMRDsimulationFPT, trajNumList)):
            state, time = result
            if state == 'A' or state == 'B':
                file.write(state + ' ' + str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
            else:
                print("Simulation " + str(index) + ", done. Failed :(")


# Run parallel code
multiprocessingHandler()

# Serial code for testing with gdb
# trajNumList = list(range(numTrajectories))
# with open(filename, 'w') as file:
#     for index in range(numTrajectories):
#         state, time = MSMRDsimulationFPT(index)
#         if state == 'A' or state == 'B':
#             file.write(state + ' ' + str(time) + '\n')
#             print("Simulation " + str(index) + ", done. Success!")
#         else:
#             print("Simulation " + str(index) + ", done. Failed :(")
