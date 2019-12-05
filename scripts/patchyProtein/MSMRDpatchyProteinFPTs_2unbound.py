import numpy as np
import pickle
import multiprocessing
from multiprocessing import Pool
import msmrd2
from msmrd2.potentials import patchyProteinMarkovSwitch
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.markovModels import msmrdMarkovModel as msmrdMSM
from msmrd2.integrators import msmrdPatchyProtein
import msmrd2.tools.particleTools as particleTools
import msmrd2.tools.quaternions as quaternions
import os

'''
Creates an MSM/RD simulation of two particle that calculates first passage times (FPTs) from a 
given bound state to any unbound configuration. This requires to input the MSM for the MSM/RD algorithm 
calculated from an MD simulation (the MSM is loaded using pickle). The data is 
written to '../data/patchyProtein/first_passage_times/MSMRDfilename_here.
'''

# Main parameters for particle and integrator
numparticles = 2
dt = 0.0001 #0.00001 #0.000005
bodytype = 'rigidbody'
particleTypes = [0, 1]
numBoundStates = 6
maxNumBoundStates = 10
initialState = 1 # 1 to 6 possible values

radialBounds = [1.25, 2.25] # must match patchyProtein discretization trajectory
minimumUnboundRadius = 2.5
numParticleTypes = 2 # num. of particle types (not states) in unbound state
numTrajectories = 5000 #10000

# Other important parameters
lagtime = 150 #50 #75 #300
boxsize = 6 #8 #6
dtMDsimulation = 0.00001
stride = 100
realLagtime = lagtime*dtMDsimulation*stride

# Discretization parameters (need to be consistent with the on used to generate the rate dictionary
numSphericalSectionsPos = 6
numRadialSectionsQuat = 4
numSphericalSectionsQuat = 6
totalnumSecsQuat = numSphericalSectionsQuat*(numRadialSectionsQuat - 1) + 1
numTransitionsStates = numSphericalSectionsPos * totalnumSecsQuat #228

# Define Patchy Protein potential parameters (This values are fixed and should match
# those used to determine metastable states in potential and trajectory.) Note here, we use
# no patches since we only want the isotropic attractive part.
sigma = 1.0
strength = 65
angularStrength = 2
patchesCoordinates1 = []
patchesCoordinates2 = []

# Set diffusion coefficients for bound states
# Parameters to define coupling Markov model for bound dynamics: couplingMSM
Dbound = 0.5*np.ones(numBoundStates)
DboundRot = np.ones(numBoundStates)

# Chooses parent directory
parentDirectory = "../../data/patchyProtein/first_passage_times/"
# parentDirectory = "/group/ag_cmb/scratch/maojrs/msmrd2_data/dimer/first_passage_times/"

# Creates parent directory if it doesn't exist already
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("First passage times directory already exists. Simulation continues.")

# Chooses filename for output file with the results of the parallel simulation
filename = parentDirectory + 'MSMRDpatchyProteinFPTs_2unbound_trajs' + str(numTrajectories) + \
           '_lagt' + str(lagtime) + '_boxsize' + str(boxsize) + '.xyz'


def generateParticleList(state, boxsize, types, couplingMSM, randomSeed = -1):
    '''
    Generate a random particle list of two particles corresponding to an initial unbound state. As each
    state has several different configurations, this functions picks one randomly. The definition of the
    states is taken form patchyProtein.cpp (trajectory)
    :param state: bound state in which the two particle list should begin in (can be a list of states)
    :param boxsize: size of simulation box, if scalar it assumes the three box edges are the same in all dimensions
    :param types: array containing the particle types. The index should correspond to that of of the particle list.
    :param couplingMSM: MSM for MSM/RD, needed to extract diffusion coefficients of particles.
    :param randomSeed: seed for python random generator. Important to specify in parallel runs. Default value of -1
    will use the default seed.
    :return: random particle list corresponding to one of the initialBoundstates.
    '''

    # Transform boxsize to vector if neccesarry.
    boxsize = boxsize - 1 # to ensure the whole compund is initially inside box (1 is the norm of relpos below)
    if np.isscalar(boxsize):
        boxsize = np.array([boxsize, boxsize, boxsize])

    # Transform state to list if neccesary
    if np.isscalar(state):
        state = [state]

    position1 = np.array([0.0, 0.0, 0.0])
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

    # Assign position and orientation depending on initial bound state (between 0 and 5)
    substate = np.random.choice(state) # if states was originally a scalar, always picks that value.
    position2 = relPos[substate - 1]
    orientation2 = quatRotations[substate - 1]
    # Obtain diffusion coefficients from unboundMSM
    Dbound = couplingMSM.D[substate - 1]
    DrotBound = couplingMSM.Drot[substate - 1]
    # Define particles
    part1 = msmrd2.particle(0, -1, Dbound, DrotBound, position1, orientation1)
    part2 = msmrd2.particle(1, -1, 0.0, 0.0, position2, orientation2)
    # Define properties of bound particles
    part1.setBoundState(substate)
    part2.setBoundState(substate)
    part1.setBoundTo(1)
    part2.setBoundTo(0)
    part1.deactivateMSM()
    part2.deactivateMSM()
    part2.deactivate()
    part2.setPosition([10000000.0, 10000000.0, 10000000.0])
    # Define and return particle list
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
    dummyTraj = msmrd2.trajectories.patchyProtein(numparticles, 1, radialBounds[0], minimumUnboundRadius)

    # Define discretization
    discretization = msmrd2.discretizations.positionOrientationPartition(radialBounds[1],
                    numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat)
    discretization.setThetasOffset(np.pi/4.0);

    # Define boundary
    boxBoundary = msmrd2.box(boxsize,boxsize,boxsize,'periodic')

    # Define isotropic potential (no patches)
    potentialPatchyProtein = patchyProteinMarkovSwitch(sigma, strength, angularStrength,
                                                       patchesCoordinates1, patchesCoordinates2)

    # Define base seed
    seed = int(-1*trajectorynum) # Negative seed, uses random device as seed

    # Load rate dicitionary
    pickle_in = open("../../data/pickled_data/MSM_patchyProtein_t1.20E+07_s100_lagt" + str(lagtime)
                     +  ".pickle","rb")
    mainMSM = pickle.load(pickle_in)
    tmatrix = mainMSM['transition_matrix']
    activeSet = mainMSM['active_set']

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

    # Set coupling MSM for MSM/RD
    seed = int(-3*trajectorynum) # Negative seed, uses random device as seed
    couplingMSM = msmrdMSM(numBoundStates, maxNumBoundStates,  tmatrix, activeSet, realLagtime, seed)
    couplingMSM.setDbound(Dbound, DboundRot)

    # Define integrator, boundary and discretization
    seed = -int(1*trajectorynum) # Negative seed, uses random device as seed
    integrator = msmrdPatchyProtein(dt, seed, bodytype, numParticleTypes, radialBounds, unboundMSMlist, couplingMSM)
    integrator.setBoundary(boxBoundary)
    integrator.setDiscretization(discretization)
    integrator.setPairPotential(potentialPatchyProtein)

# Creates random particle list
    seed = int(trajectorynum)
    partlist = generateParticleList(initialState, boxsize, particleTypes, couplingMSM, seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    bound = True
    while(bound):
        integrator.integrate(partlist)
        currentState = dummyTraj.getState(partlist[0], partlist[1])
        if (partlist[0].boundTo == -1) and (currentState == 0):
            bound = False
            return initialState, integrator.clock
        elif integrator.clock >= 15000.0:
            bound = False
            return 'failed', integrator.clock



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
            if state != 'failed':
                file.write(str(state) + ' ' + str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
            else:
                print("Simulation " + str(index) + ", done. Failed :(")


# Run parallel code
multiprocessingHandler()

# # Serial code for testing with gdb
# with open(filename, 'w') as file:
#     for index in range(numTrajectories):
#         state, time = MSMRDsimulationFPT(index)
#         if state in boundStates:
#             file.write(str(state) + ' ' + str(time) + '\n')
#             print("Simulation " + str(index) + ", done. Success!")
#         else:
#             print("Simulation " + str(index) + ", done. Failed :(")
