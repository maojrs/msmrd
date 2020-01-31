import numpy as np
import msmrd2
from msmrd2.potentials import patchyProteinMarkovSwitch
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.integrators import overdampedLangevinMarkovSwitch as odLangevinMS
import msmrd2.tools.particleTools as particleTools
from multiprocessing import Pool
import os

'''
Creates an MD simulation of two patchyProtein particle and calculates their first passage times (FPTs) from a random 
unbound configuration to any given bound state(s). The data is written to 
'../data/patchyProtein/first_passage_times/filename_here.
'''

# Main parameters for particle and integrator
numparticles = 2
dt = 0.0001 #0.0001 #0.00001 #0.000005
bodytype = 'rigidbody'
particleTypes = [0, 1]
minimumUnboundRadius = 2.5
numTrajectories = 600 #10000

# Define Patchy Protein potential parameters (This values are fixed and should match
# those used to determine metastable states in potential and trajectory.)
sigma = 1.0
strength = 65
angularStrength = 2
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
filename = parentDirectory  + 'patchyProteinFPTs_2bound_trajs' + str(numTrajectories) \
           +'_boxsize' + str(boxsize) + '.xyz'

def simulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state a bound state
    :param trajectorynum: number of trajectories on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define simulation boundaries (choose either spherical or box)
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, 'periodic')


    # Define dummy trajectory to extract bound states from python (needed to use getState function)
    radialLowerBound = 1.25
    radialUpperBound = 2.25
    dummyTraj = msmrd2.trajectories.patchyProtein(numparticles,1024, radialLowerBound, radialUpperBound)
    dummyTraj.setBoundary(boxBoundary)

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


    # Define potential
    potentialPatchyProtein = patchyProteinMarkovSwitch(sigma, strength, angularStrength,
                                                       patchesCoordinates1, patchesCoordinates2)

    # Define integrator and boundary (over-damped Langevin)
    #integrator = odLangevin(dt, seed, bodytype)
    integrator = odLangevinMS(unboundMSMlist, dt, seed, bodytype)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyProtein)

    # Creates random particle list
    seed = int(trajectorynum)
    partlist = particleTools.randomParticleMSList(numparticles, boxsize, minimumUnboundRadius,
                                                  particleTypes, unboundMSMlist, seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    while(unbound):
        integrator.integrate(partlist)
        boundState = dummyTraj.getState(partlist[0], partlist[1])
        if boundState in boundStates:
            unbound = False
            return boundState, integrator.clock
        elif integrator.clock >= 15000.0:
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
        for index, result in enumerate(pool.imap(simulationFPT, trajNumList)):
            state, time = result
            if state in boundStates:
                file.write(str(state) + ' ' + str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
                print(str(time))
            else:
                print("Simulation " + str(index) + ", done. Failed :(")


# Run parallel code
multiprocessingHandler()

# Serial code for testing with gdb
# trajNumList = list(range(numTrajectories))
# with open(filename, 'w') as file:
#     for index in range(numTrajectories):
#         state, time = simulationFPT(index)
#         if state in boundStates:
#             file.write(str(state) + ' ' + str(time) + '\n')
#             print("Simulation " + str(index) + ", done. Success!")
#         else:
#             print("Simulation " + str(index) + ", done. Failed :(")


