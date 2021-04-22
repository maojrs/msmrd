import numpy as np
import pickle
import multiprocessing
from multiprocessing import Pool
import msmrd2
from msmrd2.potentials import patchyProteinMarkovSwitch
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.markovModels import msmrdMarkovModel as msmrdMSM
from msmrd2.integrators import msmrdPatchyProtein
#from msmrd2.integrators import msmrdPatchyProtein2
import msmrd2.tools.particleTools as particleTools
import os

'''
Creates an MSM/RD simulation of two particle that calculates first passage times (FPTs) from a random 
unbound configuration to a given bound state. This requires to input the MSM for the MSM/RD algorithm 
calculated from an MD simulation (the MSM is loaded using pickle). The data is 
written to '../../data/patchyProtein/first_passage_times/MSMRDfilename_here.
'''

# Main parameters for particle and integrator
numparticles = 2
dt = 0.0001 #0.0001 #0.00001 #0.000005
bodytype = 'rigidbody'
particleTypes = [0, 1]
numBoundStates = 6
maxNumBoundStates = 10
radialBounds = [1.25, 2.25] # must match patchyProtein discretization trajectory
minimumUnboundRadius = 2.5
numParticleTypes = 2 # num. of particle types (not states) in unbound state
numTrajectories = 5000 #600 #1000 #10000

# Other important parameters
lagtime = 150 #75 #75 #50 #75 #150 #75 #300
boxsize = 6 #8 #6
dtMDsimulation = 0.00001
stride = 25
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
strength = 60 #65
angularStrength = 10 #2
patchesCoordinates1 = []
patchesCoordinates2 = []

# Set diffusion coefficients for bound states
# Parameters to define coupling Markov model for bound dynamics: couplingMSM
Dbound = np.ones(numBoundStates)
DboundRot = np.ones(numBoundStates)

# Bound states definition, needed to calculate boundstate
boundStates = [1, 2, 3, 4, 5, 6]
#boundStates = [1]

# Chooses parent directory
parentDirectory = "../../data/patchyProtein/first_passage_times/"
MSMdirectory = '../../data/patchyProtein/MSMs/'

# Creates parent directory if it doesn't exist already
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("First passage times directory already exists. Simulation continues.")

# Chooses filename for output file with the results of the parallel simulation
filename = parentDirectory + 'MSMRDpatchyProteinFPTs_2bound_trajs' + str(numTrajectories) + \
           '_lagt' + str(lagtime) + '_boxsize' + str(boxsize) + '.xyz'

def MSMRDsimulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state a bound state
    :param trajectorynum: trajectory number (index on parallel computation) on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define discretization
    discretization = msmrd2.discretizations.positionOrientationPartition(radialBounds[1],
                                                                         numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat)
    discretization.setThetasOffset(np.pi/4.0)

    # Define boundary
    boxBoundary = msmrd2.box(boxsize,boxsize,boxsize,'periodic')

    # Define isotropic potential (no patches)
    potentialPatchyProtein = patchyProteinMarkovSwitch(sigma, strength, angularStrength,
                                                       patchesCoordinates1, patchesCoordinates2)


    # Define base seed
    seed = int(-1*trajectorynum) # Negative seed, uses random device as seed

    # Load rate dicitionary
    pickle_in = open(MSMdirectory + "MSM_patchyProtein_t6.00E+06_s25_lagt" + str(lagtime) +  ".pickle","rb")
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
    partlist = particleTools.randomParticleMSList(numparticles, boxsize, minimumUnboundRadius,
                                                  particleTypes, unboundMSMlist, seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    while(unbound):
        integrator.integrate(partlist)
        currentState = partlist[0].boundState
        if currentState in boundStates:
            unbound = False
            #filenameLog = filename = "../../data/patchyProtein/eventLog_" + str(trajectorynum)
            #integrator.printEventLog(filenameLog)
            return currentState, integrator.clock
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
        for index, result in enumerate(pool.imap(MSMRDsimulationFPT, trajNumList)):
            state, time = result
            if state in boundStates:
                file.write(str(state) + ' ' + str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
                print(str(time))
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

