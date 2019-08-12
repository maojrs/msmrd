import numpy as np
import pickle
import multiprocessing
from multiprocessing import Pool
import msmrd2
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.markovModels import msmrdMarkovStateModel as msmrdMSM
from msmrd2.integrators import msmrdIntegrator
import msmrd2.tools.particleTools as particleTools

'''
Creates an MSM/RD simulation of two particle that calculates first passage times (FPTs) from a random 
unbound configuration to either of the two bound states (A or B). This requires as input a rate dictionary 
calculated from an MSM of an MD simulation (the rate dictionary is loaded using pickle). The data is 
written to '../data/dimer/first_passage_times/MSMRDfilename_here.
'''

# Main parameters for particle and integrator
numParticles = 2
partTypes = 0
dt = 0.0001 #0.002 # should be smaller than Gillespie inverse transition rates
bodytype = 'rigidbody'
numBoundStates = 8
maxNumBoundStates = 10
radialBounds = [1.4, 2.2] #2.2
relativeDistanceCutOff = radialBounds[1]
numParticleTypes = 1 # num. of particle types (not states) in unbound state
numTrajectories = 10000
# Other important parameters
lagtime = 600
boxsize = 6

# Discretization parameters (need to be consistent with the on used to generate the rate dictionary
numSphericalSectionsPos = 7
numRadialSectionsQuat = 5
numSphericalSectionsQuat = 7
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

# Create empty files to save the data in parallel algorithm
filename = '../../data/dimer/first_passage_times/MSMRDpatchyDimerFPTs_2bound_trajs' + str(numTrajectories) + \
           '_lagt' + str(lagtime) + '_boxsize' + str(boxsize) + '.xyz'

def MSMRDsimulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state a bound state
    :param trajectorynum: trajectory number (index on parallel computation) on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define discretization
    discretization = msmrd2.discretizations.positionOrientationPartition(relativeDistanceCutOff,
                                            numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat)

    # Define boundary
    boxBoundary = msmrd2.box(boxsize,boxsize,boxsize,'periodic')

    # Load rate dicitionary
    pickle_in = open("../../examples/pickled_data/ratedictionary_dimer_t2.00E+06_s25_lagt" + str(lagtime)
                     +  ".pickle","rb")
    rateDictionary = pickle.load(pickle_in)

    # Set unbound MSM
    seed = int(-2*trajectorynum) # Negative seed, uses random device as seed
    unboundMSM = ctmsm(MSMtype, ratematrix, seed)
    unboundMSM.setD(Dlist)
    unboundMSM.setDrot(Drotlist)

    # Set coupling MSM
    seed = int(-3*trajectorynum) # Negative seed, uses random device as seed
    couplingMSM = msmrdMSM(numBoundStates, numTransitionsStates, seed, rateDictionary)
    couplingMSM.setDbound(Dbound, DboundRot)
    couplingMSM.setMaxNumberBoundStates(maxNumBoundStates)

    # Define integrator, boundary and discretization
    seed = -int(1*trajectorynum) # Negative seed, uses random device as seed
    integrator = msmrdIntegrator(dt, seed, bodytype, numParticleTypes, radialBounds, unboundMSM, couplingMSM)
    integrator.setBoundary(boxBoundary)
    integrator.setDiscretization(discretization)

    # Generate random position and orientation particle list with two particles
    seed = int(trajectorynum)
    partlist = particleTools.randomParticleMSList(numParticles, boxsize,
                                                  relativeDistanceCutOff, partTypes, [unboundMSM], seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    while(unbound):
        integrator.integrate(partlist)
        currentState = partlist[0].state
        if currentState in boundStatesA:
            unbound = False
            return 'A', integrator.clock
        elif currentState in boundStatesB:
            unbound = False
            return 'B', integrator.clock
        elif integrator.clock >= 2000.0:
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
