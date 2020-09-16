import numpy as np
import pickle
import multiprocessing
from multiprocessing import Pool
import msmrd2
from msmrd2.markovModels import continuousTimeMarkovStateModel as ctmsm
from msmrd2.markovModels import msmrdMarkovModel as msmrdMSM
from msmrd2.integrators import msmrdMultiParticleIntegrator
import msmrd2.tools.particleTools as particleTools
import os

'''
Creates an MD simulation of five particles and calculates their first passage times (FPTs) 
from a random configuration to a pentamer configuration following the MSMRD algorithm. The data is written to 
'../data/pentamer/first_passage_times/filename_here. It uses the patchy dimer model of molecules
with only one stable angular configuration.
'''

# Main parameters for particle and integrator
numParticles = 5
partTypes = 0 # All particles are type 0
dt = 0.0001 #0.002 # should be smaller than Gillespie inverse transition rates
bodytype = 'rigidbody'
numBoundStates = 4
maxNumBoundStates = 10
radialBounds = [1.25, 2.25] # must match patchyDimer discretization
minimumUnboundRadius = 1.5
numParticleTypes = 1 # num. of particle types (not states) in unbound state
numTrajectories = 6000

# Other important parameters
lagtime = 40 #100 #300
boxsize = 6 #2.5 #6 #8 #6
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
boundStates = [1, 2, 3, 4]

# Chooses parent directory
parentDirectory = "../../data/pentamer/first_passage_times/"
# parentDirectory = "/group/ag_cmb/scratch/maojrs/msmrd2_data/pentamer/first_passage_times/"

# Creates parent directory if it doesn't exist already
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("First passage times directory already exists. Simulation continues.")

# Chooses filename for output file with the results of the parallel simulation
filename = parentDirectory + 'MSMRDpentamerFPTs_trajs' + str(numTrajectories) + \
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

    # Define boundary
    boxBoundary = msmrd2.box(boxsize,boxsize,boxsize,'periodic')

    # Load rate dicitionary
    pickle_in = open("../../data/pentamer/MSMs/MSM_dimer4trimer_t3.00E+06_s25_lagt" + str(lagtime)
                     +  ".pickle","rb") # Same MSM as trimer
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
    integrator = msmrdMultiParticleIntegrator(dt, seed, bodytype, numParticleTypes, radialBounds, unboundMSM, couplingMSM)
    integrator.setBoundary(boxBoundary)
    integrator.setDiscretization(discretization)

    # Generate random position and orientation particle list with two particles
    seed = int(trajectorynum)
    partlist = particleTools.randomParticleMSList(numParticles, boxsize,
                                                  minimumUnboundRadius, partTypes, [unboundMSM], seed)

    # Calculates the first passage times to a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    while(unbound):
        integrator.integrate(partlist)
        #if (partlist[0].compoundIndex > -1):
        if (partlist[0].compoundIndex != -1): #FIND SEGFAULT, PROBABLY IN joinParticleCompounds FNCTION
            #numBindings = integrator.getNumberOfBindingsInCompound(partlist[0].compoundIndex)
            numBindings = integrator.getCompoundSize(partlist[0].compoundIndex)
            print(integrator.clock, partlist[0].compoundIndex, numBindings)
            if (numBindings >= 5 and numBindings <= 10):
                unbound = False
                return "pentamer", integrator.clock
        elif integrator.clock >= 1000.0:
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
            if state == 'pentamer':
                file.write(state + ' ' + str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
            else:
                print("Simulation " + str(index) + ", done. Failed :(")


# Run parallel code
#multiprocessingHandler()

# Serial code for testing with gdb
with open(filename, 'w') as file:
    for index in range(numTrajectories):
        state, time = MSMRDsimulationFPT(index)
        if state == 'pentamer':
            file.write(state + ' ' + str(time) + '\n')
            print("Simulation " + str(index) + ", done. Success!")
        else:
            print("Simulation " + str(index) + ", done. Failed :(")
