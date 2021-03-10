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
'../../data/pentamer/first_passage_times/filename_here. It uses the patchy dimer model of molecules
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
numTrajectories = 20000 #4*6000
printAllRingFormations = False # if false, it only print the pentameric ring FPTs (not trimeric nor tetrameric)


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
Dbound = np.ones(numBoundStates)
DboundRot = 0.5*np.ones(numBoundStates)

# Complex diffusion coefficients (D,Drot) (taken from estimateDiffusionCoefficients script)
DlistCompound = np.array([0.6424, 0.2484, 0.06956, 0.0196])
DrotlistCompound = np.array([0.7234, 0.1869, 0.04041, 0.0341])
# 2 0.6423712078535424 0.7233534913946515
# 3 0.24844179777215183 0.1868690626702439
# 4 0.06956378582632722 0.04040898930201607
# 5 0.019619258022256694 0.034116943729134555

# Bound states definition, needed to calculate boundstate
boundStates = [1, 2, 3, 4]

# Chooses parent directory
parentDirectory = "../../data/pentamer/first_passage_times/"

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
    pickle_in = open("../../data/pentamer/MSMs/MSM_dimerForPentamer_t3.00E+06_s25_lagt" + str(lagtime)
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
    integrator = msmrdMultiParticleIntegrator(dt, seed, bodytype, numParticleTypes, radialBounds,
                                              unboundMSM, couplingMSM, DlistCompound, DrotlistCompound)
    integrator.setBoundary(boxBoundary)
    integrator.setDiscretization(discretization)

    # Generate random position and orientation particle list with two particles
    seed = int(trajectorynum)
    partlist = particleTools.randomParticleMSList(numParticles, boxsize,
                                                  minimumUnboundRadius, partTypes, [unboundMSM], seed)

    # Calculates the first passage times to a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    ii = 0
    while(unbound):
        ii += 1
        integrator.integrate(partlist)
        # Check status every 5000 timesteps (adds possible error of up to 5000*dt = 0.5), but
        # speeds up simulations
        if (ii % 5000 == 0):
            ringFormations = integrator.findClosedBindingLoops(partlist)
            if (3 in ringFormations):
                unbound = False
                return 'trimeric-loop', integrator.clock
            elif (4 in ringFormations):
                unbound = False
                return 'tetrameric-loop', integrator.clock
            elif (5 in ringFormations):
                unbound = False
                return "pentameric-loop", integrator.clock
            elif integrator.clock >= 400.0:
                unbound = False
                return 'Failed at:', integrator.clock
                #filenameLog = filename = "/run/media/maojrs/Mr300/Documents/Posdoc/projects/MSMRD2/" \
                #                         "msmrd2/data/pentamer/debug/eventLog_" + str(trajectorynum)
                #integrator.printEventLog(filenameLog)



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
            if state == 'trimeric-loop':
                print("Simulation " + str(index).zfill(5) + ", done. Trimeric loop formed!")
                if printAllRingFormations:
                    file.write(state + ' ' + str(time) + '\n')
            elif state == 'tetrameric-loop':
                print("Simulation " + str(index).zfill(5) + ", done. Tetrameric loop formed!!")
                if printAllRingFormations:
                    file.write(state + ' ' + str(time) + '\n')
            elif state == 'pentameric-loop':
                file.write(state + ' ' + str(time) + '\n')
                print("Simulation " + str(index).zfill(5) + ", done. PENTAMERIC loop formed!!!")
            else:
                print("Simulation " + str(index).zfill(5) + ", done. Failed :(")


# Run parallel code
multiprocessingHandler()

# # Serial code for testing with gdb
# with open(filename, 'w') as file:
#     for index in range(numTrajectories):
#         state, time = MSMRDsimulationFPT(index)
#         if state == 'pentamer':
#             file.write(state + ' ' + str(time) + '\n')
#             print("Simulation " + str(index) + ", done. Success!")
#             #input("Press the <ENTER> key to continue...")
#         elif state == 'loop':
#             print("Simulation " + str(index) + ", done. Failed by loop formation :(")
#             #input("Press the <ENTER> key to continue...")
#         else:
#             print("Simulation " + str(index) + ", done. Failed :(")
#             #input("Press the <ENTER> key to continue...")
