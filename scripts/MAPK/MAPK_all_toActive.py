import numpy as np
import msmrd2
from msmrd2.potentials import patchyProteinMAPK
from msmrd2.integrators import integratorMAPK 
import msmrd2.tools.particleTools as particleTools
import multiprocessing
from multiprocessing import Pool
import os

'''
Creates an MD simulation of an MAPK particle surrounded by a given number of kinases and
phosphotases and calculates their first passage times (FPTs) to activation from a random 
unbound configuration. The data is written to 
'../../data/MAPK/first_passage_times/filename_here.
'''

# Main parameters for particles and integrator
numMAPKs = 1 
numKinases = 3
numPhosphatases = 2
numparticles = numMAPKs + numKinases + numPhosphatases
D = 1.0
Drot = 1.0
particleTypes = [0, 1, 2]
# Main parameters for integrator
dt = 0.00001 #0.0001 #0.00001 #0.000005
bodytype = 'rigidbody'
anglePatches = np.pi/2
reactivationRateK = 1.0
reactivationRateP = 1.0
minimumUnboundRadius = 2.5
numTrajectories = 5000 #600 #10000

# Define Patchy Protein MAPK potential parameters (This values are fixed and should match
# those used to determine metastable states in potential and trajectory.)
sigma = 1.0
strength = 100 #65
angularStrength = 10 #2
patchesCoordinates1 = [np.array([np.cos(anglePatches/2), np.sin(anglePatches/2), 0.]), 
                       np.array([np.cos(-anglePatches/2), np.sin(-anglePatches/2), 0.])]
patchesCoordinates2 = [ np.array([np.cos(-anglePatches/2), np.sin(-anglePatches/2), 0.]) ]


# Define simulation boundaries parameters (choose either spherical or box)
boxsize = 5
boundaryType = 'periodic'

# Chooses parent directory
parentDirectory = "../../data/MAPK/first_passage_times/"

# Creates parent directory if it doesn't exist already
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("First passage times directory already exists. Simulation continues.")

# Create empty files to save the data in parallel algorithm
filename = parentDirectory  + 'MAPK_toActive_trajs' + str(numTrajectories) \
           +'_boxsize' + str(boxsize) + '.xyz'

def simulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state the active state
    :param trajectorynum: number of trajectories on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define simulation boundaries (choose either spherical or box)
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, 'periodic')

    # Define base seed
    seed = int(-1*trajectorynum) # Negative seed, uses random device as seed

    # Define potential
    potentialPatchyProteinMAPK = patchyProteinMAPK(sigma, strength, patchesCoordinates1, patchesCoordinates2)

    # Creates random particle list
    seed = int(trajectorynum)
    partlist, mapkIndex, kinaseIndex, phosIndex = particleTools.randomMAPKparticleList(numMAPKs, numKinases,
                                                                                       numPhosphatases, boxsize,
                                                                                       minimumUnboundRadius, D,
                                                                                       Drot, seed)
	
    # Define integrator and boundary (over-damped Langevin)
    integrator = integratorMAPK(dt, seed, bodytype, anglePatches, reactivationRateK, reactivationRateP,
                                 mapkIndex, kinaseIndex, phosIndex)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyProteinMAPK)

    # Calculates the first passage times for all MAPK particles to activate. Each trajectory is integrated until
    # this event. The output in the files is the elapsed time.
    inactive = [True] * len(mapkIndex)
    while(inactive != [False] * len(mapkIndex)):
        integrator.integrate(partlist)
        for i, index in enumerate(mapkIndex):
            if partlist[index].state == 2:
                inactive[i] = False
            else:
                inactive[i] = True
        if inactive == [False] * len(mapkIndex):
            return 'active', integrator.clock
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
            if state == "active":
                file.write(state + ' ' + str(time) + '\n')
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
#         if state == "active":
#             file.write(str(state) + ' ' + str(time) + '\n')
#             print("Simulation " + str(index) + ", done. Success!")
#         else:
#             print("Simulation " + str(index) + ", done. Failed :(")


