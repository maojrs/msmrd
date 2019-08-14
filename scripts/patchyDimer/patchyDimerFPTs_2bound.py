import numpy as np
import msmrd2
from msmrd2.potentials import patchyParticleAngular
from msmrd2.integrators import overdampedLangevin as odLangevin
import msmrd2.tools.particleTools as particleTools
import multiprocessing
from multiprocessing import Pool

'''
Creates an MD simulation of two particle and calculates their first passage times (FPTs) from a random 
unbound configuration to either of the two bound states (A or B). The data is written to 
'../data/dimer/first_passage_times/filename_here.
'''

# Main parameters for particle and integrator
numparticles = 2
dt = 0.0001 #0.00001 #0.000005
bodytype = 'rigidbody'
D = 1.0
Drot = 1.0
relativeDistanceCutOff = 2.2
numTrajectories = 10000
# Other important parameters
boxsize = 6

# Parameters of patchy Particle potential with angular dependence (make sure it is consistent with msmrd data)
sigma = 1.0
strength = 100.0
angularStrength = 10.0
angleDiff = 3*np.pi/5.0
patch1 = np.array([np.cos(angleDiff/2),np.sin(angleDiff/2),0.])
patch2 = np.array([np.cos(-angleDiff/2),np.sin(-angleDiff/2),0.])
patchesCoordinates = [patch1, patch2]

# Bound states, needed to calculate boundstate of patchydimer
boundStatesA = [1, 2, 5, 6] # U-shaped bound dimer, corresponds to A state
boundStatesB = [3, 4, 7, 8] # Zigzag-shaped bound dimer, corresponds to B state

# Create empty files to save the data in parallel algorithm
filename = '../../data/dimer/first_passage_times/patchyDimerFPTs_2bound_trajs' + str(numTrajectories) \
           +'_boxsize' + str(boxsize) + '.xyz'

def simulationFPT(trajectorynum):
    '''
    Calculates first passage time of a trajectory with random initial
    conditions and final state a bound state
    :param trajectorynum: number of trajectories on which to calculate the FPT
    :return: state, first passage time
    '''

    # Define dummy trajectory to extract bound states from python (needed to use getState function)
    dummyTraj = msmrd2.trajectories.patchyDimer(2,1)

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
    partlist = particleTools.randomParticleList(numparticles, boxsize, relativeDistanceCutOff, D, Drot, seed)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    while(unbound):
        integrator.integrate(partlist)
        boundState = dummyTraj.getState(partlist[0], partlist[1])
        if boundState in boundStatesA:
            unbound = False
            return 'A', integrator.clock
        elif boundState in boundStatesB:
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
        for index, result in enumerate(pool.imap(simulationFPT, trajNumList)):
            state, time = result
            if state == 'A' or state == 'B':
                file.write(state + ' ' + str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
            else:
                print("Simulation " + str(index) + ", done. Failed :(")


# Run parallel code
multiprocessingHandler()




