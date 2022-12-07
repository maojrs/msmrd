import numpy as np
import msmrd2
from msmrd2.integrators import integratorMoriZwanzig
from msmrd2.potentials import combinedPairPotential, WCA, pairBistable
import msmrd2.tools.particleTools as particleTools
import msmrd2.tools.analysis as analysisTools
import msmrd2.tools.trajectoryTools as trajectoryTools


import multiprocessing
from multiprocessing import Pool
from functools import partial
import random
import sys
import os

## Units

# ### Boltzman constant
# - $k_B = 1.38064852 \times 10^{-23} \frac{m^2}{s^2} \frac{kg}{K} \left(= \frac{nm^2}{ns^2}\frac{kg}{K}\right)$.
#
# ### Basic units
# - Length (l): nanometers (nm)
# - Energy ($\epsilon$) : $k_B T = g/mol \frac{nm^2}{ns^2}$
# - Mass (m): gram/mol (g/mol)
#
# ### Derived units
# - Time: $l\sqrt{m/\epsilon} = ns$
# - Temperature: $\epsilon/k_B =$ Kelvin ($K$)
# - Force: $\epsilon/l = kg \frac{nm}{ns^2}$
#
# ### Reduced quantities (dimensionless)
# - Reduced pair potential: $U^* = U/\epsilon$
# - Reduced force: $F^* = F l /\epsilon$
# - Reduced distance: $r^* = r/l$
# - Reduced density: $\rho^*=\rho l^3$
# - Reduced Temperature: $T^* = k_B T/\epsilon$
# - Reduced Pressure: $P^* = Pl^3/\epsilon$
# - Reduced friction: $\sigma^2/time$

# Main parameters
numBathParticles = 500 #500 #500
numparticles = 2 + numBathParticles #Added distinguished particle (index 0)
# D = 3.0E-2 #1.0E-3 #(nm^2/ns) Note 1.0E-3 nm^2/ns = 1 micrometer^2/s #0.1
Gamma = 0.3 #30 # Friction coefficient (units of KbT/D = mass over time (gram/mol)/ns)
particlemass = 18.0 # (g/mol) approximately mass of water
distinguishedParticleMass = 3 * particlemass # (kg)
particleDiameter = 0.5 # (nm)
separationDistance = 1 * particleDiameter # minimum separation distance for initial condition
numSimulations = 1000 #250 #500
# For computations, we assume KbT=1, thus the force F must be: F=KbT f, where f is the force computed
# from the potential. This means the plotted potential is on reduced units (not the distances though);
KbT = 1

# Main parameters for integrator
dt = 0.05 #0.005 #0.001 #0.0005 # (ns)
seed = -1 # Seed = -1 used random device as seed
bodytype = 'point'

# Define simulation boundaries (choose either spherical or box)
boxsize = 8 #(nm)
boundaryType = 'periodic'

# Parameters for WCA potential (rm=2^(1/6)sigma)
epsilon = 1.0
rm = particleDiameter # (diameter in nm)
sigma = rm * 2**(-1/6)
excludeParticleTypesPairs = [1] # Two particles of this same type wll not interact with each other.


# Parameters for pair bistable potential (will only act on distinguished particles (type 1))
x0 = 1.0*particleDiameter # location of first minima
rad = 1.0*particleDiameter # half the distance between minimas
scalefactor = 2

# Parameters for FPT calculations
transitionType = 'CO' # CO or OC, closed to open or open to closed
if transitionType == 'CO':
    initialSeparation = 1*x0 # Either first minima: x0 or second minima: 2*rad
    finalSeparation = 2*rad # Either first minima: x0 or second minima: 2*rad
else:
    initialSeparation = 2*rad # Either first minima: x0 or second minima: 2*rad
    finalSeparation = 1*x0 # Either first minima: x0 or second minima: 2*rad
minimaThreshold = 1.9*rad

# Simulation parameters
tfinal = 10000 #100 #10000 #100000
timesteps = int(tfinal/dt)
bufferSize = 100 * 1024
stride = 1 #50 #5
outTxt = False
outH5 = True
outChunked = True
trajtype = "moriZwanzigVelocity" # Samples position and velocity of distinguished particle (type 1) + raux variables
distinguishedTypes = [1]

# Parent directory location
#parentDirectory = "../../data/MoriZwanzig/bistable/"
parentDirectory = os.environ['DATA'] + 'stochasticClosure/dimer/boxsize' + str(boxsize) + '/'

# Create folder for data
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("Folder " + parentDirectory + " already exists.")
    proceed = True

# Create folder for benchmark data
foldername = "benchmarkFPTcomparison"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder stochasticClosure/dimer/boxsize" + str(boxsize) + "/" + foldername + " already exists. Previous data files might be overwritten. Continue, y/n?")
    proceed = input()
    if proceed != 'y':
        sys.exit()

# Create parameter dictionary to write to parameters reference file
parameterfilename = os.path.join(filedirectory, "parameters")
parameterDictionary = {'numFiles' : numSimulations, 'numParticles' : numparticles, 'dt' : dt, 'bodytype' : bodytype,
                       'Gamma' : Gamma, 'sigma' : sigma, 'KbT' : KbT, 'mass' : distinguishedParticleMass,
                       'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype,
                       'boxsize' : boxsize, 'boundaryType' : boundaryType}
analysisTools.writeParameters(parameterfilename, parameterDictionary)


# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simMoriZwanzig")

# Create empty files to save the data in parallel algorithm
filename = filedirectory  + '/simMoriZwanzigFPTs_' + transitionType + '_box' + str(boxsize) + '_nsims' + str(numSimulations) +'.xyz'

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list
    seed = int(simnumber)
    random.seed(seed)
    partlist = particleTools.randomLangevinParticleListTwoDistinguished(numparticles, boxsize, separationDistance,
                                                                        particlemass, seed, initialSeparation)
    # Set distinguished particle (default type is zero)
    partlist[0].setType(1)
    partlist[1].setType(1)
    partlist[0].setMass(distinguishedParticleMass)
    partlist[1].setMass(distinguishedParticleMass)

    # Define boundary
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)

    # Define potentials
    potentialWCA = WCA(epsilon, sigma, excludeParticleTypesPairs)
    potentialWCA.setForceCapValue(100.0)
    potentialPairBistable = pairBistable(x0, rad, distinguishedTypes, scalefactor)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
    integrator = integratorMoriZwanzig(dt, seed, bodytype, Gamma)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPairBistable)
    integrator.setAuxPairPotential(potentialWCA)  # Aux so it can be saved into aux variable
    integrator.setDistinguishedTypes(distinguishedTypes)
    integrator.setKbT(KbT)

    # WORKING HEEEREEEEEEEE WIPPPP
    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    while(unbound):
        integrator.integrate(partlist)
        relativePosition = trajectoryTools.relativePosition(partlist[0].position, partlist[1].position,
                                                            boundaryType, boxsize)
        relDistance = np.linalg.norm(relativePosition)
        if np.abs(relDistance - initialSeparation) >= minimaThreshold:
            unbound = False
            return 'success', integrator.clock
        elif integrator.clock >= 10000.0:
            unbound = False
            return 'Failed at:', integrator.clock
    print("Simulation " + str(simnumber) + ", done.")

def multiprocessingHandler():
    '''
    Handles parallel processing of simulationFPT and writes to same file in parallel
    '''
    # Runs several simulations in parallel
    num_cores = multiprocessing.cpu_count() - 1
    pool = Pool(processes=num_cores)
    trajNumList = list(range(numSimulations))
    with open(filename, 'w') as file:
        for index, result in enumerate(pool.imap(runParallelSims, trajNumList)):
            status, time = result
            if status == 'success':
                file.write(str(time) + '\n')
                print("Simulation " + str(index) + ", done. Success!")
            else:
                print("Simulation " + str(index) + ", done. Failed :(")

# Run parallel code
multiprocessingHandler()

## Run serial for debugging
#for i in range(numSimulations):
#    runParallelSims(i)
