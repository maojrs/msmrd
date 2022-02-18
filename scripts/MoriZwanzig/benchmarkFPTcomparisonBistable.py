import numpy as np
import msmrd2
from msmrd2.integrators import integratorMoriZwanzig
from msmrd2.potentials import WCA, bistable2
import msmrd2.tools.particleTools as particleTools
import msmrd2.tools.analysis as analysisTools
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
numparticles = 1 + numBathParticles #Added distinguished particle (index 0)
# D = 3.0E-2 #1.0E-3 #(nm^2/ns) Note 1.0E-3 nm^2/ns = 1 micrometer^2/s #0.1
Gamma = 0.3 #30 # Friction coefficient (units of KbT/D = mass over time (gram/mol)/ns)
particlemass = 18.0 # (g/mol) approximately mass of water
distinguishedParticleMass = 3 * particlemass # (kg)
particleDiameter = 0.3 # (nm)
separationDistance = 2 * particleDiameter # minimum separation distance for initial condition
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

# Parameters for external potential (will only acts on distinguished particles (type 1))
# potential of the form ax^4 - bx^2 + c y^2 + d z^2
xminima = 2.5
a = 0.025
b = 2 * a * xminima**2
c = 1 * b
d = 1 * b
potentialParams = [a, b, c, d]
scalefactor = 1

# Parameters for FPT calculations
initialPosition = np.array([-1.0*minimaDist, 0., 0.])
finalPosition = np.array([1.0*minimaDist, 0., 0.])
minimaThreshold = 0.3

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
parentDirectory = os.environ['DATA'] + 'stochasticClosure/bistable/'

# Create folder for data
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("Folder stochasticClosure/bistable/ already exists.")
    proceed = True

# Create folder for benchmark data
foldername = "benchmarkFPTcomparison"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder stochasticClosure/bistable/" + foldername + " already exists. Previous data files might be overwritten. Continue, y/n?")
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
filename = filedirectory  + '/simMoriZwanzigFPTs_' + str(numSimulations) + '.xyz'

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list
    seed = int(simnumber)
    random.seed(seed)
    partlist = particleTools.randomLangevinParticleList(numparticles, boxsize, separationDistance,
                                                        particlemass, seed, distinguishedParticleOrigin=True)
    # Set distinguished particle (default type is zero)
    partlist[0].setType(1)
    partlist[0].setMass(distinguishedParticleMass)

    # Define boundary
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)

    # Define pair potential
    potentialWCA = WCA(epsilon, sigma)
    potentialWCA.setForceCapValue(100.0)

    # Define external potential
    externalPotential = bistable2(potentialParams, distinguishedTypes, scalefactor)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
    integrator = integratorMoriZwanzig(dt, seed, bodytype, Gamma)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialWCA)
    integrator.setExternalPotential(externalPotential)
    integrator.setDistinguishedTypes(distinguishedTypes)
    integrator.setKbT(KbT)

    # Calculates the first passage times for a given bound state. Each trajectory is integrated until
    # a bound state is reached. The output in the files is the elapsed time.
    unbound = True
    while(unbound):
        integrator.integrate(partlist)
        distanceFinalState = np.linalg.norm(partlist[0].position - finalPosition)
        if distanceFinalState <= minimaThreshold:
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