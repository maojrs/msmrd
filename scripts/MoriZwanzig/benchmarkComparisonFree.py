import numpy as np
import msmrd2
from msmrd2.integrators import integratorMoriZwanzig
from msmrd2.potentials import WCA
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
particleDiameter = 0.5 # (nm)
separationDistance = 1 * particleDiameter # minimum separation distance for initial condition
numSimulations = 100 #250 #500
# For computations, we assume KbT=1, thus the force F must be: F=KbT f, where f is the force computed
# from the potential. This means the plotted potential is on reduced units (not the distances though);
KbT = 1

# Main parameters for integrator
dt = 0.05 #0.005 #0.001 #0.0005 # (ns)
seed = -1 # Seed = -1 used random device as seed
bodytype = 'point'

# Define simulation boundaries (choose either spherical or box)
boxsize = 5 #(nm)
boundaryType = 'periodic'

# Parameters for WCA potential (rm=2^(1/6)sigma)
epsilon = 1.0
rm = particleDiameter # (diameter in nm)
sigma = rm * 2**(-1/6)

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
equilibrationSteps = 10000


# Parent directory location
#parentDirectory = "../../data/MoriZwanzig/free/"
parentDirectory = os.environ['DATA'] + 'stochasticClosure/free/boxsize' + str(boxsize) + '/'


# Create folder for data
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("Folder " + parentDirectory + " already exists.")
    proceed = True

# Create folder for benchmark data
foldername = "benchmarkComparison"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder " + foldername + " already exists. Previous data files might be overwritten. Continue, y/n?")
    proceed = input()
    if proceed != 'y':
        sys.exit()

# Create parameter dictionary to write to parameters reference file
parameterfilename = os.path.join(filedirectory, "parameters")
parameterDictionary = {'numFiles' : numSimulations, 'numParticles' : numparticles, 'dt' : dt, 'bodytype' : bodytype,
                       'Gamma' : Gamma, 'sigma' : sigma, 'KbT' : KbT, 'mass' : distinguishedParticleMass,
                       'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype,
                       'boxsize' : boxsize, 'boundaryType' : boundaryType, 'equilibrationSteps' : equilibrationSteps}
analysisTools.writeParameters(parameterfilename, parameterDictionary)


# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simMoriZwanzig")

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list
    seed = int(simnumber)
    random.seed(seed)
    partlist = particleTools.randomLangevinParticleList(numparticles, boxsize, separationDistance,
                                                        particlemass, seed, distinguishedParticleOrigin=True)

    print('Particle collocation of simulation ' + str(simnumber) + ' done.')

    # Set distinguished particle (default type is zero)
    partlist[0].setType(1)
    partlist[0].setMass(distinguishedParticleMass)

    # Define boundary
    boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)

    # Define pair potential
    potentialWCA = WCA(epsilon, sigma)
    potentialWCA.setForceCapValue(100.0)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
    integrator = integratorMoriZwanzig(dt, seed, bodytype, Gamma)
    integrator.setBoundary(boxBoundary)
    integrator.setAuxPairPotential(potentialWCA) # Aux so it can be saved into aux variable
    integrator.setDistinguishedTypes(distinguishedTypes)
    integrator.setKbT(KbT)

    # Creates simulation
    sim = msmrd2.simulation(integrator)

    # Output filename definitionconj
    filename = basefilename + "_{:04d}".format(simnumber)

    # Runs simulation
    sim.setEquilibrationSteps(equilibrationSteps)
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked, trajtype)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
num_cores = multiprocessing.cpu_count() - 1
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)

## Run serial for debugging
#for i in range(numSimulations):
#    runParallelSims(i)