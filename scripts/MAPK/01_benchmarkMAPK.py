import numpy as np
import msmrd2
from msmrd2.potentials import patchyProteinMAPK
from msmrd2.integrators import integratorMAPK
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
# - Energy ($\epsilon$) : $k_B T = kg \frac{nm^2}{ns^2}$
# - Mass (m): kilogram (kg)
#
# ### Derived units
# - Time: $l\sqrt{m/\epsilon} = nm$
# - Temperature: $\epsilon/k_B =$ Kelvin ($K$)
# - Force: $\epsilon/l = kg \frac{nm}{ns^2}$
#
# ### Reduced quantities (dimensionless)
# - Reduced pair potential: $U^* = U/\epsilon$
# - Reduced distance: $r^* = r/l$
# - Reduced density: $\rho^*=\rho l^3$
# - Reduced Temperature: $T^* = k_B T/\epsilon$
# - Reduced Pressure: $P^* = Pl^3/\epsilon$

# Main parameters
numMAPKs = 1
numKinases = 1
numPhosphatases = 0
numparticles = numMAPKs + numKinases + numPhosphatases
sigma = 5.0 #(nanometers) particles diameters
boxsize = 15 #(nanometers) #2.3
D = 1.0E-3 #(nm^2/ns) Note 1.0E-3 nm^2/ns = 1 micrometer^2/s
Drot = 1.6E-4 #(rad^2/ns) Note 1.6E-4 rad^2/ns = 1.6E5 rad^2/s
kB = 1.38064852E-23 #(Boltzmann constant (nm^2/ns^2 * kg/K))
Temp = 300 #(Kelvin)
particleTypes = [0, 1, 2]

# Main parameters for integrator
dt = 0.1 #(nanoseconds)
bodytype = 'rigidbody'
anglePatches = np.pi/2
trel = 1000 # (nanoseconds) Note 1000 ns = 1 micro second. Can be varied up to 10 miliseconds
reactivationRateK = np.log(2)/trel # not relevant for MSMRD parametrization
reactivationRateP = np.log(2)/trel # not relevant for MSMRD parametrization
minimumUnboundRadius = 1.25 * sigma
numSimulations = 4 #500

# Simulation parameters
timesteps = 100000000 #(0.01 second) #3000000 #3000000
bufferSize = 1024
stride = 500
outTxt = False
outH5 = True
outChunked = True
trajtype = "MAPK" # "trajectoryPositionOrientationState"

# Define Patchy Protein MAPK potential parameters (This values are fixed and should match
# those used to determine metastable states in potential and trajectory.)
strength = 100 #100 #65
patchesCoordinates1 = [np.array([np.cos(anglePatches/2), np.sin(anglePatches/2), 0.]),
                       np.array([np.cos(-anglePatches/2), np.sin(-anglePatches/2), 0.])]
patchesCoordinates2 = [ np.array([np.cos(-anglePatches/2), np.sin(-anglePatches/2), 0.]) ]
potentialPatchyProteinMAPK = patchyProteinMAPK(sigma, strength, patchesCoordinates1, patchesCoordinates2)

# Define simulation boundaries (choose either spherical or box)
boundaryType = 'periodic'
boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)

# Parent directory location
parentDirectory = "../../data/MAPK/"

# Create folder for data
try:
    os.mkdir(parentDirectory)
except OSError as error:
    print("Folder MAPK already exists.")
    proceed = True

# Create folder for benchmark data        
foldername = "benchmark"
if numMAPKs == 1 and numKinases == 1 and numPhosphatases == 0:
    foldername = "benchmark_kinase"
elif numMAPKs == 1 and numKinases == 0 and numPhosphatases == 1:
    foldername = "benchmark_phosphatase"
filedirectory =  os.path.join(parentDirectory, foldername)
try:
    os.mkdir(filedirectory)
except OSError as error:
    print("Folder MAPK/benchmark already exists. Previous data files might be overwritten. Continue, y/n?")
    proceed = input()
    if proceed != 'y':
        sys.exit()

# Create parameter dictionary to write to parameters reference file
parameterfilename = os.path.join(filedirectory, "parameters")
parameterDictionary = {'numFiles' : numSimulations, 'numParticles' : numparticles, 'dt' : dt, 'bodytype' : bodytype,
                       'D' : D, 'Drot' : Drot, 'timesteps' : timesteps, 'stride' : stride, 'trajtype' : trajtype,
                       'boxsize' : boxsize, 'boundaryType' : boundaryType, 'potentialStrength' : strength,
                       'potentialAngularStrength' : angularStrength, 'anglePatches' : anglePatches}
analysisTools.writeParameters(parameterfilename, parameterDictionary)


# Provides base filename (folder must exist (and preferably empty), otherwise H5 might fail)
basefilename = os.path.join(filedirectory, "simMAPK")

# Simulation wrapper for parallel runs
def runParallelSims(simnumber):

    # Define particle list
    seed = int(simnumber)
    random.seed(seed)
    defaultOrientvector = patchesCoordinates2[0]
    partlist, mapkIndex, kinaseIndex, phosIndex = particleTools.randomMAPKparticleList(numMAPKs, numKinases,
            numPhosphatases, boxsize, minimumUnboundRadius, D, Drot, seed, defaultOrientvector)

    # Integrator definition
    seed = int(-1*simnumber) # random seed (negative and different for every simulation, good for parallelization)
    integrator = integratorMAPK(dt, seed, bodytype, anglePatches, reactivationRateK, reactivationRateP,
                                sigma, mapkIndex, kinaseIndex, phosIndex)
    integrator.setBoundary(boxBoundary)
    integrator.setPairPotential(potentialPatchyProteinMAPK)
    integrator.setKbT(kB * Temp)
    integrator.disableDeactivation()

    # Creates simulation
    sim = msmrd2.simulation(integrator)

    # Output filename definitionconj
    filename = basefilename + "_{:04d}".format(simnumber)

    # Runs simulation
    sim.run(partlist, timesteps, stride, bufferSize, filename, outTxt, outH5, outChunked, trajtype)
    print("Simulation " + str(simnumber) + ", done.")

# Runs several simulations in parallel
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)
