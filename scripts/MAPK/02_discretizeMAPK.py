import msmrd2
import msmrd2.tools.analysis as analysisTools
import numpy as np

# Given an MD simulation of patchyProteins, generate discrete trajectories, using patchyProtin discrete trajectory.

benchmark_type = '_kinase' # '_phosphatase' # ''

# Load parameters from parameters file (from original MD simulation)
parentDirectory = '../../data/MAPK/benchmark' + benchmark_type + '/'
parameterDictionary = analysisTools.readParameters(parentDirectory + "parameters")
nfiles = 2 #parameterDictionary['numFiles']
numParticles = parameterDictionary['numParticles']
dt = parameterDictionary['dt']
sigma = parameterDictionary['sigma']
stride = parameterDictionary['stride']
totalTimeSteps = parameterDictionary['timesteps']
boxsize = parameterDictionary['boxsize']
boundaryType = parameterDictionary['boundaryType']
anglePatches = parameterDictionary['anglePatches']

# Calculated parameters
effectivetimeSteps = int(totalTimeSteps/stride)
fnamebase = parentDirectory + 'simMAPK_'
bufferSize = effectivetimeSteps

# Set trajectory discretizator
radialLowerBound = sigma*1.25 #default values
radialUpperBound = sigma*2.25 #default values
numSphericalSectionsPos = 6 #default values
numSphericalSectionsOrientvec = 6 #default values
discretizator = msmrd2.trajectories.MAPK(numParticles, bufferSize, anglePatches,
                                         numSphericalSectionsPos, numSphericalSectionsOrientvec,
                                         radialLowerBound, radialUpperBound)
discretizator.setTolerances(0.6,0.6)

# Set boundary (important for discretizer)
boxBoundary = msmrd2.box(boxsize, boxsize, boxsize, boundaryType)
discretizator.setBoundary(boxBoundary)

# Loads H5 files and generates discrete trajectories at once directly on c++.
dtrajs = []
for i in range(nfiles):
    filename = fnamebase + str(i).zfill(4) + '.h5'
    dtraj = discretizator.discretizeTrajectoryH5(filename)
    dtrajs.append(dtraj)
    print("Loading file ", i+1, " of ", nfiles, " done.", end="\r")
print("\nDone loading files")

# Write discrete trajectory to xyz files
for i, dtraj in enumerate(dtrajs):
    datafile  = open(fnamebase + str(i).zfill(4) + '_discrete.xyz', 'w')
    for j in range(len(dtraj)):
        datafile.write(str(dtraj[j]) + '\n')
    datafile.close()
    print("Writing discrete trajectory ", i+1, " of ", nfiles, " done.", end="\r")
print("\nDone writing discrete trajectories")
