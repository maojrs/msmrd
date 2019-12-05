import pickle
import numpy as np
import pyemma
import msmrd2.tools.trajectoryTools as trajectoryTools
import msmrd2.tools.analysis as analysisTools

# Given discrete data, creates MSMs for a list of lagtimes. To see
# implied time scales use notebook version of this script.

# Load parameters from parameters file (from original MD simulation)
parentDirectory = '../../data/patchyProtein/benchmark/'
parameterDictionary = analysisTools.readParameters(parentDirectory + "parameters")
nfiles = parameterDictionary['numFiles']
dt = parameterDictionary['dt']
stride = parameterDictionary['stride']
totalTimeSteps = parameterDictionary['timesteps']

# Calculated parameters
dtEffective = dt*stride # needed to obtain rate dictionary
fnamebase = parentDirectory + 'simPatchyProtein_'
#fnamebase = parentDirectory + 'fullSphereQuatOffset/simPatchyProtein_'

# Parameters for MSM generation
numBoundStates = 6
lagtimes = [100, 150, 175, 200, 250] #[50, 75, 100, 150] # [75]
reversible = True #False
stitching = True


# Load discrete trajectories
dtrajs = []
fnamesuffix = '_discrete' #'_discrete_test' #'_discrete_python' # '_discrete' # '_discrete_new'
filetype = 'xyz' # 'h5' or 'xyz'
for i in range(nfiles):
    dtraj = trajectoryTools.loadDiscreteTrajectory(fnamebase, i, fnamesuffix, filetype)
    dtrajs.append(dtraj)
    print("Loading discrete trajectory ", i+1, " of ", nfiles, " done.", end="\r")
print("Loaded all trajectories!                      ")


# Slice trajectories getting rid of the unbound state 0
unboundStateIndex = 0
slicedDtrajs = trajectoryTools.splitDiscreteTrajs(dtrajs, unboundStateIndex)
# Stitch trajectories if wanted
print("All trajectories sliced, now stitching trajectories ...")
if stitching:
    minLength = 5000
    finalTrajs = trajectoryTools.stitchTrajs(slicedDtrajs, minLength)
else:
    finalTrajs = slicedDtrajs
print("Stitched all trajectories!")


# Loop over all desired lagtime to obatin MSMs
for i, lagtime in enumerate(lagtimes):
    print("Generating MSM ", i, "of ", len(lagtimes),  " at lagtime ", lagtime, end="\r")

    # Create MSM between transision states and bound states without stitching
    mainmsm = pyemma.msm.estimate_markov_model(finalTrajs, lagtime, reversible=reversible)
    # The active set keep track of the indexes used by pyemma and the ones used to describe the state in our model.
    activeSet = mainmsm.active_set

    # Pickle MSM transition matrix and active set as a dictionary
    MSM = {'transition_matrix' : mainmsm.transition_matrix, 'active_set': mainmsm.active_set}
    pickle_out = open(parentDirectory + "MSM_patchyProtein_t" + "{:.2E}".format(totalTimeSteps ) +
                      "_s" + "{:d}".format(stride) + "_lagt" + "{:d}".format(lagtime) + ".pickle","wb")
    pickle.dump(MSM, pickle_out)
    pickle_out.close()
print("Generated all MSMs.")
