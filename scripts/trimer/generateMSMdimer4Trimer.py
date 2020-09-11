import pickle
import numpy as np
import pyemma
import pyemma.plots as mplt
import matplotlib.pyplot as plt
import msmrd2.tools.trajectoryTools as trajectoryTools
import msmrd2.tools.analysis as analysisTools

# Given discrete data, creates MSMs for a list of lagtimes. To see
# implied time scales use notebook version of this script.

# Plot implied time scales
plotImpliedTimescalesDraft = True
plotImpliedTimescalesPaperVersion = False

# Load parameters from parameters file (from original MD simulation)
parentDirectory = '../../data/trimer/benchmark/'
#parentDirectory = '/group/ag_cmb/scratch/maojrs/msmrd2_data/trimer/benchmark/'
parameterDictionary = analysisTools.readParameters(parentDirectory + "parameters")
nfiles = parameterDictionary['numFiles']
dt = parameterDictionary['dt']
stride = parameterDictionary['stride']
totalTimeSteps = parameterDictionary['timesteps']

# Calculated parameters
dtEffective = dt*stride # needed to obtain rate dictionary
fnamebase = parentDirectory + 'simDimer4Trimer_'

# Parameters for MSM generation
numBoundStates = 4
lagtimes = [100, 150, 175, 200, 250, 300] #[50, 75, 100, 150] # [75]
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
print("\nLoaded all trajectories!")


# Slice trajectories getting rid of the unbound state 0
unboundStateIndex = 0
slicedDtrajs = trajectoryTools.splitDiscreteTrajs(dtrajs, unboundStateIndex)
# Stitch trajectories if wanted
print("\nSliced all trajectories!")
if stitching:
    minLength = 5000
    finalTrajs = trajectoryTools.stitchTrajs(slicedDtrajs, minLength)
else:
    finalTrajs = slicedDtrajs
print("\nStitched all trajectories!")


# Loop over all desired lagtime to obatin MSMs
for i, lagtime in enumerate(lagtimes):
    print("Generating MSM ", i+1, "of ", len(lagtimes),  " at lagtime ", lagtime, end="\r")

    # Create MSM between transision states and bound states without stitching
    mainmsm = pyemma.msm.estimate_markov_model(finalTrajs, lagtime, reversible=reversible)
    # The active set keep track of the indexes used by pyemma and the ones used to describe the state in our model.
    activeSet = mainmsm.active_set

    # Pickle MSM transition matrix and active set as a dictionary
    MSM = {'transition_matrix' : mainmsm.transition_matrix, 'active_set': mainmsm.active_set}
    pickle_out = open(parentDirectory + "MSM_dimer4trimer_t" + "{:.2E}".format(totalTimeSteps ) +
                      "_s" + "{:d}".format(stride) + "_lagt" + "{:d}".format(lagtime) + ".pickle","wb")
    pickle.dump(MSM, pickle_out)
    pickle_out.close()
print("\nGenerated all MSMs.")

# Generate implied timescales plots (draft verision)
if plotImpliedTimescalesDraft:
    print("Generating implied timescales plots (draft version)")
    maxlagtime = 300 #100 #200 #300
    its = pyemma.msm.its(finalTrajs, maxlagtime, reversible=reversible)
    nits = 20
    fig, ax = plt.subplots(figsize=(10, 7))
    mplt.plot_implied_timescales(its, ax = ax, nits = nits, ylog=True, units='steps', linewidth=2, dt=1)
    plt.ylabel(r"log(timescale/steps)", fontsize = 18)
    plt.xlabel(r"lag time/steps", fontsize = 18)
    plt.savefig(parentDirectory + 'its_draft_dimer4trimer' + "{:.2E}".format(totalTimeSteps) + '.pdf')
    # No log version
    fig, ax = plt.subplots(figsize=(10, 7))
    mplt.plot_implied_timescales(its, ax = ax, nits = nits, ylog=False, units='steps', linewidth=2, dt=1)
    plt.ylabel(r"timescale/steps", fontsize = 24)
    plt.xlabel(r"lag time/steps", fontsize = 24)
    plt.savefig(parentDirectory + 'its_draft_nolog_dimer4trimer' + "{:.2E}".format(totalTimeSteps) + '.pdf')
    print("Finished draft plots.")

# Generate implied timescales plots (paper verision w/error bars)	
if plotImpliedTimescalesPaperVersion:
    print("Generating implied timescales plots (paper version)")
    maxlagtime = 300
    its = pyemma.msm.its(finalTrajs, maxlagtime, reversible=reversible, errors='bayes')
    nits = 20
    fig, ax = plt.subplots(figsize=(10, 7))
    mplt.plot_implied_timescales(its, ax = ax, nits = nits, ylog=True, units='steps', linewidth=2, dt=1, 
    	                         show_mean=False, markersize=0, confidence=0.95)
    plt.ylabel(r"log(timescale/steps)", fontsize = 24)
    plt.xlabel(r"lag time/steps", fontsize = 24)
    plt.xticks(fontsize = 28)
    plt.yticks(fontsize = 28)
    #plt.ylim([10.0,100000])
    plt.savefig(parentDirectory + 'its_paper_dimer4trimer' + "{:.2E}".format(totalTimeSteps) + '.pdf')
    # No log version
    fig, ax = plt.subplots(figsize=(10, 7))
    mplt.plot_implied_timescales(its, ax = ax, nits = nits, ylog=False, units='steps', linewidth=2, dt=1, 
    	                         show_mean=False, markersize=0, confidence=0.95)
    plt.ylabel(r"timescale/steps", fontsize = 24)
    plt.xlabel(r"lag time/steps", fontsize = 24)
    plt.xticks(fontsize = 28)
    plt.yticks(fontsize = 28)
    #plt.ylim([10.0,4500])
    plt.savefig(parentDirectory + 'its_paper_nolog_dimer4trimer' + "{:.2E}".format(totalTimeSteps) + '.pdf')
    print("Finished paper plots.")
