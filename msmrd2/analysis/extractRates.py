import h5py
import numpy as np
import msmrd2.tools as msmrdtls
import itertools
import pyemma
import msmtools.analysis


# Functions to help extract rates from data and/or MSMs


def loadTrajectory(fnamebase, fnumber):
    '''
    Reads data from discrete trajectory and returns a simple np.array of
    integers representing the discrete trajectory
    :param fnamebase, base of the filename
    :param fnumber, filenumber
    :return: array of arrays representing the trajectory
    '''
    filename = fnamebase + str(fnumber).zfill(4) + '.h5'
    f = h5py.File(filename, 'r')

    # Get the data
    a_group_key = list(f.keys())[0]
    data = np.array(f[a_group_key])

    return data


def loadDiscreteTrajectory(fnamebase, fnumber):
    '''
    Reads data from discrete trajectory and returns a simple np.array of
    integers representing the discrete trajectory
    :param fnamebase, base of the filename
    :param fnumber, filenumber
    :return: array with integers representing the discrete trajectory
    '''
    filename = fnamebase + str(fnumber).zfill(4) + '_discrete.h5'
    f = h5py.File(filename, 'r')

    # Get the data
    a_group_key = list(f.keys())[0]
    data = np.array(f[a_group_key])

    # Transform data into simple 1D array
    array = np.zeros(len(data), dtype = int)
    for i in range(len(data)):
        array[i] = data[i][0]
    return array


def createStatesDictionaries(boundstates, orientations):
    '''
    Given number of bound states and number of orientation transition states,
    creates two dictionaries to count events and accumulated time for all relevant transitions.
    This is required to obtain the rates from the discrete trajectories.
    :param boundstates: number of bound states. They will be labelled as b0, b1, ....
    :param orientations: number of discrete orientations of one particle, i.e the number of partitions in the
    sphere around one particle, which determine the number of possible transition states between
    relative orientations (e.g. for 3 => 11, 12, 13, 22, 23, 33).
    :return: {timecountDict, eventcountDict} empty dictionaries with defined keys that map state label
    to accumulated time to transition and to number of events found for that particular transition.
    The dictionary keys have the form stateA->stateB
    '''
    combinationList = list(itertools.combinations_with_replacement(np.arange(orientations)+1, 2))
    timecountDict = {}
    eventcountDict = {}
    # Create bound states keys
    for i in range(boundstates):
        for j in range(boundstates):
            if j != i:
                timecountDict['b' + str(i+1) + '->b' + str(j+1)] = 0.0
                eventcountDict['b' + str(i+1) + '->b' + str(j+1)] = 0
    # Create orientational transition states keys
    # To bound state
    for i in range(boundstates):
        for j in range(len(combinationList)):
            state1, state2 = combinationList[j]
            stateStr = str(10*state1 + state2) + '->b' + str(i+1)
            timecountDict[stateStr] = 0.0
            eventcountDict[stateStr] = 0
    # From bound state
    for i in range(boundstates):
        for j in range(len(combinationList)):
            state1, state2 = combinationList[j]
            stateStr = 'b' + str(i+1) + '->' + str(10*state1 + state2)
            timecountDict[stateStr] = 0.0
            eventcountDict[stateStr] = 0
    return timecountDict, eventcountDict


def extractRatesMSM(discreteTrajectories, lagtime, boundstates, dtEffective, stitching = False, fullDictionary = False):
    '''
    Construct an MSM and extract rates from it using pyemma and msmtools. These rates differ from the rates in a
    rate transition matrix due to the fact that this takes all possible transmission pathways within a discrete
    time MSM (see transition path theory (TPT) to understand one possible way of obtaing these mfpts).
    :param discreteTrajectories: raw discrete trajectories including the unbound state.
    :param lagtime: preferred lagtime chose for MSM (might need tweaking for each case)
    :param boundstates: numbe of boundstates
    :param dtEffective takes into account the effective time step from the data, usually equalt to dt*stride of
    the simulation.
    :return rateDict: dictionary that maps transition with rates. The dictionary keys have the form
    "stateA->stateB"; if the state is a bound state, the key is a "b" followed by the state number.
    For the transition states composed by two integers, they correspond to the two closest orientational
    discrete states of each of the two particles (touched by a line between the two centers of mass). Also
    note rates need to be scale by dt.
    '''
    unboundStateIndex = 0
    # Slice trajectories getting rid of the unbound state 0
    slicedDtrajs = splitDiscreteTrajs(discreteTrajectories, unboundStateIndex)
    if stitching:
        finalTrajs = stitchTrajs(slicedDtrajs, 500)
    else:
        finalTrajs = slicedDtrajs
    # Create MSM between transision states and bound states without stitching
    mainmsm = pyemma.msm.estimate_markov_model(finalTrajs, lagtime, reversible=False)
    # The active set keep track of the indexes used by pyemma and the ones used to describe the state in our model.
    activeSet = mainmsm.active_set
    rateDict = {}
    # Loop over active set, where "i" is the index in pyemma for state "state".
    for i, originState in enumerate(activeSet):
        for j, transitionState in enumerate(activeSet):
            if i != j:
                originKey = str(originState)
                transitionKey = str(transitionState)
                if (originState in list(range(1,boundstates + 1))):
                    originKey = 'b' + originKey
                if (transitionState in list(range(1,boundstates + 1))):
                    transitionKey = 'b' + transitionKey
                key = originKey + '->' + transitionKey
                # These two are equivalent statements
                #mfpt = msmtools.analysis.mfpt(mainmsm.transition_matrix, i, j, tau = lagtime)
                mfpt = mainmsm.mfpt(i, j) # already takes lagtime into consideration
                rateDict[key] = 1.0/(dtEffective*mfpt) # Scaling by dt*stride included here.
    outputDict = dict(rateDict) #shallow copy
    if not fullDictionary:
        for key in rateDict:
            if 'b' not in key:
                outputDict.pop(key)
    return outputDict



def MSMtoRateDictionary(MarkovStateModel, numBoundStates, dtEffective, fullDictionary = False):
    '''
    Given an MSM calculated by PyEMMA, calculates the rate dictionary
    :param MarkovStateModel: MSM obtained by PyEMMA from the MD trajectories.
    :param numBoundStates number of bound states in the model.
    :param dtEffective takes into account the effective time step from the data, usually equalt to dt*stride of
    the simulation.
    :param fullDictionary: if true outputs rates for all possible transitions. Otherwise
    only the ones that include a bound state.
    :return: :return rateDict: dictionary that maps transition with rates. The dictionary keys have the form
    "stateA->stateB"; if the state is a bound state, the key is a "b" followed by the state number.
    For the transition states composed by two integers, they correspond to the two closest orientational
    discrete states of each of the two particles (touched by a line between the two centers of mass). Also
    note rates need to be scale by dt.
    '''
    # The active set keep track of the indexes used by pyemma and the ones used to describe the state in our model.
    activeSet = MarkovStateModel.active_set
    rateDict = {}
    # Loop over active set, where "i" is the index in pyemma for state "state".
    for i, originState in enumerate(activeSet):
        for j, transitionState in enumerate(activeSet):
            if i != j:
                originKey = str(originState)
                transitionKey = str(transitionState)
                if (originState in list(range(1,numBoundStates + 1))):
                    originKey = 'b' + originKey
                if (transitionState in list(range(1,numBoundStates + 1))):
                    transitionKey = 'b' + transitionKey
                key = originKey + '->' + transitionKey
                # The next two are equivalent statements
                #mfpt = msmtools.analysis.mfpt(mainmsm.transition_matrix, i, j, tau = lagtime)
                mfpt = MarkovStateModel.mfpt(i, j) # already takes lagtime into consideration
                rateDict[key] = 1.0/(dtEffective*mfpt) # Scaling by dt*stride included here.
    outputDict = dict(rateDict) #shallow copy
    if not fullDictionary:
        for key in rateDict:
            if 'b' not in key:
                outputDict.pop(key)
    return outputDict




def listIndexSplit(inputList, *args):
    '''
    Function that splits inputList into smaller list by slicing in the indexes given by *args.
    :param inputList:
    :param args: int indexes where list should be splitted (Note to convert a
    list "mylist" into *args just do: *mylist)
    :return: list of sliced lists
    If extra arguments were passed prepend the 0th index and append the final
    # index of the passed list, in order toa v check for oid checking the start
    # and end of args in the loop. Also, add one in args for correct indexing.
    '''
    if args:
        args = (0,) + tuple(data+1 for data in args) + (len(inputList)+1,)
    # Slice list and return list of lists.
    slicedLists = []
    for start, end in zip(args, args[1:]):
        slicedLists.append(inputList[start:end-1])
    if slicedLists == []:
        slicedLists.append(inputList)
    return slicedLists

def splitDiscreteTrajs(discreteTrajs, unboundStateIndex = 0):
    '''
    Splits trajectories into smaller trajectories by cutting out
    all the states unboundStateindex (0)
    :param discreteTrajs: list of discrete trajectories
    :param unboundStateIndex: index of the unbound state used to
    decide where to cut the trajectories, normally we choose it to be
    zero.
    :return: List of sliced trajectories
    '''
    slicedDtrajs = []
    trajnum = 0
    for dtraj in discreteTrajs:
        # Slice trajectory using zeros as reference point
        indexZeros = np.where(dtraj==unboundStateIndex)
        slicedlist = listIndexSplit(dtraj, *indexZeros[0])
        # Remove the empty arrays
        for array in slicedlist:
            if array.size > 1:
                slicedDtrajs.append(array)
        trajnum += 1
        print("Trajectory ", trajnum, " of ", len(discreteTrajs), " done.", end="\r")
    return slicedDtrajs

def stitchTrajs(slicedDtrajs, minlength = 1000):
    '''
    Joins splitted trajectories into long trajectories of at least minlength. The trajectories that cannot
    be joined are left as they were.
    :param slicedDtrajs: list of discrete trajectories. Each discrete trajectory is a numpy array
    :param minlength: minimum length of stitched trajectory if any stititching is possible
    :return: list of stitched trajectories
    '''
    myslicedDtrajs = slicedDtrajs.copy()
    stitchedTrajs = []
    # Stitch trajectories until original sliced trajectory is empty
    while len(myslicedDtrajs) > 0:
        traj = myslicedDtrajs[0]
        del myslicedDtrajs[0]
        # Keep trajectories under a certain length min length
        while traj.size <= minlength:
            foundTrajs = False
            for i in reversed(range(len(myslicedDtrajs))):
                # If end point and start point match, join trajectories
                if traj[-1] == myslicedDtrajs[i][0]:
                    foundTrajs = True
                    traj = np.concatenate([traj, myslicedDtrajs[i]])
                    del myslicedDtrajs[i]
            # If no possible trajectory to join is found, save trajectory and continue.
            if foundTrajs == False:
                break;
        stitchedTrajs.append(traj)
    return stitchedTrajs