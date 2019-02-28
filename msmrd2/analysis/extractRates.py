import h5py
import numpy as np
import msmrd2.tools as msmrdtls
import itertools


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
    :param orientations: number of discrete orientations (e.g. 3), which determine the number of possible
    transition states between relative orientations (e.g. 11, 12, 13, 22, 23, 33).
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


def extractRates(discreteTrajectories, timecountDict, eventcountDict):
    '''
    Calculates rates from discrete trajectories. Verify convention used with corresponding trajectory class,
    in this case: 0-unbound, 1-first bound state, 2-second bound state, 3 is a transition state, not relevant for
    this calculation, ij orientation transition state (patchyDimer)
    :param discreteTrajectories: list of discrete trajectories, i.e. each element is one discrete trajectory
    :return: {rateDict, timecountDict, eventcountDict} the second two dictionaries map state label to accumulated time
    to transition and to number of events found for that particular transition. The first dictionary returns the
    rates corresponding to each transition. The dictionary keys have the form stateA->stateB; if the state is a bound
    state, the key is a "b" followed by the state number. For the transision states composed by two integers,
    correspond to the two closest orientational discrete states of each of the two particles (touched by a line
    between the two centers of mass).
    '''
    for dtraj in discreteTrajectories:
        # Loop over one trajectory values
        for i in range(len(dtraj)-1):
            # Make sure there is a transition and that neither the current state nor the endstate are zero
            state = dtraj[i]
            endstate = dtraj[i+1]
            if (state != endstate) and (state != 0) and (endstate != 0) and (state != 3) and (endstate != 3):
                # If neither the current state nor the end state is one, don't store transition (skip one cycle in loop).
                if (state != 1 and endstate != 1 and state != 2 and endstate !=2):
                    continue
                # Count how many timesteps the trajectory remained in current state before transitioning to endstate
                prevstate = state
                tstep = 1
                while (prevstate == state):
                    prevstate = dtraj[i-tstep]
                    if (prevstate == state):
                        tstep += 1
                # Update number of events and accumulated time (timesteps) in dictionaries
                if (state == 1 or state == 2) :
                    dictkey1 = 'b' + str(state)
                else:
                    dictkey1 = str(state)
                if (endstate ==1 or endstate == 2):
                    dictkey2 = 'b' + str(endstate)
                else:
                    dictkey2 = str(endstate)
                dictkey = dictkey1 + '->' + dictkey2
                timecountDict[dictkey] += tstep
                eventcountDict[dictkey] += 1
    # Scale by the number of events and save in rate Dictionary (note rates need to be scaled by dt)
    rateDict = {}
    for key in timecountDict:
        if eventcountDict[key] != 0:
            rateDict[key] = timecountDict[key]/eventcountDict[key]

    return rateDict, timecountDict, eventcountDict


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
    return slicedLists

def splitDiscreteTrajs(discreteTrajs, unboundStateIndex = 0):
    '''
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