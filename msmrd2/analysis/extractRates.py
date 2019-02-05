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
    in this case: 0-unbound, 1-first bound state, 2-second bound state, ij orientation transition state (patchyDimer)
    :param discreteTrajectories: list of discrete trajectories, i.e. each element is one discrete trajectory
    :return: {timecountDict, eventcountDict} dictionaries that maps state label to accumulated time
    to transition and to number of events found for that particular transition. The dictionary keys
    have the form stateA->stateB
    '''
    for dtraj in discreteTrajectories:
        # Loop over one trajectory values
        for i in range(len(dtraj)-1):
            
            if (dtraj[i+1] == 1 and dtraj[i] != 1 and dtraj[i] != 0):
                state = dtraj[i]
                prevstate = state
                tstep = 1
                while (prevstate == state):
                    prevstate = dtraj[i-tstep]
                    if (prevstate == state):
                        tstep += 1
                # Update dictionaries
                if (state == 1 or state == 2) :
                    dictkey = 'b' + str(state)
                else:
                    dictkey = str(state)
                timecountDict[dictkey] += tstep
                eventcountDict[dictkey] += 1
    return timecountDict, eventcountDict