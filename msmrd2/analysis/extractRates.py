import h5py
import numpy as np
import msmrd2.tools as msmrdtls
import itertools


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
    data = list(f[a_group_key])

    # Transform data into simple 1D array
    array = np.zeros(len(data))
    for i in range(len(data)):
        array[i] = data[i][0]
    return array


def createStatesDisctionaries(boundstates, orientations):
    '''
    Given number of bound states and number of orientation transition states,
    creates two dictionaries to count events and accumulated time for given transition.
    This is required to obtain the rates from the discrete trajectories.
    :param boundstates: number of bound states. They will be labelled as b0, b1, ....
    :param orientations: number of discrete orientations (e.g. 3), which determine the number of possible
    transition states between relative orientations (e.g. 11, 12, 13, 22, 23, 33).
    :return: {timecountDict, eventcountDict} dictionaries that maps state label to accumulated time
    to transition and to number of events found for that particular transition.
    '''
    combinationList = list(itertools.combinations_with_replacement(np.arange(orientations)+1, 2))
    timecountDict = {}
    eventcountDict = {}
    for i in range(boundstates):
        timecountDict['b' + str(i)] = 0.0
        eventcountDict['b' + str(i)] = 0
    for i in range(len(combinationList)):
        state1, state2 = combinationList[i]
        stateStr = str(10*state1 + state2)
        timecountDict[stateStr] = 0.0
        eventcountDict[stateStr] = 0
    return timecountDict, eventcountDict


def extractRates(discreteTrajectories, timecountDict, eventcountDict):
    '''
    Calculates rates from discrete trajcetories. Verify convention used with corresponding trajectory class,
    in this case: 0-unbound, 1-first bound state, 2-second bound state, ij orientation transition state (patchyDimer)
    :param discreteTrajectories: list of discrete trajectories, i.e. each element is one discrete trajectory
    :return: to be determined
    '''
    for dtraj in discreteTrajectories:
        for i in range(len(dtraj)-1):
            if (dtraj[i+1] == 1 and dtraj[i] != 1 and dtraj[i] != 0):
                state = dtraj[i]
                prevstate = state
                tstep = 1
                while (prevstate == state):
                    prevstate = dtraj[i-tstep]
                    if (prevstate == state):
                        tstep += 1
                timecountDict[str(state)] += tstep
                eventcountDict[str(state)] += tstep
    return timecountDict, eventcountDict