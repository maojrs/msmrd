

# Tools to extract rates from Markov state models (MSMs) and other tools related to MSMs

def MSMtoRateDictionary(MarkovStateModel, numBoundStates, dtEffective, fullDictionary = False):
    '''
    Given an MSM calculated by PyEMMA, calculates the rate dictionary
    :param MarkovStateModel: MSM obtained by PyEMMA from the MD trajectories. If not calculated by pyemma,
    it should be an object that contains an active_set (keep tracks of the how indexes in the MSM relate to
    state index in the model), and it should have a method mfpt(i,j), that calculates the mean first passage
    time from state i to state j (i, j, MSM indexes). This should calculate it the TPT way, i.e. passing
    through any of the other states.
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