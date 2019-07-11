import numpy as np

# Functions to help analyze data

def bootstrappingMFPTs(fptList, numBootstrapSamples):
    '''
    Computes the bootsrapping statistic of a list of first passage time.
    :param fptList: list or array of first passage times
    :param numBootstrapSamples: number of bottstrapped samples used to calculate the
    mean first passage time and its standard deviation
    :return: mean first passage time and its standard deviation
    '''
    numValues = len(fptList)
    mfpts = np.zeros(numBootstrapSamples)
    # Calculate boostrapped samples
    for i in range(numBootstrapSamples):
        bootstrapSample = np.random.choice(fptList, size = numValues, replace = True)
        mfpts[i] = bootstrapSample.mean()
    # Compute mean and standard deviation over bootstrapped samples
    mfpt = mfpts.mean()
    fptStdDev = np.std(mfpts)
    return mfpt, fptStdDev

