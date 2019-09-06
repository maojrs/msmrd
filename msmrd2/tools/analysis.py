import numpy as np
import pickle

# Functions to help analyze data

def bootstrapping(valuesList, numBootstrapSamples, numValues = None):
    '''
    Computes the bootsrapping statistic of a list of arbitrary values.
    :param valuesList: list or array of values to calculate bootstrap statistics on
    :param numBootstrapSamples: number of bootstrapped samples used to calculate the
    mean value and its standard deviation
    :return: mean value and its standard deviation
    '''
    if numValues == None:
        numValues = len(valuesList)
    values = np.zeros(numBootstrapSamples)
    # Calculate boostrapped samples
    for i in range(numBootstrapSamples):
        bootstrapSample = np.random.choice(valuesList, size = numValues, replace = True)
        values[i] = bootstrapSample.mean()
    # Compute mean and standard deviation over bootstrapped samples
    meanValues = values.mean()
    valuesStdDev = np.std(values)
    return meanValues, valuesStdDev

def writeParameters(filename, parameterDictionary):
    '''
    Writes parameters into .dat file for easy reading and into a .pickle file for
    easy loading.
    :param filename: filename (and location) to write parameters
    :param parameterDictionary: dictionary of parameters to be written
    '''
    with open(filename + ".dat", 'w') as file:
        for key, value in parameterDictionary.items():
            file.write(key + ' ' + str(value) + '\n')

    pickle_out = open(filename + ".pickle","wb")
    pickle.dump(parameterDictionary, pickle_out)
    pickle_out.close()


def readParameters(filename):
    pickle_in = open(filename + ".pickle","rb")
    parameterDictionary = pickle.load(pickle_in)
    return parameterDictionary
