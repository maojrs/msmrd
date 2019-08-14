import h5py
import numpy as np

# Functions to load trajectories and manipulate them

def loadTrajectory(fnamebase, fnumber, fastload = False):
    '''
    Reads data from discrete trajectory and returns a simple np.array of
    integers representing the discrete trajectory. Assumes h5 file.
    :param fnamebase, base of the filename
    :param fnumber, filenumber
    :param fastload if true loads the H5 data, if false it converts the data to numpy.
    this however makes the loading very slow.
    :return: array of arrays representing the trajectory
    '''
    filename = fnamebase + str(fnumber).zfill(4) + '.h5'
    f = h5py.File(filename, 'r')

    # Get the data
    a_group_key = list(f.keys())[0]
    if fastload:
        data = f[a_group_key]
    else:
        data = np.array(f[a_group_key])
        #data = f[a_group_key][:] # equivalent

    return data



def loadDiscreteTrajectory(fnamebase, fnumber, fnamesuffix = '_discrete', filetype = 'h5'):
    '''
    Reads data from discrete trajectory and returns a simple np.array of
    integers representing the discrete trajectory. The file can be in the
    h5 or xyz format.
    :param fnamebase, base of the filename
    :param fnumber, filenumber
    :param fnamesuffix, suffix added at end of filename before the extension
    :param filetype, string indicating which format, h5 or xyz, is the file
    :return: array with integers representing the discrete trajectory
    '''
    if filetype == 'h5':
        filename = fnamebase + str(fnumber).zfill(4) + fnamesuffix + '.h5'
        f = h5py.File(filename, 'r')

        # Get the data
        a_group_key = list(f.keys())[0]
        data = np.array(f[a_group_key])

        # Transform data into simple 1D array
        array = np.zeros(len(data), dtype = int)
        for i in range(len(data)):
            array[i] = data[i][0]
        return array

    if filetype == 'xyz':
        filename = fnamebase + str(fnumber).zfill(4) + fnamesuffix + '.xyz'
        file = open(filename, "r")

        # Read file and save to array
        filelines = file.readlines()
        array = np.zeros(len(filelines), dtype = int)
        for i, line in enumerate(filelines):
            array[i] = int(float(line))
        return array



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