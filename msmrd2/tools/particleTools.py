''' Library of functions to do operations with particles '''

import numpy as np
import msmrd2
import warnings



def randomPartList(numparticles, boxsize, separationDistance, D, Drot):
    '''
    :param numparticles: number of particle sin list
    :param boxsize: size of simulation box, if scalar it assumes the three box edges are the same in all dimensions
    :param separationDistance: separation distance between particle to avoid overlapping
    :param D and Drot: diffusion coefficients of particles.
    :return:
    Generates particle list with uniformly random position and orientation. It also enforces the particles don't
    overlap over a distance given by the relative distance cut off.
    '''

    # Warning in case of not ideal parameters.
    if separationDistance > boxsize:
        warnings.warn("The separation distance between particles should be much smaller than the "
                      "boxsize. Otherwise risk of infinite loop.")

    # Transform boxsize to vector if neccesarry.
    if np.isscalar(boxsize):
       boxsize = np.array([boxsize, boxsize, boxsize])

    # Create non-overlapping particle list
    pyPartlist = []
    for i in range(numparticles):
        overlap = True
        while overlap:
            position = np.array([boxsize[0]*np.random.rand()-0.5*boxsize[0],
                                 boxsize[1]*np.random.rand()-0.5*boxsize[1],
                                 boxsize[2]*np.random.rand()-0.5*boxsize[2]])
            overlap = False
            for j in range(len(pyPartlist)):
                if np.linalg.norm(position - pyPartlist[j].position) < separationDistance:
                    overlap = True
                    continue

        orientation = np.array([np.random.rand(),np.random.rand(),np.random.rand(),np.random.rand()])
        orientation = orientation/np.linalg.norm(orientation)
        part = msmrd2.particle(D, Drot, position, orientation)
        pyPartlist.append(part)

    partlist = msmrd2.integrators.particleList(pyPartlist)
    return partlist