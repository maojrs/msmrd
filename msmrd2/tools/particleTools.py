''' Library of functions to do operations with particles '''

import numpy as np
import msmrd2
import warnings



def randomParticleList(numparticles, boxsize, separationDistance, D, Drot):
    '''
    :param numparticles: number of particles in list
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


def randomParticleMSList(numparticles, boxsize, separationDistance, types, unboundMSMs):
    '''
    :param numparticles: number of particles in list
    :param boxsize: size of simulation box, if scalar it assumes the three box edges are the same in all dimensions
    :param separationDistance: separation distance between particle to avoid overlapping
    :param types types of particles (can be array of integers or integer if all particles are the same)
    :param unboundMSMs: list of unboundMSM, needed to extract diffusion coefficients of particles.
    :param state initial state of particles (can be array of integers or integer if all particles are the same)
    :return:
    Generates list of particleMS with uniformly random position and orientation. It also enforces the particles don't
    overlap over a distance given by the relative distance cut off.
    '''

    # Warning in case of not ideal parameters.
    if separationDistance > boxsize:
        warnings.warn("The separation distance between particles should be much smaller than the "
                      "boxsize. Otherwise risk of infinite loop.")

    # Transform boxsize, diffusion coefficients, type and state to vectors if neccesarry.
    if np.isscalar(boxsize):
        boxsize = np.array([boxsize, boxsize, boxsize])
    if np.isscalar(types):
        newtypes = np.ones(numparticles, dtype = 'int')
        types = (types*newtypes).astype(int)

    # Obtain initial state randomly from all available states
    states = np.zeros(numparticles, dtype = 'int')
    for i in range(numparticles):
        maxstates = len(unboundMSMs[types[i]].D)
        states[i] = np.random.randint(maxstates)

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
        D = unboundMSMs[types[i]].D[states[i]]
        Drot = unboundMSMs[types[i]].Drot[states[i]]
        part = msmrd2.particleMS(types[i], states[i], D, Drot, position, orientation)
        # Deactivate inner unbound MSM if there is only one state.
        if len(unboundMSMs[types[i]].D) == 1:
            part.deactivateMSM()
        pyPartlist.append(part)

    partlist = msmrd2.integrators.particleMSList(pyPartlist)
    return partlist