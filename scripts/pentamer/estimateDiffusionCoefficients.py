import numpy as np
import matplotlib.pyplot as plt
import msmrd2
import msmrd2.tools.quaternions as quats
from msmrd2.potentials import patchyParticleAngular2
from msmrd2.integrators import overdampedLangevin as odLangevin

'''
Estimates translational and rotational diffusion coefficients of pentamer or parts of the pentamer with a 
given number of particles (between 1 and 5) Note that if ran with numparticles=1, we recover (with some error)
the diffusion coefficients of the particles.
'''

# Define global parameters
boxsize = 3
D = 1.0
Drot = 1.0
np.random.seed(seed=1) # seed 1 good for pentamer

# pentamer parameters
th = 2*np.pi/5.0
thextra = np.pi

# potential parameters
sigma = 1.0
strength = 160 #100# 60 #200.0
angularStrength = 20 #10 #200.0
angleDiff = 3*np.pi/5.0
patch1 = np.array([np.cos(angleDiff/2),np.sin(angleDiff/2),0.])
patch2 = np.array([np.cos(-angleDiff/2),np.sin(-angleDiff/2),0.])
patchesCoordinates = [patch1, patch2]

# Integration paramaters
dt = 0.00001
timesteps = 1000000 # CAREFUL, if simulation time too large, particles might unbind
stride = 10 #25 #250 #1000
# Define reference vectors for plotting reference system
refVec4 = np.array([0., 0., 3.])
refVec5 = np.array([0., 0., -3.])
# output for vmd plotting to check configurations
outputVMD = False

def calculateDiffusionCoefficeints(numparticles):
    '''
    numparticles = 1 # It
    :param numparticles: must have a value between 1 and 5
    :return: translational diffusion coefficeint and rotational diffusion coefficient
    '''
    simulationSuccess = False

    # Define particle list with between two and five particles to estimate the diffusion
    # of dimers/ trimers / cuatrimers? and pentamer
    pyPartlist = []

    # Position List for pentamer IC
    for i in range(numparticles):
        position = 0.85 * np.array([np.cos(th*i),np.sin(th*i),0.0])
        orientation = np.array([np.cos(0.5*(thextra - th*i)),0,0,np.sin(0.5*(thextra - th*i))])
        part = msmrd2.particle(D, Drot, position, quats.conjugate(orientation))
        pyPartlist.append(part)

    # Create list of particles that can be read from msmrd
    # Note the particles in this list will be independent from the python list.
    partlist = msmrd2.integrators.particleList(pyPartlist)

    # Over-damped Langevin integrator definition
     #0.000005
    seed = -1 #1 #-1 # Negative seed, uses random device as seed
    bodytype = 'rigidbody'
    integrator = odLangevin(dt, seed, bodytype)

    # Define Patchy Particle potential
    if numparticles > 1:
        potentialPatchyParticleAngular2 = patchyParticleAngular2(sigma, strength, angularStrength, patchesCoordinates)
        integrator.setPairPotential(potentialPatchyParticleAngular2)

    # Define arrays to calculate autocorrelation functions
    pentamerPositionArray = []
    pentamerOrientationArray = []

    #Integrate particle list and print only positions
    cross = [None]*2
    if outputVMD:
        datafile  = open('../../data/vmd/pentamerTest.xyz', 'w')
    for i in range(timesteps):
        if i%stride == 0:
            if outputVMD:
                datafile.write(str(3*len(partlist) + 3) + '\n')
                datafile.write(str(0) + '\n')
            if numparticles == 1:
                pentamerPositionArray.append(partlist[0].position)
                pentamerOrientationArray.append(partlist[0].orientation)
            else:
                # Calculate and plot center of pentamer position
                relpos = partlist[1].position - partlist[0].position
                relposnorm = np.linalg.norm(relpos)
                if (relposnorm >= 1.5):
                    return simulationSuccess, 0, 0
                relpos = relpos/relposnorm
                patch1 = 0.5*sigma*quats.rotateVec(patchesCoordinates[0], partlist[0].orientation)
                patch2 = 0.5*sigma*quats.rotateVec(patchesCoordinates[1], partlist[0].orientation)
                cross[0] = np.cross(relpos,patch1)
                cross[1] = np.cross(relpos,patch2)
                maxIndex = np.argmax(np.linalg.norm(cross, axis=1))
                rotAxis = cross[maxIndex]/np.linalg.norm(cross[maxIndex])
                rotation = 3*np.pi*rotAxis/10.0
                quatRotation = quats.angle2quat(rotation)
                pentamerCenter = 0.85*quats.rotateVec(relpos,quatRotation)
                pentamerCenter = pentamerCenter + partlist[0].position
                pentamerPositionArray.append(pentamerCenter)
                # Calculate and plot orientation of pentamer (using only the orientation of particle 0)
                orientation0 = np.array([np.cos(0.5*(thextra)),0,0,np.sin(0.5*(thextra))])
                vec1 = 0.85 * np.array([1.,0.,0.]) + 0.5*sigma*quats.rotateVec(patchesCoordinates[0],orientation0)
                vec2 = 0.85 * np.array([1.,0.,0.]) + 0.5*sigma*quats.rotateVec(patchesCoordinates[1],orientation0)
                rotVec1 = partlist[0].position + 0.5*sigma*quats.rotateVec(patchesCoordinates[0],partlist[0].orientation)
                rotVec2 = partlist[0].position + 0.5*sigma*quats.rotateVec(patchesCoordinates[1],partlist[0].orientation)
                pentamerOrientation = quats.recoverRotationFromVectors(np.array([0.,0.,0.]), vec1, vec2, \
                                                                       pentamerCenter, rotVec1, rotVec2)
                pentamerOrientationArray.append(pentamerOrientation)
                if outputVMD:
                    # Calculate cross reference
                    v0 = pentamerCenter
                    v1 = v0 + 0.25*quats.rotateVec(refVec4, pentamerOrientation)
                    v2 = v0 + 0.25*quats.rotateVec(refVec5, pentamerOrientation)
                    # Plot cross reference
                    datafile.write('type_2' + ' ' + ' '.join(map(str, v0)) + '\n')
                    datafile.write('type_3' + ' ' + ' '.join(map(str, v1)) + '\n')
                    datafile.write('type_3' + ' ' + ' '.join(map(str, v2)) + '\n')
        if outputVMD:
            for j, part in enumerate(partlist):
                if i%stride == 0:
                    v0 = part.position
                    v1 = v0 + 0.5*sigma*quats.rotateVec(patchesCoordinates[0], part.orientation)
                    v2 = v0 + 0.5*sigma*quats.rotateVec(patchesCoordinates[1], part.orientation)
                    datafile.write('type_0' + ' ' + ' '.join(map(str, v0)) + '\n')
                    datafile.write('type_1' + ' ' + ' '.join(map(str, v1)) + '\n')
                    datafile.write('type_1' + ' ' + ' '.join(map(str, v2)) + '\n')
        integrator.integrate(partlist)
        if i%10000 == 0:
            print("Percentage complete: ", 100*i/timesteps, "%", end="\r")
    if outputVMD:
        datafile.close()
        # Generate TCL script to visualize with VMD
        msmrdvis.generateTCL_pentamerTest(numparticles = numparticles,
                                          outfname = 'pentamerTest',
                                          tclfname = '../../data/vmd/pentamerTest_2vmd.tcl')
    print("Simulation complete: ", 100, " % w/", numparticles, " particles")

    # Calculate auto-correlation/mean-square displacement for several lagtimes
    tsteps = len(pentamerPositionArray)
    lagtimesIndexes = np.arange(100,2000,100) #np.arange(1,20,1)
    lagtimes = np.zeros(len(lagtimesIndexes) + 1)
    MSD_3D = np.zeros(len(lagtimesIndexes) + 1)
    MSD_3D_orientation = np.zeros(len(lagtimesIndexes) + 1)
    for i, lagtimeIndex in enumerate(lagtimesIndexes):
        lagtimes[i+1] = dt * lagtimeIndex * stride
        MSD = 0.0
        MSD_orientation = 0.0
        for j in range(tsteps-lagtimeIndex):
            dr = pentamerPositionArray[j+lagtimeIndex] - pentamerPositionArray[j]
            dq = pentamerOrientationArray[j+lagtimeIndex] - pentamerOrientationArray[j]
            MSD += dr*dr
            MSD_orientation += dq*dq
        MSD = MSD/(tsteps - lagtimeIndex + 1)
        MSD_3D[i+1] = sum(MSD) # D = sum(MSD)/(6*lagtime)
        MSD_orientation = MSD_orientation/(tsteps - lagtimeIndex + 1)
        MSD_3D_orientation[i+1] = sum(MSD_orientation[1:])

        # Least square approximation with numpy
        A = np.vstack([lagtimes, np.ones(len(lagtimes))]).T
        slope,b = np.linalg.lstsq(A, MSD_3D/6, rcond=None)[0]
        #print(slope,b)
        # Least square approximation with numpy (rotation)
        y = - np.log(1-4*MSD_3D_orientation/3.0)/2
        slope2,b2 = np.linalg.lstsq(A, y, rcond=None)[0]
        #print(slope2,b2)
    print("Estimated diffusion coefficients")

    # Plot MSD against lagtime and fits
    plt.plot(lagtimes, MSD_3D/6, 'o', label = 'position')
    plt.plot(lagtimes, - np.log(1-4*MSD_3D_orientation/3.0)/2, 'o', label='orientation')
    plt.plot(lagtimes, slope*lagtimes + b, '-', label = 'position fit')
    plt.plot(lagtimes, slope2*lagtimes + b2, '-', label = 'orientation fit')
    plt.legend()
    plt.savefig('../../data/pentamer/diffusion_coefficients/diffusionFitPentamer_' + str(numparticles) + 'particles.pdf', bbox_inches='tight')
    plt.close()
    print("Finished plots")

    Dapprox = slope
    DrotApprox = slope2
    simulationSuccess = True
    return simulationSuccess, Dapprox, DrotApprox

datafile2  = open('../../data/pentamer/diffusion_coefficients/diffusionCoefficients.dat', 'w')
for i in range(5):
    succesful = False
    while (not succesful):
        succesful, D, Drot = calculateDiffusionCoefficeints(i+1)
        if succesful == False:
            print("Unbinding event, repeating simulation")
        else:
            print('Diffusion coefficients for ', i+1, ' particles: ', D, Drot)
            datafile2.write(str(i+1) + ' ' + str(D) + ' ' + str(Drot) + '\n')
datafile2.close()