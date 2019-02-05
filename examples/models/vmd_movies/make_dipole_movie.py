import msmrd2
import msmrd2.visualization as msmrdvis
import numpy as np
import os

# Define continuous-time MSM
MSMtype = 1
ratematrix = np.array([[-3.0,3.0],[1.5,-1.5]])
seed = 0 # Seed = -1 used random device as seed
ctmsm = msmrd2.ctmsm(MSMtype, ratematrix, seed)
Dlist = np.array([0.0, 0.0])
Drotlist = np.array([1.0, 8.0])
ctmsm.setD(Dlist)
ctmsm.setDrot(Drotlist)

# Particle definition
bodytype = 'rod'
orientation = np.array([1,0,0,0])
particlelist = []
# Create a grid of particles
for i in range(6):
    for j in range(6): 
        position = np.array([i-2.5,j-2.5,0])
        state = np.random.randint(0,2)
        part = msmrd2.particleMS(MSMtype, state, Dlist[state], Drotlist[state], bodytype, position, orientation)   
        orientVector = np.random.normal(0,1,3)
        orientVector = orientVector/np.linalg.norm(orientVector)
        part.setOrientVector(orientVector)
        particlelist.append(part)
# Define list object that can be read by pybound functions
partlist = msmrd2.integrators.particleMSList(particlelist)

# Over-damped Langevin integrator with Markovian Switch definition
dt = 0.001
seed = -1 # Seed = -1 used random device as seed
rotation = True
intg = msmrd2.integrators.odLangevinMarkovSwitch(ctmsm, dt, seed, rotation)

# Define and set potential
scalefactor = 30.0
Efieldvector = np.array([0,1,0])
dipolePot = msmrd2.potentials.dipole(scalefactor, Efieldvector)
intg.setExternalRodPotential(dipolePot)

# Integrate the particles and save to .dat to plot into VMD
os.chdir("../../data/vmd/")
datafile  = open('odLangevinDipole.xyz', 'w')
timeIters = 100
prevstate = [None] * len(partlist)
for i in range(timeIters):
    datafile  = open('odLangevinDipole.xyz', 'w')
    datafile.write(str(2*len(partlist)) + '\n')
    datafile.write(str(0) + '\n')
    for j, part in enumerate(partlist):
        if part.state != prevstate[j] and prevstate[j] != None:
            part.setOrientVector(-part.orientvector)
        if part.state == 0:
            v0 = part.position
            v1 = v0 + 0.35*part.orientvector
            v2 = v0 - 0.35*part.orientvector
            datafile.write('type_0' + ' ' + ' '.join(map(str, v1)) + '\n')
            datafile.write('type_0' + ' ' + ' '.join(map(str, v2)) + '\n')
        elif part.state == 1:
            v0 = part.position
            v1 = v0 + 0.35*part.orientvector
            v2 = v0 - 0.35*part.orientvector
            datafile.write('type_1' + ' ' + ' '.join(map(str, v1)) + '\n')
            datafile.write('type_1' + ' ' + ' '.join(map(str, v2)) + '\n')
        prevstate[j] = 1*part.state
    datafile.close()
    msmrdvis.generateTCL_dipole(i)
    os.system("vmd -dispdev text -e odLangevinDipole2vmd.tcl")
    intg.integrateList(partlist)

