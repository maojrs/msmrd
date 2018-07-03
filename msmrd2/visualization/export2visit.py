import numpy as np
from pyevtk.hl import imageToVTK
# Create file from trajectory data to visualize 3D density plot in visit

def createDensityData(trajs, filename = "densityData", radius=2.0):
    MSMradius = radius
    X = np.arange(-MSMradius, MSMradius, 0.05)
    Y = np.arange(-MSMradius, MSMradius, 0.05)
    Z = np.arange(-MSMradius, MSMradius, 0.05)
    Zfull = np.zeros([X.shape[0]-1, Y.shape[0]-1, Z.shape[0]-1])
    for traj in trajs:
        hist = np.histogramdd(traj, bins = (X, Y, Z))
        Zfull += hist[0]
    
    # Dimensions 
    nx, ny, nz = X.shape[0]-1, Y.shape[0]-1, Z.shape[0]-1 
    ncells = nx * ny * nz 
    npoints = (nx + 1) * (ny + 1) * (nz + 1) 

    # Variables 
    filename = "../data/visit/" + filename
    imageToVTK(filename, cellData = {"Density" : Zfull})
