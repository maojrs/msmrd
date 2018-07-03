import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def partitionSphere(num_partitions):
    ''' Calculate equal area partition of unit sphere with "num_partitions" partitions. '''

    def angle_to_cap_area(phi):
        return 4*np.pi*(np.sin(phi/2.0))**2

    def cap_area_to_angle(area):
        return 2.0*np.arcsin(np.sqrt(area/(4.0*np.pi)))

    # Calculate areas of each state and polar caps angle (phi0 and pi-phi0)
    state_area = 4*np.pi/num_partitions
    phi0 = cap_area_to_angle(state_area)
    # Calculate the number of collars between the polar caps
    ideal_collar_angle = np.sqrt(state_area)
    ideal_num_collars = (np.pi - 2*phi0)/ideal_collar_angle
    num_collars = int(max(1,np.round(ideal_num_collars)))
    if num_partitions == 2:
        num_collars = 0
    collar_angle = (np.pi - 2*phi0)/num_collars
    # Initialize variables for number of regions in each collar
    ideal_regionsPerCollar = np.zeros(num_collars)
    phis = np.zeros(num_collars + 1)
    regionsPerCollar = np.zeros(num_collars)
    thetas = []
    a = [0]
    # Iterate over each collar to get right number of regions per collar
    # and correct location of phi angles of each collar
    for i in range(num_collars):
        # Calculate num of regions in collar i
        cap_area_phi1 = angle_to_cap_area(phi0 + i*collar_angle)
        cap_area_phi2 = angle_to_cap_area(phi0 + (i+1)*collar_angle)
        ideal_regionsPerCollar[i] = (cap_area_phi2 - cap_area_phi1)/state_area
        regionsPerCollar[i] = np.round(ideal_regionsPerCollar[i] + a[i])
        # Correct values of phi around collar i
        suma = 0
        for j in range(i+1):
            suma = suma + ideal_regionsPerCollar[j] - regionsPerCollar[j]
        a.append(suma)
        summ = 1
        for j in range(i):
            summ = summ + regionsPerCollar[j]
        phis[i] = cap_area_to_angle(summ*state_area)
        phis[-1] = np.pi - phi0
        # Obtain list of thetas for a given collar
        thetasi = np.zeros(int(regionsPerCollar[i]))
        dth = 2.0*np.pi/regionsPerCollar[i]
        for j in range(int(regionsPerCollar[i])):
            thetasi[j] = j*dth
        thetas.append(thetasi)
    phis = np.insert(phis,0,0)
    regionsPerCollar = np.append(regionsPerCollar,1)
    regionsPerCollar = np.insert(regionsPerCollar,0,1)
        # return number of regions for all collars, 
        # phi angles of collars and theta angles for each collar
    return regionsPerCollar.astype(np.int32), phis, thetas


def plotPartitionedSphere(numPartitions = None, save = None, plotState = None, coord = None):
    '''
    Plots partitioned sphere function for a given numPartitions,
    save: set to True saves the output figure to file,
    plotState: choose a state between 1 and numPartitions to mark that state. Note when
    plotState is set, no inner sphere is plotted, and full lines are shown.
    coord: takes list/array of coordinates and plots the instersection between the unit sphere and the ray from
    the origin to coord. Used to test that correct states are returned.
    '''
    if numPartitions == None:
        numPartitions == 20
    if save == None:
        save = False
    if plotState != None:
        if plotState > numPartitions or plotState <= 0:
            print("State to plot out of range")
            return
    numRegionsCollar, phis, thetas = partitionSphere(numPartitions)

    fig = plt.figure(figsize=plt.figaspect(0.95)*1.5)
    ax = fig.gca(projection='3d')
    ax._axis3don = False
    # Adjust for proper viewing with mplot3d (lines behind sphere not plotted)
    if plotState == None:
        minth = 0 
        maxth = np.pi
    else:
        minth = 0 
        maxth = 2*np.pi

    # Plot inner white sphere
    r=1
    if plotState == None:
        u = np.linspace(0, 2 * np.pi, 400)
        v = np.linspace(0, np.pi, 400)
        xx = r * np.outer(np.cos(u), np.sin(v))
        yy = r * np.outer(np.sin(u), np.sin(v))
        zz = r * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(xx, yy, zz, color='white', linewidth=0, antialiased=False, alpha = 1.0)

    # Plot collars
    for phi in phis:
        theta = np.linspace(minth,maxth, 50)
        x = r * np.cos(theta) * np.sin(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(phi) 
        ax.plot(x, y, z, '-k')
    # Plot divisions in every collar
    for i in range(len(numRegionsCollar)):
        numDiv = int(numRegionsCollar[i])
        if numDiv > 1:
            dth = 2 * np.pi /numDiv
            phi = np.linspace(phis[i],phis[i+1],10)
            for j in range(numDiv):
                theta = j*dth
                if (theta >= minth) and (theta <= maxth):
                    x = r * np.cos(theta) * np.sin(phi)
                    y = r * np.sin(theta) * np.sin(phi)
                    z = r * np.cos(phi)
                    ax.plot(x, y, z, '-k')
    
    # Find location of state in plotState (counts from one onwards)
    # First find collar (phis)
    if plotState != None:
        collar = 0
        while (sum(numRegionsCollar[0:collar+1])) < plotState:
            collar += 1
        phi1 = phis[collar]
        if collar+1 < len(numRegionsCollar):
            phi2 = phis[collar+1]
        else:
            phi2 = np.pi
        # Find thetas
        prevStates = sum(numRegionsCollar[0:collar])
        if ((prevStates == 0) or (prevStates == numPartitions -1)):
            theta1 = 0
            theta2 = 2*np.pi
        else:
            thetasCollar = thetas[collar-1]
            statesInCollar = plotState - prevStates
            theta1 = thetasCollar[statesInCollar-1]
            if statesInCollar == len(thetasCollar):
                theta2 = 2*np.pi
            else:
                theta2 = thetasCollar[statesInCollar]
        # Plot special state
        if min(phi1,phi2) == 0 or max(phi1,phi2) == np.pi:
            thetas = np.linspace(theta1,theta2, 50)
            if phi1 == 0:
                phis = phi2
            else:
                phis = phi1
            xx = r * np.cos(thetas) * np.sin(phis)
            yy = r * np.sin(thetas) * np.sin(phis)
            zz = r * np.cos(phis)
            ax.plot(xx, yy, zz, '-r', lw=3)
        else:
            thetas = np.linspace(theta1,theta2, 50)
            phis = np.linspace(phi1,phi2, 50)
            x1 = r * np.cos(thetas) * np.sin(phi1)
            y1 = r * np.sin(thetas) * np.sin(phi1)
            z1 = r * np.cos(phi1)
            x2 = r * np.cos(thetas) * np.sin(phi2)
            y2 = r * np.sin(thetas) * np.sin(phi2)
            z2 = r * np.cos(phi2)
            x3 = r * np.cos(theta1) * np.sin(phis)
            y3 = r * np.sin(theta1) * np.sin(phis)
            z3 = r * np.cos(phis)
            x4 = r * np.cos(theta2) * np.sin(phis)
            y4 = r * np.sin(theta2) * np.sin(phis)
            z4 = r * np.cos(phis)
            ax.plot(x1, y1, z1, '-r', lw=3)
            ax.plot(x2, y2, z2, '-r', lw=3)
            ax.plot(x3, y3, z3, '-r', lw=3)
            ax.plot(x4, y4, z4, '-r', lw=3)
    
    # plot point of coord
    if coord != None:
        coordnorm = np.linalg.norm(coord)
        xc = coord[0]/coordnorm
        yc = coord[1]/coordnorm
        zc = coord[2]/coordnorm
        ax.plot([xc],[yc],[zc],'or')
                         
    # Plot the surface
    #ax.set_aspect('equal')
    ax.view_init(0,270)
    ax.dist = 5.65
    if save:
        plt.savefig('spherePartion_' + str(numPartitions) + '.png')
        


def getSectionNumber(coords, numPartitions = 20):
    '''
    Given some coordinates, returns corresponding section number in
    spherical partition (from 1 to numPartitions)
    '''
    numRegionsCollar, phis, thetas = partitionSphere(numPartitions)
    theta = np.arctan2(coords[1], coords[0])
    if theta < 0:
        theta += 2*np.pi
    r = np.linalg.norm(coords)
    phi = np.arccos(coords[2]/r)
    currentCollarIndex = sum(phis<=phi) - 1
    if currentCollarIndex == 0:
        sectionNum = 1
        return sectionNum
    if currentCollarIndex == len(numRegionsCollar) - 1:
        sectionNum = numPartitions
        return sectionNum
    collarThetas = thetas[currentCollarIndex - 1]
    currentThetaIndex = sum(collarThetas<=theta)
    sectionNum = sum(numRegionsCollar[0:currentCollarIndex]) + currentThetaIndex
    return sectionNum


def getAngles(secNumber, numPartitions = None):
    '''
    Return angles [phi1, phi2] and [theta1, theta2] of a given section number
    in the spherical partition
    '''
    if numPartitions == None:
        numPartitions == 20
    if secNumber > numPartitions:
        print("Error: section number is larger than number of partitions")
        return
    numRegionsCollar, phis, thetas = partitionSphere(numPartitions)
    collar = 0
    while (sum(numRegionsCollar[0:collar+1])) < secNumber: 
        collar += 1
    phi1 = phis[collar]
    if collar+1 < len(numRegionsCollar):
        phi2 = phis[collar+1]
    else:
        phi2 = np.pi
    # Find thetas
    prevStates = sum(numRegionsCollar[0:collar])
    if ((prevStates == 0) or (prevStates == numPartitions -1)):
        theta1 = 0
        theta2 = 2*np.pi
    else:
        thetasCollar = thetas[collar-1]
        statesInCollar = secNumber - prevStates
        theta1 = thetasCollar[statesInCollar-1]
        if statesInCollar == len(thetasCollar):
            theta2 = 2*np.pi
        else:
            theta2 = thetasCollar[statesInCollar]
    phiInterval = [phi1,phi2]
    thetaInterval = [theta1, theta2]
    return phiInterval, thetaInterval
            
