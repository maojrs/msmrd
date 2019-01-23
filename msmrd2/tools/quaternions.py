''' Library of functions to operate with quaternions '''

import numpy as np

def rotateVec(vect, rotQuat):
    ''' 
    Rotate vector (np array) or list of vectors (list of np arrays) using a quaternion (4d-nparray).
    '''
    r0 = rotQuat
    r0m = -1.0*rotQuat
    r0m[0] = rotQuat[0]
    v1 = np.insert(vect,0,0)
    v1 = multiply(v1,r0m)
    v1 = multiply(r0,v1)
    rotVect = v1[1:]/np.linalg.norm(v1[1:])
    return rotVect


def phiTransform(dphi,rotQuat):
    ''' Calculate Phi function (4x3 matrix). '''
    s = rotQuat[0]
    p1 = rotQuat[1]
    p2 = rotQuat[2]
    p3 = rotQuat[3]
    phiMatrix = 0.5*np.array([[-p1, -p2, -p3],
                              [s, -p3, -p2],
                              [-p3, s, -p1],
                              [-p2, -p1, s]])
    dquat = np.dot(phiMatrix,dphi)
    return dquat

def multiply(q1,q2):
    ''' Multiply two quaternions. '''
    s1 = q1[0]
    p1 = q1[1:]
    s2 = q2[0]
    p2 = q2[1:]
    s = s1*s2 - np.dot(p1,p2)
    p = s1*p2 + s2*p1 + np.cross(p1,p2)
    qout = np.insert(p,0,s)
    return qout

def conjugate(q):
    '''Return the conjugate of a quaternion'''
    s = q[0]
    p = -q[1:]
    qout = np.insert(p,0,s)
    return qout

def relativeOrientation(q1,q2):
    '''Return relative quaternion between q1 and q2, measured from q1'''
    q1conj = conjugate(q1)
    qout = multiply(q2,q1conj)
    return qout

def angle2quat(phi):
    ''' Transform from angle-axis representation (w*dt) to quaternion representation. '''
    phinorm = np.linalg.norm(phi)
    phiunit = phi/phinorm
    s = np.cos(0.5*phinorm)
    p = np.sin(0.5*phinorm)*phiunit
    qout = np.insert(p,0,s)
    return qout

def quat2angle(q):
    ''' Transform from quaternion representation to angle-axis representation. '''
    phi = q[1:]
    qnorm = np.linalg.norm(phi)
    phi = phi/qnorm
    theta = 2*np.arctan2(qnorm,q[0])
    phi = phi*theta
    return phi

def renormalize(q):
    ''' Renormalize quaternion. '''
    qnorm = np.linalg.norm(q)
    newq = q/qnorm
    return newq

