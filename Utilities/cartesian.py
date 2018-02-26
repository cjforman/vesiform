'''
Created on 14 Dec 2017

@author: chris forman
'''
import numpy as np
import coordSystems as coords

# ***** cartesian vector geometry functions ******
def AlignTwoBlocksOfVectors(BlockOne, BlockTwo):
    ''' Repositions and orients BlockTWo such that the tangent between first 
        two elements of block Two is parallel with the tangent between the last 
        two elements of block one. Then translates block two so that its first element 
        is placed at the penultimate position of block one. If you use a staple 
        (a list of 3 vectors containing a bond angle) and two calls to this function 
        you can join two separate blocks of vectors so they form a single block along 
        a chain with the spacing between the end points defined by the spacing in the 
        staple.'''
    director = BlockOne[-1] - BlockOne[-2]
    AlignedBlockTwo = alignBlockWithDirector(director, BlockTwo, BlockTwo[1]-BlockTwo[0])
    return translateBlock(BlockOne[-2], AlignedBlockTwo)

def translateBlock(basePosition, listOfVecs):
    return [ v - listOfVecs[0] + basePosition for v in listOfVecs]

def translateBlockRelRef(basePosition, listOfVecs, refVec):
    return [ v - refVec + basePosition for v in listOfVecs]

def alignBlockWithDirector(director, listOfVecs, testVec):
    ''' Figure out the rotation that rotates between testVec and director when they are both at the origin,
     and then apply that rotation to the list of Vecs.'''
    
    # normalise the testVector
    tNorm = np.linalg.norm(testVec)
    tHat = testVec/tNorm
    
    # normalise director just in case it isn't already
    dNorm = np.linalg.norm(director)
    dHat = director/dNorm
    
    # get angle between director and test vector
    rotAngle = np.arccos(np.dot(tHat, dHat))
    
    # check for co-linearity vector between test and director (if they are colinear, the job is already done...!).
    if abs(rotAngle) > 0.00001:
        nHat = np.cross(tHat, dHat)
        return [ rotPAboutAxis(v, nHat, rotAngle) for v in listOfVecs ]
    else:
        return listOfVecs

def getCentreOfMass(listOfVectors):
    sumVal=np.array([0,0,0])
    for v in listOfVectors:
        sumVal= sumVal + v
    com = sumVal/float(len(listOfVectors))
    return np.around(com, 12) # round the floats to 12th DP. cleans up machine precision noise 

def getPrincipalAxis(listOfVectors):

    iTensor = np.array([[0,0,0], [0,0,0], [0,0,0]])

    for v in listOfVectors:
        iTensor[0][0] += v[1]**2 + v[2]**2
        iTensor[1][1] += v[0]**2 + v[2]**2
        iTensor[2][2] += v[0]**2 + v[1]**2        
        iTensor[0][1] += v[0]*v[1]
        iTensor[0][2] += v[0]*v[2]
        iTensor[1][2] += v[1]*v[2]
    iTensor[1][0] = iTensor[0][1]
    iTensor[2][0] = iTensor[0][2]
    iTensor[2][1] = iTensor[1][2]

    # get the eigen values and vectors of the inertia tensor
    eig = np.linalg.eig(iTensor)
    
    # Take the eigenvector which corresponds to the smallest eigenvalue
    return  eig[1][np.argmax(eig[0]),:]

def rotPAboutAxisAtPoint(p, p0, n, angle):
    # rotates the vector P about an axis n through the point p0 by an angle.
    if not coords.isEqual(p, p0, 1e-10):
        r = p - p0
        nHat = n/np.linalg.norm(n) # ensure we are using a normalised vector
        try:
            r_Rot = r*np.cos(angle) + np.cross(nHat, r)*np.sin(angle) + nHat*np.dot(nHat, r)*(1- np.cos(angle))
            outVec = r_Rot + p0
        except ValueError:
            print "valueError"
            outVec = p
    else:
        outVec = p
    return outVec

def rotPAboutAxis(p, n, angle):
    nHat = n/np.linalg.norm(n) # ensure we are using a normalised vector
    return p*np.cos(angle) + np.cross(nHat, p)*np.sin(angle) + nHat*np.dot(nHat, p)*(1- np.cos(angle))

def vectorBetweenTwoPoints(p1, p2 = np.array([0.0, 0.0, 0.0]), direction='TwoToOne'):
    # Returns a unit vector pointing from p2 to p1, or point1 to point2 if 
    # direction is 'OneToTwo'. If point2 is not specified just returns a unit vector 
    # from origin to point1. (i.e. normalises point 1).
    director = p1 - p2
    if direction=='OneToTwo':
        director = p2 - p1
    d = np.linalg.norm(director)
    return director/d
