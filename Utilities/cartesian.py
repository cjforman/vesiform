'''
Created on 14 Dec 2017

@author: chris forman
'''
import numpy as np
import Utilities.coordSystems as coords
from numpy import inf

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

    COM = getCentreOfMass(listOfVectors)

    for u in listOfVectors:
        # make sure that we calculate the PrincipleComponents in the COM frame
        v = u - COM
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
    
    # Take the eigenvector which corresponds to the largest eigenvalue
     
    # makes zeros zero.
    if np.abs(eig[0][0]) < 1e-8:
        eig[0][0] = 0.0
    if np.abs(eig[0][1]) < 1e-8:
        eig[0][1] = 0.0
    if np.abs(eig[0][2]) < 1e-8:
        eig[0][2] = 0.0        

    
    # if there are three distinct eigenvalues then take evector corresponding to smallest    
    if len(np.unique(eig[0]))==3:
        pVec = eig[1][np.argmin(eig[0]), :]

    #if two of the eigenvalues are the same then pick the unique one.
    # symmetric object so pick the eigen value that is orthogonal to the symmetry.        
    if len(np.unique(eig[0]))==2:
        if eig[0][0]==eig[0][1]:
            pVec = eig[1][2, :]
        if eig[0][0]==eig[0][2]:
            pVec = eig[1][1, :]
        if eig[0][1]==eig[0][2]:
            pVec = eig[1][0, :]

    # completely symmetric so any axis is fine.            
    if len(np.unique(eig[0]))==1:       
        pVec = eig[1][np.argmin(eig[0]), :]
        
    return  pVec 

def rotPAboutAxisAtPoint(p, p0, n, angle):
    # rotates the vector P about an axis n through the point p0 by an angle.
    if not coords.isEqual(p, p0, 1e-10):
        r = p - p0
        nHat = n/np.linalg.norm(n) # ensure we are using a normalised vector
        try:
            r_Rot = r*np.cos(angle) + np.cross(nHat, r)*np.sin(angle) + nHat*np.dot(nHat, r)*(1- np.cos(angle))
            outVec = r_Rot + p0
        except ValueError:
            print("valueError")
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

def clamp(val, valMin, valMax):
    retVal = val
    if val < valMin:
        retVal = valMin
    if val> valMax:
        retVal = valMax
    return retVal

def closestApproachTwoLineSegmentsSquared(p1, q1, p2, q2, returnVec=False):
    # finds the square of the distance between the closest points of 
    # approach on two line segments between p1 and q1 and p2 and q2.
    # Works by computing s and t - the displacement parameters along the 
    # line's equations:
    # r(s) = p1 + s * (q1 - p1)
    # r(t) = p2 + t * (q2 - p2)  
    
    d1 = q1 - p1 # line directors
    d2 = q2 - p2 # line directors
    r = p1 - p2 # vector between lines
    a = np.dot(d1, d1) # squared length of segment 1
    e = np.dot(d2, d2) # square length of segment 2
    f = np.dot(d2, r)
    epsilon = 1e-10
    
    # check if both segments are of zero length (degenerate into points)
    if (a <= epsilon) and (e <= epsilon):
        s = t = 0.0
 
    if (a <= epsilon):
        # first segment is a point
        s = 0.0
        t = f/e # s=0 => t = (b*s + f)/e = f/e
        t = clamp(t, 0, 1)
    else:
        c = np.dot(d1, r)
        if (e <= epsilon):
            # second segment degenerates into a point
            t = 0.0
            s = clamp(-c/a, 0.0, 1.0) # t = 0 => s = (b*t - c)/a = -c/a
        else:
            # general non-degenerate case starts here
            b = np.dot(d1, d2)
            denom = a*e - b*b  # always non-negative
    
            # if segments not parallel compute closest point on L1 to L2 and 
            # clamp to segment S1. Else pick arbitrary s (here 0).
            if (denom > epsilon):
                s = clamp((b * f - c * e)/denom, 0.0, 1.0)
            else:
                s = 0.0
    
            # compute point on L2 cloest to s1(s) using
            # t = dot( (p1 + D1 * s) - P2, D2 ) / dot(d2, d2) = (b*s + f)/e
            t = (b * s + f)/e
            
            # if t in 0, 1 we are done. else clamp t, recompute s for new
            # value of t using s = dot((P2 + d2*t) - p1, d1) / dot(d1, d1) = 
            # (t *b -c)/a and clamp s to [0,1]
            if (t < 0.0):
                t = 0.0
                s = clamp( - c /a, 0.0, 1.0)
            elif (t > 1.0):
                t = 1.0
                s = clamp( (b-c) /a, 0, 1.0)
                
    c1 = p1 + d1 * s
    c2 = p2 + d2 * t
    
    # optionally output the director between closest points of approach as well as the distance squared depending on the flag
    director = c1 - c2
    outVals = np.dot(director, director)
    if returnVec:
        outVals = (outVals, c1 - c2) 
     
    return outVals

def closestApproachPointToLineSegmentSquared(p1, p2, q2, returnVec=False):
        return closestApproachTwoLineSegmentsSquared(p1, p1, p2, q2, returnVec=returnVec)

def distanceOfLineToAPlane(linePoint, lineDirector, planePoint, planeDirector):
    # computes the line parameter on the line where the line intersects a plane. 
    # Returns np.inf if the line is parallel or embedded within the plane. 
    # the sign of the displacement is determined by the sense of direction of the lineDirector.
    
    # ensure the directors are normalised
    lineDirector = lineDirector/np.linalg.norm(lineDirector)
    planeDirector = planeDirector/np.linalg.norm(planeDirector)    
    
    nDotD = np.dot(planeDirector, lineDirector)
    # assume there is no solution to start with
    outValue = np.inf
    if np.abs(nDotD) > 1e-6:
        outValue = (np.dot(planeDirector, planePoint) - np.dot(planeDirector, linePoint))/nDotD
    return outValue
    
if __name__ == "__main__":
    
#    linePoint = np.array([0.0, 0.0, 0.0])
#    lineDirector = np.array([0.0, 1.0, 01.0])
#    planePoint = np.array([1.0, 1.0, 0.0])
#    planeDirector = np.array([1.0, 1.0, 1.0])
    
#    d = distanceOfLineToAPlane(linePoint, lineDirector, planePoint, planeDirector)
    
#    print(d)
    
    
    
    
# listOfVectors = [ np.array([1.0, 1.0, 0.0]), np.array([-1.0, -1.0, 0.0]), np.array([-1.0, 1.0, 0.0]), np.array([1.0, -1.0, 0.0])]
# pAxis = getPrincipalAxis(listOfVectors)
# print(pAxis)
 
 
 p1 = np.array([-0.0, 0.0, 0.0])
 q1 = np.array([ -0.0, 0.0, 0.0])
 p2 = np.array([ -5.0, -0.0, 2.0])
 q2 = np.array([ -0.0, -5.0, 2.0])
 
 d = closestApproachTwoLineSegmentsSquared(p1, q1, p2, q2, returnVec=True)
 print(d)
  
