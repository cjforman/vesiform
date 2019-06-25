'''
Created on 14 Dec 2017

@author: chris forman
'''
import numpy as np
import random as rnd
import Utilities.cartesian as cart
import Utilities.fileIO as FIO

# ***** Co-ordinate system specific functions ****** 

def sphericalPolar2XYZ(pos):
    # takes [r, theta, phi] numpy array. Computes the x,y,z coords on unit sphere and 
    # scales it to radius r thus returning x, y, z position of spherical polar input
    unitSpherePos = polarToUnitSphereXYZ(pos[1], pos[2]) 
    return pos[0] * unitSpherePos   

def bondAngleDihedral2TNB(pos): # given r, beta, alpha in terms of bond angle (zero to pi) and dihedral(-pi to pi)
    # return a cartesian vector in the TNB convention.  X->B, y-> N, and Z->T.
    cartVect = sphericalPolar2XYZ([pos[0], pos[1] - np.pi/2, pos[2]])
    return np.array([cartVect[2], cartVect[1], cartVect[0]]) 

def polarToUnitSphereXYZ(theta, phi):
    # take theta and phi and computes the x, y and z position of the point on unit sphere.
    # this (and inverse XYZ2SphericalPolar) is the only place where this transformation is defined.
    # theta is from -pi/2 (south pole) to pi/2 (north pole), measured from xy plane.
    # phi is from -pi to pi. zero at +ve x axis.  
    return np.array([np.cos(phi) * np.cos(theta), np.sin(phi) * np.cos(theta), np.sin(theta)]) 

def XYZ2SphericalPolar(pos): # x, y, z to r, theta (el), phi (az)
    r = np.linalg.norm(pos)
    return np.array([r, np.arcsin(pos[2]/r), np.arctan2(pos[1], pos[0])])

def cyl2XYZ(pos): # rho, phi (az), z to x, y, z
    return np.array([pos[0] * np.cos(pos[1]), pos[0] * np.sin(pos[1]), pos[2]])

def XYZ2Cyl(pos):# x, y, z to rho, phi, z
    # assume cylindrical axis is oriented along z axis
    zVal = np.dot(pos, np.array([0.0, 0.0, 1.0]))
    
    # get the vector in the radial direction
    rhoVec = pos - zVal * np.array([0.0, 0.0, 1.0])
    
    # establish the rho coordinate as magnitude of radial vector
    rho = np.linalg.norm(rhoVec)
    
    if rho > 0.0000001:
        # compute the unit vector in the radial direction
        rhoHat = rhoVec/rho
    else:
        # if rho is zero then phi is undefined so just pick a direction for radial vector.
        rhoHat= np.array([1.0, 0.0, 0.0])
        
    # phi is azimuthal angle between rhoHat and x-axis, 
    phi = np.arccos( np.dot( rhoHat, np.array( [1.0, 0.0 , 0.0] ) ) )
    
    # check sign - positive phi is in positive y half-space and vice versa
    if np.dot(np.array([0.0, 1.0, 0.0]), rhoHat) < 0.0 :
        phi = -phi
    
    return np.array([rho, phi, zVal])

def azToRhoCyl(az, r_x, r_y): 
    # returns the distance from origin to an ellipse at a particular azimuth (in radians), 
    # when the ellipse has semi axis r_x on x axis and r_y on y_axis.zero az is along x_axis.
    a = np.sin(az)/r_x
    b = np.cos(az)/r_y
    p = a**2 + b**2
    return np.power(1/p, 0.5)
    
def ellipsoidalPolarUVW(theta, phi, r_x, r_y, r_z):
    # returns the radial vector along theta and phi in XYZ and also returns the distance along that 
    # vector from the origin to the surface of an ellipsoid with semi major axes, r_x. r_y, and r_z
    # with each semi-major axis aligned with the given xyz axis. 
    
    # takes theta and phi and computes XYZ position on unit sphere.
    unitSpherePos = polarToUnitSphereXYZ(theta, phi)

    # calculate the distance to the surface of the ellipsoid at the given angle theta and phi.
    r =np.power(   ( unitSpherePos[0] * unitSpherePos[0])/(r_x * r_x) 
                 + ( unitSpherePos[1] * unitSpherePos[1])/(r_y * r_y)
                 + ( unitSpherePos[2] * unitSpherePos[2])/(r_z * r_z), -0.5)         
    
    # return the r and u, v, ws
    return [r, unitSpherePos] 

def ellipsoidalPolarToXYZ(theta, phi, r_x, r_y, r_z):
    # computes the XYZ position of the point on the ellispoids surface
    # which is pointed to by a vector from the origin at angle theta, and phi 

    # takes theta and phi and computes XYZ position on unit sphere.
    unitSpherePos = polarToUnitSphereXYZ(theta, phi)

    # calculate the distance to the surface of the ellipsoid at the given angle theta and phi.
    r =np.power(   ( unitSpherePos[0] * unitSpherePos[0])/(r_x * r_x) 
                 + ( unitSpherePos[1] * unitSpherePos[1])/(r_y * r_y)
                 + ( unitSpherePos[2] * unitSpherePos[2])/(r_z * r_z), -0.5)         
    
    # convert the r, theta, phi positions to XYZ format
    return sphericalPolar2XYZ(np.array([r, theta, phi])) 

def TNB2XYZ(TNBframe, posTNB):
    # given a TNB frame defined in XYZ coords, and a TNB vector, posTNB, defined in the TNB frame, 
    # return the XYZ coords of posTNB 
    return np.inner(np.transpose(TNBframe), posTNB)

def isColinear(p1, p2, epsilon):
    # checks to see if two vectors are colinear (cross product is zero). 
    # Zero vector defined as having a magnitude less than epsilon. 
    
    # assume the vectors are colinear
    colinear = True 
    
    n = np.cross(p1, p2)  
    nMag =  np.linalg.norm(n)
    
    if ( abs(nMag - 0.0) > epsilon):
        colinear = False
        
    return colinear 

def isEqual(p1, p2, epsilon):
    # checks to see if two vectors are equal. 
    # Zero vector defined as having a magnitude less than epsilon. 
    equal = True # assume vectors are equal
    n = p1 - p2
    nMag =  np.linalg.norm(n)
    
    if ( abs(nMag - 0.0) > epsilon):
        equal = False
        
    return equal


def isZero(p1, epsilon):
    # checks to see if a vectors is zero. 
    # Zero vector defined as having a magnitude less than epsilon. 
    zero = True # assume vector is zero
    nMag =  np.linalg.norm(p1)
    
    if ( abs(nMag - 0.0) > epsilon):
        zero = False
        
    return zero

def axisFromHelix(xyzList):
    ''' returns normalised vector in direction of mean of cross product of adjacent bivectors 
    for a sequence of points in space. If the points are in a helix, this returns the axis of the helix. 
    to understand this better:  The tangent vectors are the vectors between successive points. 
    Any two adjacent tangents form a plane.  The bivector bisects the internal angle between the tangent
    in the plane of the two tangents.  The cross product of two adjacent bivectors should yield a 
    vector that is approximately aligned with the axis of a helix if the original points are on a helix. 
    We take the mean of all those vectors as defined an axis for a set of points. Returns the unit 
    vector in that direction. '''
    
    if len(xyzList)==0 or len(xyzList)==1:
        return np.array([0.0, 0.0, 1.0])
    
    if len(xyzList)==2:
        return (xyzList[1] - xyzList[0])/np.linalg.norm((xyzList[1] - xyzList[0])) 
    
    # compute tangent unit vectors of backbone
    TangentVectors = [ x- y for x, y in zip(xyzList[1:], xyzList[0:-1])]
    TangentVectorsHat = [T/np.linalg.norm(T) for T in TangentVectors]
    
    # compute adjacent bivectors of tangents (points to central axis if polymer is a helix)
    bivectors = [ x-y for x, y in zip(TangentVectorsHat[1:], TangentVectorsHat[0:-1]) ]
     
    # take cross product of adjacent bivectors (if they are truly radial would be vector parallel to axis
    NorVectors = [ np.cross(y, x) for x, y in zip(bivectors[0:-1], bivectors[1:]) ]

    # take mean of all these normal vectors as an approximation to the central axis of the backbone (if there is one)        
    meanNorVectors = sum(NorVectors)/len(NorVectors)
     
    axisVec = -1 * meanNorVectors/(np.linalg.norm(meanNorVectors))
    return np.around(axisVec, 12) # round the floats to 12th DP. Cleans up machine precision noise on zeros.


def transformFromBlockFrameToLabFrame(labDirector, labRefPoint, labRotation, blockDirector, blockRefPoint, blockXYZVals):
    # Returns the input xyzVals rotated to a specific orientation and position.
    # Three distinct operations are performed on the blockXYZCoords
     
    # 1) Rotate block coords about the blockRefPoint so that blockFrame director and labFrame directors 
    #    placed at the blockRefPoint coincide.
    # 2) Rotate the resulting block coords about the labDirector axis through the blockRefPoint by labRotRadians  
    # 3) Translate the block so the blockRefPoint is at the labRefPoint  
    # These are the exact reverse of the transformFromLabFrameToBlockFrame
    
    labDirectorHat = labDirector/np.linalg.norm(labDirector)
    blockDirectorHat = blockDirector/np.linalg.norm(blockDirector)
    
    # Compute the rotation axis about which to rotate the backbone Director to match the given director
    rotAxis = np.cross(blockDirectorHat, labDirectorHat)
                       
    # calculate the angle to rotate the director to match the given director                            
    rotAngle = np.arccos(np.dot(labDirectorHat, blockDirectorHat))
    
    # rotate each point in the building block by the calculated angle about the 
    # calculated axis through the reference point. This way the reference point
    # won't move with the rotation. 
    if rotAngle > 1e-5:
        xyzVals = [ cart.rotPAboutAxisAtPoint(pos, blockRefPoint, rotAxis, rotAngle) for pos in blockXYZVals]
    else:
        xyzVals = blockXYZVals
    
    # tranform the building block to the lab refPoint
    xyzVals = [pos - blockRefPoint + labRefPoint for pos in xyzVals]
    
    # now rotate by "labRotRadians" about the labDirector through the lab Reference point.
    return [ cart.rotPAboutAxisAtPoint(pos, labRefPoint, labDirectorHat, labRotation) for pos in xyzVals]

def transformFromLabFrameToBlockFrame(labDirector, labRefPoint, labRotation, blockDirector, blockRefPoint, labXYZVals):
    # Returns the input xyzVals rotated to a specific orientation and position.
    # Three distinct operations are performed on the labXYZCoords in this order: 

    # 1) Rotate the lab coords about the labDirector axis through the labRefPoint by -labRotRadians  
    # 2) Translate the lab coords so they are relative to block Ref point rather than lab refpoint. 
    # 3) Rotate the resulting coords about the blockRefPoint so that blockFrame director and 
    #    labFrame directors placed at the blockRefPoint coincide.
    # These are the exact reverse of the transformFromBlockFrameToLabFrame function 
  

    # normalise lab and block directors        
    labDirectorHat = labDirector/np.linalg.norm(labDirector)
    blockDirectorHat = blockDirector/np.linalg.norm(blockDirector)
    
    # Rotate by "-labRotation" about the labDirector through the lab Reference point.
    xyzVals = [ cart.rotPAboutAxisAtPoint(pos, labRefPoint, labDirectorHat, -labRotation) for pos in labXYZVals]
    
    # The input coords are relative to the labRefPoint, transform them so they are relative to the blockRefPoint.
    xyzVals= [pos - labRefPoint + blockRefPoint for pos in xyzVals]

    # Compute the rotation axis about which to rotate the labDirector to match the blockDirector
    rotAxis = np.cross(blockDirectorHat, labDirectorHat)
                       
    # calculate the angle to rotate the labDirector to match the blockDirector                            
    rotAngle = np.arccos(np.dot(labDirectorHat, blockDirectorHat))
    
    # rotate each point in the xyzCoords by the calculated angle about the 
    # calculated axis through the block Reference point. 
    return [ cart.rotPAboutAxisAtPoint(pos, blockRefPoint, rotAxis, -rotAngle) for pos in xyzVals]

def convertTriplesToEllipsoids(blockNames, xyzVals, aAxis, cAxis):
    eAtomNames = [ grainName for grainName in blockNames if grainName in ['P', 'B']]  
    ePositions = [ pos for pos, grainName in zip(xyzVals, blockNames) if grainName in ['P', 'B']]
    eSizes = [ [2*aAxis, np.linalg.norm(xyzVals[i - 1] - xyzVals[i + 1]), 2*cAxis] for i, grainName in enumerate(blockNames[0:-1]) if grainName in ['P', 'B']] 
    eRs = [ XYZPointsToOrientationMat(xyzVals[i - 1], xyzVals[i + 1]) for i, grainName in enumerate(blockNames[0:-1]) if grainName in ['P', 'B']]
    eRotVecs = [ eR[1] for eR in eRs]
    return eAtomNames, ePositions, eSizes, eRs, eRotVecs

def XYZPointsToOrientationMat(pos1, pos2):
        # Computes the unit vector between two points, and two further orthogonal vectors.
        # The first orthogonal vector is a vector orthogonal to the unitVector in the XYPlane.
        # This is found by projecting the first unit into the xy plane, and then finding a 
        # second vector in the plane, and doing a gram-Schmit orthogonalisation. 
        # The second orthogonal vector is the cross product of the first two unit vectors.
        bVec = (pos1-pos2)
        bVecHat = bVec/np.linalg.norm(bVec)
        xyPlaneVec = np.array([bVecHat[0], bVecHat[1], 0.0])
        xyPlaneVecHat = xyPlaneVec/np.linalg.norm(xyPlaneVec)
        if np.abs(xyPlaneVecHat[1]) < 1e-10:
            otherXYPlaneVec = xyPlaneVecHat + np.array([1.0, 1.0, 0,0])
        else:
            otherXYPlaneVec = xyPlaneVecHat + np.array([1.0, 0.0, 0.0])
        otherXYPlaneVecOrth = otherXYPlaneVec - np.dot(otherXYPlaneVec, xyPlaneVecHat) * xyPlaneVec
        otherXYPlaneVecOrthHat = otherXYPlaneVecOrth/np.linalg.norm(otherXYPlaneVecOrth) 
        
        finalVec = np.cross(otherXYPlaneVecOrthHat, bVecHat)
        return np.array([otherXYPlaneVecOrthHat, bVecHat, finalVec])

    
def constructTNBFrame(p1, p2, p3):
    # given three arbitrary points defined in the lab XYZ frame, this function constructs 
    # three orthogonal vectors (T, N and B). T points from the second seed point in the list to the third seed point.
    # The N axis is normal to the plane formed by the three points, and is therefore orthogonal to T.
    # B is the cross product of T and N. The three unit vectors in these directions are returned as a numpy array of numpy arrays (a matrix)
    # returns None if p1, p2 and p3 are colinear
    
    # assume failure
    outVecs=None
    
    T = p3 - p2  # equivalent of b2 in bond angle and dihedral functions
    tHat = T/np.linalg.norm(T)

    U = p2 - p1 # equivalent of b1 in bond angle and dihedral functions
    uHat = U/np.linalg.norm(U)
    
    n = np.cross(uHat, tHat)  # equiv of n1 = b1 x b2 in dihedrals
    nMag =  np.linalg.norm(n)
    
    if ( abs(nMag - 0.0) > 1e-10):  # check for colinearity (cross product of T and U is zero)
        # now we're not colinear can divide by norm(n)
        nHat= n/nMag
        
        b = np.cross(nHat, tHat)  # equiv of m1 = n1Hat x b2Hat in dihedral reconstruction 
        bHat = b/np.linalg.norm(b) # just in case
        outVecs = np.array([tHat, nHat, bHat])
    
    return outVecs

def bondAngle(p1, p2, p3):
    # returns the size of the interior angle at P2.
    # This is the angle between the vectors from P2 to the other points.
    # Returns a value in the range 0 to pi
    U = p1 - p2
    V = p3 - p2
    UHat = U/np.linalg.norm(U)
    VHat = V/np.linalg.norm(V)
    return np.arccos(np.clip(np.dot(UHat, VHat), -1.0, 1.0))
    # construct the TNB frame
    #TNB = constructTNBFrame(p1, p2, p3)

    #if TNB.all()==None:
        # co linear points so bond angle is either 0 or 180.
        # if line segment 1 to 3 is great than 1 to 2 then we have 180 else 0
    #    lineSeg1To3 = np.linalg.norm(p3-p1)
    #    lineSeg1To2 = np.linalg.norm(U)
    #    bondAngle = 0.0
    #    if lineSeg1To3>lineSeg1To2:
    #        bondAngle = np.pi
    #else:
        # the bond angle is the angle between U and T in the T, B plane
        # U's coords in the T, B plane are U.THat = p and U.BHat = q
        # Atan2(q,p) gives the internal angle between U and T (b1 and b2) 
        # thus beta (bond angle is pi minus this.
    #    q = np.dot(U, TNB[2])
    #    p = np.dot(U, TNB[0])
    #    bondAngle = np.pi - abs(np.arctan2(q, p)) # atan2(u.bHat, u.tHat)
    
    #return bondAngle

def Dihedral(p1, p2, p3, p4):
    # returns dihedral angle of four points in range -pi to pi
    
    # not the most efficient algorithm but it guarantees to use the same definition as the TNB system.
    
    # generate the TNB frame 
    # b1, b2, b3 are the vectors from p1 to p2, p2 to p3 and p3 to p4 
    # THat = b2Hat, 
    # NHat = n1Hat = b1 x b2
    # BHat = m1Hat = n1Hat x b2Hat = nHat x THat)
    dihedral = None #assume failure
    TNB1 = constructTNBFrame(p1, p2, p3)
    if TNB1.all()!=None:
        TNB2 = constructTNBFrame(p2, p3, p4)
        if TNB2.all()!=None:
            # the dihedral is the angle between N2 (from TNB2) and N1 from TNB1, in the N1, B1 plane.
            dihedral = np.arctan2(-np.dot(TNB1[2], TNB2[1]), np.dot(TNB1[1], TNB2[1])) 
    
    if TNB1.any()==None or TNB2.any()==None:
        dihedral=0.0
    
    return dihedral

def measureAnglesAtConnection(s0, s1, s2, m2, m1, m0):
    # measure d, alpha1, beta1, alpha2, beta2 and alpha3 angles for a six
    # atom chain. Allows characterisation and reconstruction of a connection between two building blocks.
    # We use the same codes for the atoms in the connector that we use in the building blocks to aid 
    # comprehension. This doesn't belong in the building block function because the atoms may come
    # from anywhere. and the alphas and betas are input values into the buildingblock construction.
    d = np.linalg.norm(s2 - m2)
    alpha1 = Dihedral(s0, s1, s2, m2)
    beta1 = bondAngle(s1, s2, m2)
    alpha2 = Dihedral(s1, s2, m2, m1)
    beta2 = bondAngle(s2, m2, m1)
    alpha3 = Dihedral(s2, m2, m1, m0)
    return d, alpha1, beta1, alpha2, beta2, alpha3

def computeDistBetweenM1AndS3(angle, *data):
    # Used for the hairPinResidueToAtoms function.
    # Compute the difference between bondLength and the distance between M1 and S3 for a given 
    # value of dihedral around the S2 to M0 axis
    TNB_S1_S2_S3_M0 = data[0]
    beta = data[1]
    S3 = data[2]
    bondLength = data[3]
    M1 = generateTNBVecXYZ(TNB_S1_S2_S3_M0, beta, angle)
    return np.abs(np.linalg.norm(S3 - M1) - bondLength) 

def computeDistBetweenPoints(angles, *data):
    # Used for fitting angles to a connection to yield a specific distance.
    # Function computes the difference between a given value and the distance between two
    # points. The two points are defined by their alpha and beta in the given TNB 
    # frames when the frames are placed at p1 and p2 respectively. In this way we 
    # can vary the alphas and get a single error for the bondlength.  
    # If that error number is zero then we have solved the constraint problem.
    
    # Can use this function for both the two point or one point final end states 
    # by setting the target distance between the points as zero or bondLength.   
    # alphaA = angles[0]
    # alphaB = angles[1]
    # betaA = angles[2]
    # betaB = angles[3]
    
    # p1 = data[0]
    # p2 = data[1]
    # TNBA = data[2]
    # TNBB = data[3]
    # betaA = data[4]
    # betaB = data[5]
    # bondLength = data[6]
    # distBetweenPoints = data[7]

    # compute the distance - bondlength between the two new points defined by the angles
    # in the given TNB frames
    pointA = data[0] + data[6] * generateTNBVecXYZ(data[2], data[4], angles[0])
    pointB = data[1] + data[6] * generateTNBVecXYZ(data[3], data[5], angles[1])
    
    return np.linalg.norm(pointA - pointB) - data[7]


def cpoa(p1List, p2List, p3List, p4List):
    '''Computes the point of closest approach between the vector between p1 and p2 and the vector between p3 and p4.
    p1List, p2List, p3List, and p4List are lists of numpy arrays of the same length
    
    Returns Lists of: 
    pa - the vector on the line p2-p1 which is at the point of closest approach
    mua - distance from p1 to pa.
    pb - the vector on the line p4-p3 which is at the point of closest approach
    mub - distance from p1 to pb.'''
    p13List=[p1-p3 for p1,p3 in zip(p1List,p3List)]
    p43List=[(p4-p3) / np.linalg.norm(p4-p3) for p4,p3 in zip(p4List,p3List)]
    p21List=[(p2-p1) / np.linalg.norm(p2-p1) for p2,p1 in zip(p2List,p1List)]
    
    d1343List=[np.dot(p13,p43) for p13,p43 in zip(p13List,p43List)]
    d4321List=[np.dot(p43,p21) for p43,p21 in zip(p43List,p21List)]
    d1321List=[np.dot(p13,p21) for p13,p21 in zip(p13List,p21List)]
    d4343List=[np.dot(p43,p43) for p43 in p43List]
    d2121List=[np.dot(p21,p21) for p21 in p21List]

    denomList=[(d2121*d4343) - (d4321 * d4321) for d2121,d4343,d4321 in zip(d2121List,d4343List,d4321List)]
    numerList=[(d1343*d4321) - (d1321 * d4343) for d1343,d4321,d1321,d4343 in zip(d1343List,d4321List,d1321List,d4343List)]

    mubList=[]
    muaList=[]
    pbList=[]
    paList=[]
    for denom,numer,d1343,d4321,d4343,p1,p21,p3,p43 in zip(denomList,numerList,d1343List, d4321List,d4343List,p1List,p21List,p3List,p43List):
        mua=0
        mub=0
        pa=p1
        pb=p3
        if (abs(denom)>1e-14):
            mua=numer/denom
            mub= (d1343+ d4321*(numer/denom))/d4343
            pa=p1+mua*p21
            pb=p3+mub*p43
    
        muaList.append(mua)
        mubList.append(mub)
        paList.append(pa)
        pbList.append(pb)
        
    return [muaList, paList, mubList, pbList]


def CylArea(phiMin, phiMax, radius, length):
    return (phiMax - phiMin) * radius * length

def SphereArea(phiMin, phiMax, thetaMin, thetaMax, radius):
    return radius*(phiMax - phiMin)*(np.sin(thetaMax) - np.sin(thetaMin))

def FrustumRadiusAtGivenZ(testZ, z1, z2, r1, r2):
    # Returns the radius of a circular frustum at a given intermediate z height. The 
    # central axis of the frustum is parallel with the z-axis and it is 
    # defined by two radii at distinct z heights. The function interpolates the radius of the frustum at 
    # an intermediate z height. 
    return (( testZ - z1) * (r2 - r1 ) / ( z2 - z1 ) )  +  r1

def FrustumZAtZeroRadius(z1, z2, r1, r2):
    retZ = None
    DZ = np.abs(z2-z1)
    DR = np.abs(r2-r1)
    
    if DZ<1e-10:
        retZ = z1
    else:
        if (DR > 1e-10):
            # take +ve gradient branch
            m = np.abs(DR/DZ)
            c = r1 - z1 * m            
            retZ = -c/m
    if retZ==None:
        print("Warning: attempt to find zero point of cone that is actually a cylinder.")
    
    return retZ
        
def checkPointInFrustum(pos, Z1, R1, Z2, R2, verbose):
    # check to see if a position is inside a given frustum
    # returns true if it is.
    ZMin = min(Z1, Z2)
    ZMax = max(Z1, Z2)
    
    # principal axis of envelope is aligned with Z axis. 
    # Compute radial distance from z axis and z height
    zComponent = pos[2] 
    r = np.sqrt(pos[0]**2 + pos[1]**2)
        
    if np.abs(zComponent - 0.0) < 1e-10:
        zComponent = 0.0
        
    inBounds = True # assume it's inside the frustum
    if  zComponent > ZMax or zComponent < ZMin:
        inBounds = False
        if verbose==1:
            print("frustum Z Violation")

    if inBounds:     
        if r > FrustumRadiusAtGivenZ(zComponent, Z1, Z2, R1, R2):
            inBounds = False
            if verbose==1:
                print("frustum R Violation") 
    
    return inBounds
    
    
def rotatePointInChain(points, angle):
    # takes a list of three points and rotates point 2 about an axis between points 1 and 3 
    # by the given angle. Conserves the distance between the points and bond angle but changes the bond angles in a chain.
    return cart.rotPAboutAxisAtPoint( points[1], 
                                      points[0], 
                                      points[2] - points[0],
                                      angle)
    
def rotateDihedral(points, angle):
    # takes a list of three points and rotates point 2 about an axis between points 0 and 1
    # with the given angle. preserves the distance between the points and bond angle between the points. 
    return cart.rotPAboutAxisAtPoint( points[2], 
                                      points[1], 
                                      points[1]-points[0],
                                      -angle)      # +ve in rotaAxis is different from typical dihedral def so *-1.     
    
      
def pickRandomPointInFrustumAtGivenZHeightXYZ(z, z1, z2, r1, r2):
    # returns a random point from a given layer inside the frustum defined by z1, z2, r1 and r2.
    
    # select the first point at random inside the unit circle. return azimuthal angle -pi to pi and 0 <= r <=1
    phi, rUnit = pickRandomPointInUnitCircle()
        
    # compute the max value at the given z height, and use this to scale the rUnit radius.
    rMax = FrustumRadiusAtGivenZ(z, z1, z2, r1, r2)
        
    return cyl2XYZ(np.array([rMax * rUnit, phi, z]))

def pickRandomPointInFrustumXYZ(z1, z2, r1, r2):
    # returns a random point inside the frustum define by z1, z2, r1 and r2.
    
    insideFrustum = False
    
    while not insideFrustum: 
    
        # get a point inside a unit cylinder and scale it to Z and the larger of r vals
        cylPos = pickRandomPointInUnitCylinder()
    
        cylPos[0] = cylPos[0] * max(r2, r1)
        
        cylPos[2] = cylPos[2] * (z2 - z1) + z1
        
        radAllow = FrustumRadiusAtGivenZ(cylPos[2], z1, z2, r1, r2)
        
        if cylPos[0] < radAllow:
            insideFrustum = True
    
    return cyl2XYZ(cylPos)

def generateTNBVecXYZ(genTNB, beta, alpha):
    # This function returns a unit vector pointing in the direction 
    # in the TNB frame defined by alpha and beta.
    # The vector is returned in the XYZ coords of the frame in 
    # which the TNB vectors are defined. 
     
    # Alpha is the azimuthal angle about the T vector (first in TNB list). 
    # Alpha = 0 is in the direction of the B vector (third in list)
    # Beta is the elevation angle which is zero in the B-N plane.
    # When the TNB frame is defined by the last three monomers in a polymer chain, 
    # Then alpha is the dihedral angle defined by the three monomers and the newly picked vector.
    # Beta is the angle between the new "bond" and the penultimate "bond".
    posTNB = bondAngleDihedral2TNB(np.array([1.0, beta, alpha]))
    return TNB2XYZ(genTNB, posTNB)

def pickRandomTNBDirectionInAngRangeXYZ(TNB, beta1, beta2, alpha1, alpha2):
    # This function returns a unit vector pointing in a random direction 
    # in the TNB frame that is somewhere within the specified alpha and beta angular ranges.
    # where alpha is dihedral (-pi < alpha < pi) and beta is bond angle (0<beta<pi)
    # The vector is returned in the XYZ coords of the frame in which the TNB vectors are defined. 
     
    # Alpha is the azimuthal angle about the T vector (first in TNB list). 
    # Alpha = 0 is in the direction of the B vector (third in list)
    # alpha +ve is clock wise looking along the T vector.
    # Beta is the elevation angle (bond angle) which is at pi/2 in the BN plane.
    # When the TNB frame is defined by the last three monomers in a polymer chain, 
    # Then alpha is the dihedral angle defined by the three monomers and the newly picked vector.
    # Beta is the angle between the new "bond" and the penultimate "bond".
    # random sphere picker returns a point -pi/2 < theta < pi/2 but we
    # want result between 0 and pi.
    beta, alpha = pickRandomPointOnUnitSphereInAngRange(beta1 - np.pi/2, beta2 - np.pi/2, alpha1, alpha2)
    posTNB = bondAngleDihedral2TNB(np.array([1.0, beta + np.pi/2, alpha]))
    return TNB2XYZ(TNB, posTNB)
    
def pickRandomPointOnUnitSphereInAngRange(theta1, theta2, phi1, phi2):
    ''' Picks a random direction in spherical polars with uniform probablity.  
    Returns theta and phi:
    -pi/2 <= Theta <= pi/2 is the elevation/latitude with -pi/2 at southern pole and pi/2 at northern pole of sphere.
    -pi <= phi <= pi is the azimuth. '''
    # pick a random direction in U, V space:
    v1 = ( np.sin( theta1 ) + 1 ) / 2 
    v2 = (np.sin( theta2 ) + 1 ) / 2 
    u = rnd.uniform( 0.0, 1.0)
    v = rnd.uniform( v1, v2 )

    # compute theta and phi from U and V (This method ensures a uniform selection of
    # angle without biasing poles, due to variable density of states in polar coords yada yada).
    theta = np.arcsin((2 * v) - 1) # theta1 to theta2
    phi = (phi2 - phi1)*u + phi1  # phi1 to phi2

    return theta, phi

def pickRandomPointOnUnitSphere():
    ''' Picks a random direction in spherical polars with uniform probablity.  
    Returns theta and phi:
    -pi/2 <= Theta <= pi/2 is the elevation/latitude with -pi/2 at southern pole and pi/2 at northern pole of sphere.
    -pi <= phi <= pi is the azimuth. '''
    # pick a random direction in U, V space:
    u = rnd.uniform(0.0, 1.0)
    v = rnd.uniform(0.0, 1.0)

    # compute theta and phi from U and V (This method ensures a uniform selection of
    # angle without biasing poles, due to variable density of states in polar coords yada yada).
    theta = np.arcsin((2 * v) - 1) # -pi/2 to pi/2
    phi = 2*np.pi*u - np.pi  # -pi to pi.

    return theta, phi

def pickRandomPointUniformlyOnEllipsoid(a, b, c):
    # picks a random point on the ellipsoid aligned with the x, y and z axes
    # where a is the semi axis along x, b along y and c along z.
    # uses a method by williamson http://iopscience.iop.org/article/10.1088/0031-9155/32/10/009/pdf
    # much understanding gleaned from chen and glotzer: https://arxiv.org/ftp/cond-mat/papers/0701/0701125.pdf
    # and stack exchange conversation: https://math.stackexchange.com/questions/973101/how-to-generate-points-uniformly-distributed-on-the-surface-of-an-ellipsoid
    
    # set up the test variables 
    testVal = 1.0
    gOverGMax = 0.0

    # precompute squares and squares of squares of semi-axes. 
    aSq = float(a) * float(a)
    bSq = float(b) * float(b)
    cSq = float(c) * float(c)
    a4 = aSq * aSq
    b4 = bSq * bSq
    c4 = cSq * cSq
    
    # compute 1/gmax (easy to see from eqn 17 in williamson.
    # g = R^4 * (u^2/a^4 + v^2/b^4 + w^2/c^4)^0.5
    # Max R is obviously at longest semi-axis so that's easy.
    # max val occurs at poles where u, v or w are 1 or -1. But all mutually exclusive.
    # so must find biggest among (1/a^4)^0.5, (1/b^4)^0.5 or (1/c^4)^0.5.
    # i.e. smallest a^2, b^2, or c^2.
    minSq = float(min([aSq, bSq, cSq]))
    max4 = float(max([a4,b4,c4]))
    
    gMaxInv = minSq/max4

    while not (testVal <= gOverGMax): 
        # pick random number between 0 and 1
        testVal = rnd.uniform(0,1) 
               
        # pick a point on our domain (unit sphere)
        theta, phi = pickRandomPointOnUnitSphere()

        # calculate XYZ position of random point on surface (this distorts the 
        # uniformity of sphere picking
        # surfacePointXYZ = ellipsoidalPolarToXYZ(theta, phi, a, b, c)
        r, surfacePointUVW = ellipsoidalPolarUVW(theta, phi, a, b, c)
        
        # calculate g/gmax for that point and loop to test inequality - eqn from williamson
        gOverGMax = gMaxInv * np.power(r,4) * np.power(  np.power( surfacePointUVW[0]/ aSq , 2 )  
                                                       + np.power( surfacePointUVW[1]/ bSq , 2 )
                                                       + np.power( surfacePointUVW[2]/ cSq , 2 ), 0.5 )
    
    return np.array([r, theta, phi]) # return a pos array in spherical polar coords

def pickRandomPointUniformlyOnEllipsoidInAngRange(a, b, c, theta1, theta2, phi1, phi2):
    # picks a random point on the ellipsoid aligned with the x, y and z axes
    # where a is the semi axis along x, b along y and c along z.
    # uses a method by williamson http://iopscience.iop.org/article/10.1088/0031-9155/32/10/009/pdf
    # much understanding gleaned from chen and glotzer: https://arxiv.org/ftp/cond-mat/papers/0701/0701125.pdf
    # and stack exchange conversation: https://math.stackexchange.com/questions/973101/how-to-generate-points-uniformly-distributed-on-the-surface-of-an-ellipsoid
    
    # set up the test variables 
    testVal = 1.0
    gOverGMax = 0.0

    # precompute squares and squares of squares of semi-axes. 
    aSq = float(a) * float(a)
    bSq = float(b) * float(b)
    cSq = float(c) * float(c)
    a4 = aSq * aSq
    b4 = bSq * bSq
    c4 = cSq * cSq
    
    # compute 1/gmax (easy to see from eqn 17 in williamson.
    # g = R^4 * (u^2/a^4 + v^2/b^4 + w^2/c^4)^0.5
    # Max R is obviously at longest semi-axis so that's easy.
    # max val occurs at poles where u, v or w are 1 or -1. But all mutually exclusive.
    # so must find biggest among (1/a^4)^0.5, (1/b^4)^0.5 or (1/c^4)^0.5.
    # i.e. smallest a^2, b^2, or c^2.
    gMaxInv = (min([aSq, bSq, cSq]))/max([a4,b4,c4])

    while not (testVal <= gOverGMax): 
        # pick random number between 0 and 1
        testVal = rnd.uniform(0,1) 
               
        # pick a point on our domain (unit sphere)
        theta, phi = pickRandomPointOnUnitSphereInAngRange(theta1, theta2, phi1, phi2)

        # calculate XYZ position of random point on surface (this distorts the 
        # uniformity of sphere picking
        # surfacePointXYZ = ellipsoidalPolarToXYZ(theta, phi, a, b, c)
        r, surfacePointUVW = ellipsoidalPolarUVW(theta, phi, a, b, c)
        
        # calculate g/gmax for that point and loop to test inequality - eqn from williamson
        gOverGMax = gMaxInv * np.power(r,4) * np.power(  np.power( surfacePointUVW[0]/ aSq , 2 )  
                                                       + np.power( surfacePointUVW[1]/ bSq , 2 )
                                                       + np.power( surfacePointUVW[2]/ cSq , 2 ), 0.5 )
    
    return np.array([r, theta, phi]) # return a pos array in spherical polar coords

def pickRandomPointInUnitCircleInAngRange(phi1, phi2):
    ''' Picks a random point in cylindrical coordinates at z = 0 '''
    # pick a random direction in U, V space:
    r = rnd.uniform(0.0, 1.0)
    phi = pickRandomPointOnUnitCircleInAngRange(phi1, phi2)
    return phi, r

def pickRandomPointOnUnitCircleInAngRange(phi1, phi2):
    phiMin = min([phi1, phi2])
    phiMax = max([phi1, phi2])
    if (phiMax - phiMin) > 1e-5:
        v = rnd.uniform(0.0, 1.0)
        phi = (phiMax - phiMin) *v + phiMin  # phi1 to phi2.
    else:
        print("Phi Range too small")
        phi = None
    return phi    

def pickRandomPointInAnnulus(r1, r2, phi1, phi2):
    ''' Picks a random point in a sector of an annulus defined by concentric circles of radius r1 and r2.
     and phi1 and phi2.'''

    rMin = min([r1, r2])
    rMax = max([r1, r2])
    
    if (rMax-rMin)> 1e-5:    
        # pick a random direction
        phi = pickRandomPointOnUnitCircleInAngRange(phi1, phi2)
        
        # pick a random r in the large circle - with a scaling for that fact that there are more points at larger r.
        # reject any points that are within the smaller circle
        r = rMin - 1
        while r < rMin:
            r = rMax * np.sqrt(rnd.uniform(0, 1))
    else:
        # if the two r's are too close in value then return None
        print("Annulus too thin to pick a random point. r1 = r2")
        r=None
        phi=None

    return r, phi


def pickRandomPointInUnitCircle():
    ''' Picks a random point in cylindrical coordinates at z = 0 '''
    # pick a random direction in U, V space:
    r = np.sqrt(rnd.uniform(0.0, 1.0))
    phi = pickRandomPointOnUnitCircle()
    return phi, r

def pickRandomPointInEllipsoid(rx, ry, rz):
    r = np.sqrt(rnd.uniform(0.0, 1.0))
    phi = pickRandomPointOnUnitCircle() # -pi to pi
    theta = pickRandomPointOnUnitCircle()/2 # -pi/2 to pi/2
    x = rx * r * np.cos(theta)*np.cos(phi)
    y = ry * r * np.cos(theta)*np.sin(phi)
    z = rz * r * np.sin(theta)
    return np.array([x, y, z])

def pickRandomPointInSphericalRange(r1, r2, theta1, theta2, phi1, phi2):
    rMin = min(r1, r2)
    rMax = max(r1, r2)
    r = np.sqrt(rnd.uniform(rMin/rMax, 1.0))
    
    phiMin = min(phi1, phi2)
    phiMax = max(phi1, phi2)
    phi = pickRandomPointOnUnitCircle()
    while phi<phiMin or phi > phiMax:
        phi = pickRandomPointOnUnitCircle() # -pi to pi
    thetaMin = min(theta1, theta2)
    thetaMax = max(theta1, theta2)
    theta = pickRandomPointOnUnitCircle()/2.0 # -pi/2 to pi/2
    while theta < thetaMin or theta > thetaMax:
        theta = pickRandomPointOnUnitCircle()/2.0 # -pi/2 to pi/2

    x = r * np.cos(theta)*np.cos(phi)
    y = r * np.cos(theta)*np.sin(phi)
    z = r * np.sin(theta)
    return np.array([x, y, z])


def pickRandomPointInEllipsoidRange(rx, ry, rz, theta1, theta2, phi1, phi2):
    r = np.sqrt(rnd.uniform(0.0, 1.0))
    
    phiMin = min(phi1, phi2)
    phiMax = max(phi1, phi2)
    phi = pickRandomPointOnUnitCircle()
    while phi<phiMin or phi > phiMax:
        phi = pickRandomPointOnUnitCircle() # -pi to pi
    thetaMin = min(theta1, theta2)
    thetaMax = max(theta1, theta2)
    theta = pickRandomPointOnUnitCircle()/2.0 # -pi/2 to pi/2
    while theta < thetaMin or theta > thetaMax:
        theta = pickRandomPointOnUnitCircle()/2.0 # -pi/2 to pi/2

    x = rx * r * np.cos(theta)*np.cos(phi)
    y = ry * r * np.cos(theta)*np.sin(phi)
    z = rz * r * np.sin(theta)
    return np.array([x, y, z])


def pickRandomPointOnUnitCircle():
    v = rnd.uniform(0.0, 1.0)
    phi = 2*np.pi*v - np.pi  # -pi to pi.
    return phi

def pickRandomPointOnUnitCylinderInAngRange(phi1, phi2):
    ''' Picks a random point in cylindrical coordinates '''
    # pick a random direction in U, V space:
    z = rnd.uniform(0.0, 1.0)
    phi = pickRandomPointOnUnitCircleInAngRange(phi1, phi2)
    return phi, z

def pickRandomPointOnUnitCylinder():
    ''' Picks a random point in cylindrical coordinates '''
    # pick a random direction in U, V space:
    z = rnd.uniform(0.0, 1.0)
    phi = pickRandomPointOnUnitCircle()
    return phi, z

def pickRandomPointInUnitCylinder():
    ''' Picks a random point in cylindrical coordinates '''
    # pick a random point on z axis
    z = rnd.uniform(0.0, 1.0)
    # pick a random point in the unit circle
    phi, r = pickRandomPointInUnitCircle()
    return np.array([r, phi, z])

def pickRandomPointOnUnitSquare():
    ''' Picks a random point in a 2D square '''
    # pick a random direction in U, V space:
    x = rnd.uniform(0.0, 1.0)
    y = rnd.uniform(0.0, 1.0)
    return x, y

if __name__== '__main__':

    rx = 10
    ry = 20
    pos=[]
    for az in np.linspace(-np.pi, np.pi, num=360): 
        rho = azToRhoCyl(az, rx, ry)
        pos.append(cyl2XYZ(np.array([rho, az, 0.0])))
                   
    FIO.saveXYZ(pos, 'C', "ellispoid.xyz")
    
    print("LabToBlockTransformation")
    labRefPoint = np.array([10.0, 10.0, 10.0])
    labDirector = np.array([1.0, 1.0, 1.0])
    labRotation = 20 * np.pi/180.0
    blockDirector = np.array([0.0, 0.0, 1.0])
    blockRefPoint = np.array([0.0, 0.0, 0.0])


    points = [ np.array([-10.0, 0.0, 0.0]), 
               np.array([-10.0, 1.0, 0.0]),
               np.array([-10.0, 0.0, 2.0]),
               np.array([10.0, -1.0, 0.0]),
               np.array([10.0, 0.0, -2.0]),
               np.array([10.0, 0.0, 0.0]),
               1.0 * blockDirector,
               2.0 * blockDirector,
               3.0 * blockDirector, 
               4.0 * blockDirector,
               blockRefPoint] 
    names = ['B', 'B', 'B', 'B', 'B', 'B', 'Ca', 'Ca', 'Ca', 'Ca', 'O']
    FIO.saveXYZList(points, names, 'initBlockPoints.xyz')
    
    labPoints = transformFromBlockFrameToLabFrame(labDirector, labRefPoint, labRotation, blockDirector, blockRefPoint, points)

    FIO.saveXYZList(labPoints, names, 'initlabPoints.xyz')

    newBlockPoints = transformFromLabFrameToBlockFrame(labDirector, labRefPoint, labRotation, blockDirector, blockRefPoint, labPoints) 
    
    FIO.saveXYZList(newBlockPoints, names, 'newBlockPoints.xyz')

    newLabPoints = transformFromBlockFrameToLabFrame(labDirector, labRefPoint, labRotation, blockDirector, blockRefPoint, newBlockPoints)

    FIO.saveXYZList(newLabPoints, names, 'newLabPoints.xyz')

    print("TNB Testing")
    # check TNB frame generation

    #v1 = np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1),  0])  # a random vector in the xyz plane
    #v2 = v1 + np.array([np.random.uniform(-10, -1), np.random.uniform(-10, -1),  0]) # add a random -ve nonzero displacements along the x-axis and y-axis to 
                                                            # yield a second vector that cannot be colinear with first but is still in xy plane

    #v3 = v2 + np.array([0, 0, np.random.uniform(1, 2) ]) # add a random z vector to v2 which has +ve z component 
    
    
    # position vectors for 3 particles
    v1 = np.array([0, 0, 0])
    v2 = np.array([1, 1, 0])
    v3 = 3 * v2 # construct colinear array
    # construct TNB frame for the end of this polymer of 3 particles
    TNB = constructTNBFrame(v1, v2, v3) # tangent should be llel with Z axis, B 
    print(TNB)
    
    # position vectors for 3 particles
    v1 = np.array([0, 0, 0])
    v2 = np.array([1, 1, 0])
    v3 = np.array([1, 1, 1])
    
    # construct TNB frame for the end of this polymer of 3 particles
    TNB = constructTNBFrame(v1, v2, v3) # tangent should be llel with Z axis, B 
    print(TNB)
    print("tHat.zHat = 1: ", np.dot(TNB[0], np.array([0,0,1]))) # should be a 1
    print( "tHat.tHat= 1: ", np.dot(TNB[0], TNB[0]) )
    print( "nHat.nHat= 1: ", np.dot(TNB[1], TNB[1]) )
    print( "bHat.bHat= 1: ", np.dot(TNB[2], TNB[2]) )
    print( "tHat.nHat= 0: ", np.dot(TNB[0], TNB[1]) )
    print( "nHat.bHat= 0: ", np.dot(TNB[1], TNB[2]) )
    print( "tHat.bHat= 0: ", np.dot(TNB[0], TNB[2]) )
    
    # test bond angle
    print("Bond Angle")
    angles = [-180, -135, -90, -45, 0, 45, 90, 135, 180, 225, 270, 315, 360]
    p1 = np.array([0,0,0])
    p2 = np.array([1,0,0])
    for angle in angles:
        p3 = np.array([1 - np.cos(angle * np.pi/180), np.sin(angle * np.pi/180), 0]) 
        print(angle,  bondAngle(p1, p2, p3) * 180.0/np.pi )    

    # test dihedral
    print("Dihedral Angle")    
    angles = [-180, -135, -90, -45, 0, 45, 90, 135, 180, 225, 270, 315, 360]
    p1 = np.array([0,0,0])
    p2 = np.array([1,0,0])
    p3 = np.array([1,0,1])
    
    for angle in angles:
        p4 = np.array([1 - np.cos(angle * np.pi/180), -np.sin(angle * np.pi/180), 1.0]) 
        print(angle,  Dihedral(p1, p2, p3, p4) * 180.0/np.pi )    
    

    # test rotatePointInChain
    points = [ np.array([-1,-1,0]), np.array([1, 1, 2]), np.array([2, 2, 0])]
    angle = 90 * np.pi /180
    newPoint = rotatePointInChain(points, angle)
    # preserve lengths and internal angle
    print(np.linalg.norm(points[1] - points[0]), np.linalg.norm(newPoint - points[0] ) )
    print(np.linalg.norm(points[2] - points[1]), np.linalg.norm(points[2] - newPoint ) )
    print("Original bond angle: " + str(bondAngle(points[0], points[1],  points[2]) * 180/np.pi ) )
    print("New bond angle: " + str( bondAngle(points[0], newPoint, points[2] )  * 180/np.pi ) )
    #export    
    points.append(newPoint)
    FIO.saveXYZ(points, 'C', "rotatePointsInChain.xyz")
    print("Done rotatePointsInChain")        
        
    # test rotate dihedral
    angle = 45 * np.pi/180
    changeDihedral = [ np.array([-2,-2,-2]),
                       np.array([2, -2, -2]),
                       np.array([2, -2, 2]), 
                       np.array([6, -2, 2])]
                                
    originalDihedral = Dihedral( changeDihedral[0],
                                 changeDihedral[1],
                                 changeDihedral[2],
                                 changeDihedral[3]) * 180/np.pi
    newPos = rotateDihedral(changeDihedral[1:4], angle)
    newDihedral = Dihedral( changeDihedral[0],
                            changeDihedral[1],
                            changeDihedral[2],
                            newPos) * 180/np.pi

    # dihedral should have changed by angle    
    print("original Dihedral: " + str(originalDihedral))
    print("new Dihedral: " + str(newDihedral))
    print("new Dihedral + angle: " + str(newDihedral + angle * 180/np.pi))
    
    print("BondLength: " + str(np.linalg.norm(changeDihedral[3] - changeDihedral[2])) +" " + str(np.linalg.norm(newPos - changeDihedral[2])))
    
    changeDihedral.append(newPos)
    FIO.saveXYZ(changeDihedral, 'C', "rotateDihedral.xyz")
    print("Done rotateDihedral")        
    
    
    print("Testing generateTNBVecXYZ")
    p1 = np.array([0,0,0])
    p2 = np.array([1,0,0])
    p3 = np.array([1,0,1])
    TNB_genTNBXYZ = constructTNBFrame(p1, p2, p3)
    
    # construct a 4th particle using the given alpha and beta which are the conventional dihedral and bond angle defining the 4th point.
    betas = [0, 45, 90, 135, 180] # beta is bond angle in range 0 to 180.
    alphas = [-180, -135, -90, -45, 0, 45, 90, 135, 180] # alpha is dihedral in range -180 to 180
    for alpha in alphas:
        for beta in betas:
            alphaP4 = alpha * np.pi/180 # dihedral
            betaP4 = beta * np.pi/180.0  # bond angle 
            p4 = p3 + generateTNBVecXYZ(TNB_genTNBXYZ, betaP4, alphaP4)
            newDihedral = Dihedral(p1, p2, p3, p4)
            if (newDihedral == None):
                newDihedral = "undefined"
            else:
                newDihedral = newDihedral  * 180 / np.pi
            
            newBondAngle = bondAngle(p2, p3, p4)
            if (newBondAngle==None):
                newBondAngle = "undefined"
            else:
                newBondAngle = newBondAngle * 180 / np.pi
            print("alpha: ", alpha, newDihedral, "beta: ", beta, newBondAngle, p4) 

    
    # testing Spherical Polar
    
    a = np.array( [2, 2, 2] )
    b = np.array( [2, 2, -2] )
    c = np.array( [2, -2, 2] )
    d = np.array( [2, -2, -2] )
    e = np.array([-2,2,2])
    f = np.array([-2,2,-2])
    g = np.array([-2,-2,2])
    h = np.array([-2,-2,-2])
    
    aSP = XYZ2SphericalPolar(a)
    bSP = XYZ2SphericalPolar(b)
    cSP = XYZ2SphericalPolar(c)
    dSP = XYZ2SphericalPolar(d)
    eSP = XYZ2SphericalPolar(e)
    fSP = XYZ2SphericalPolar(f)
    gSP = XYZ2SphericalPolar(g)
    hSP = XYZ2SphericalPolar(h)

    aXYZ = sphericalPolar2XYZ(aSP)
    bXYZ = sphericalPolar2XYZ(bSP)
    cXYZ = sphericalPolar2XYZ(cSP)
    dXYZ = sphericalPolar2XYZ(dSP)
    eXYZ = sphericalPolar2XYZ(eSP)
    fXYZ = sphericalPolar2XYZ(fSP)
    gXYZ = sphericalPolar2XYZ(gSP)
    hXYZ = sphericalPolar2XYZ(hSP)
    
    print("Spherical" )
    print(a, aSP, aXYZ )
    print( b, bSP, bXYZ )
    print( c, cSP, cXYZ )
    print( d, dSP, dXYZ )
    print( e, eSP, eXYZ )
    print( f, fSP, fXYZ )
    print( g, gSP, gXYZ )
    print( h, hSP, hXYZ )
    

    aCyl = XYZ2Cyl(a)
    bCyl = XYZ2Cyl(b)
    cCyl = XYZ2Cyl(c)
    dCyl = XYZ2Cyl(d)
    eCyl = XYZ2Cyl(e)
    fCyl = XYZ2Cyl(f)
    gCyl = XYZ2Cyl(g)
    hCyl = XYZ2Cyl(h)

    aXYZCyl = cyl2XYZ(aCyl)
    bXYZCyl = cyl2XYZ(bCyl)
    cXYZCyl = cyl2XYZ(cCyl)
    dXYZCyl = cyl2XYZ(dCyl)
    eXYZCyl = cyl2XYZ(eCyl)
    fXYZCyl = cyl2XYZ(fCyl)
    gXYZCyl = cyl2XYZ(gCyl)
    hXYZCyl = cyl2XYZ(hCyl)

    print( "Cylindrical" )

    print( a, aCyl, aXYZCyl )
    print( b, bCyl, bXYZCyl )
    print( c, cCyl, cXYZCyl )
    print( d, dCyl, dXYZCyl )
    print( e, eCyl, eXYZCyl )
    print( f, fCyl, fXYZCyl )
    print( g, gCyl, gXYZCyl )
    print( h, hCyl, hXYZCyl )
    
  
    
    alpha1 = -60 * np.pi/180.0
    alpha2 = 60 * np.pi/180.0
    
    beta1 = 100 * np.pi/180.0
    beta2 = 180 * np.pi/180.0
    
    tnbPoints = [ pickRandomTNBDirectionInAngRangeXYZ(TNB, beta1, beta2, alpha1, alpha2) for a in range(1000) ]
    tnbPointsXYZ = [ v3 + 10.0 * tnbPoint for tnbPoint in tnbPoints ] 
    tnbPointsXYZ.append(v1)
    tnbPointsXYZ.append(v2)
    tnbPointsXYZ.append(v3)
    FIO.saveXYZ(tnbPointsXYZ, 'C', "testTNBSelection.xyz")
    print("Done TNB Selection")
    
    # testing the random selection algorithms
    frustumPointsXYZ = [ pickRandomPointInFrustumXYZ(-10.0, 10.0, 5.0, 20.0) for a in range(5000) ]
    FIO.saveXYZ(frustumPointsXYZ, 'C', "testFrustumFull.xyz")
    print("Done Frustum Full")

    frustumPointsZHeightsXYZ = [ pickRandomPointInFrustumAtGivenZHeightXYZ(float(z), -10.0, 10.0, 5.0, 20.0) for z in range(-10, 12, 2) for b in range(100) ]
    FIO.saveXYZ(frustumPointsZHeightsXYZ, 'O', "testFrustumSelectedZ.xyz")
    print("Done Frustum Selected Z")


    ellipsoidPoints = [ pickRandomPointUniformlyOnEllipsoid(10.0, 15.0, 20.0) for a in range(1000) ]
    ellipsoidPointsXYZ = [ sphericalPolar2XYZ( point ) for point in ellipsoidPoints ]  
    FIO.saveXYZ(ellipsoidPointsXYZ, 'C', "testEllipsoidFull.xyz")
    print("Done Ellipsoids Full")
     
    ellipsoidPointsRange = [ pickRandomPointUniformlyOnEllipsoidInAngRange(10, 15, 20, -45 * np.pi/180, 45 * np.pi/180, -30 * np.pi/180, 30 * np.pi/180) for a in range(1000) ]
    ellipsoidPointsRangeXYZ = [ sphericalPolar2XYZ( point ) for point in ellipsoidPointsRange ]  
    FIO.saveXYZ(ellipsoidPointsRangeXYZ, 'C', "testEllipsoidRange.xyz")
    print("Done Ellipsoids range")
    
    spherePoints = [ pickRandomPointOnUnitSphere() for a in range(1000) ]
    spherePointsXYZ = [ sphericalPolar2XYZ( np.array( [10, point[0], point[1]] ) ) for point in spherePoints ]  
    FIO.saveXYZ(spherePointsXYZ, 'C', "testSphereFull.xyz")
    print("Done Sphere Full")
        
    spherePointsRange = [ pickRandomPointOnUnitSphereInAngRange(-45 * np.pi/180, 45 * np.pi/180, -30 * np.pi/180, 30 * np.pi/180) for a in range(1000) ]
    spherePointsRangeXYZ = [ sphericalPolar2XYZ( np.array( [10, point[0], point[1]] ) ) for point in spherePointsRange ]  
    FIO.saveXYZ(spherePointsRangeXYZ, 'C', "testSphereRange.xyz")
    print("Done Sphere Range")

    circlePoints = [ pickRandomPointOnUnitCircle() for a in range(1000) ]
    circlePointsXYZ  = [ np.array([10 * np.cos(point), 10*np.sin(point), 0.0]) for point in circlePoints ] 
    FIO.saveXYZ(circlePointsXYZ, 'C', "testCircleFull.xyz")
    print("Done Circle Full")
    
    circlePointsRange = [ pickRandomPointOnUnitCircleInAngRange(0 * np.pi/180, 90 * np.pi/180) for a in range(1000) ]
    circlePointsRangeXYZ  = [ np.array([10 * np.cos(point), 10*np.sin(point), 0.0]) for point in circlePointsRange ] 
    FIO.saveXYZ(circlePointsRangeXYZ, 'C', "testCircleRange.xyz")
    print("Done Circle Range")
    
    print("done")
    
    