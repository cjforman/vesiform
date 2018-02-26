'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB
import sys
import numpy as np
from Utilities import coordSystems as coords
from Utilities import cartesian as cart
from Utilities import fileIO as fIO
import random as rnd


class RandomPolymerPackNBB(NBB):
    '''
    Over loads the NPack object to give a packing which results in a random polymeric chain. 
    # Uses the NPack algorithm to randomly choose a dihedral angle for a new monomer on the end of a chain
    # within certain angular constraints.
    # also ensures that the polymer does not self intersect, and does not exceed an outer envelope
    # which is a cylindrical frustum. 
    # the frustum is parallel with the z-axis has an upper and low radii and a defined z length. 
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
   
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Random Polymer Pack"
            sys.exit()                    

    def generateBuildingBlock( self, 
                               numPoints, 
                               rInner, 
                               rOuter, 
                               z1,
                               z2, 
                               alpha1, 
                               alpha2,
                               beta1,
                               beta2, 
                               minDist,
                               bondLength,
                               seedZ='lower',
                               envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):

        self.radiusInner = rInner
        self.radiusOuter = rOuter

        # Set length of outer frustum
        self.z1 = z1
        self.z2 = z2
        if (self.z2 < self.z1):
            self.z2 = self.z1
            self.z1 = z2

        # compute the length of the frustum
        self.zLength = abs(self.z2 - self.z1) 

        # define which plane the start point is in: lower, center, upper, random
        self.seedZ = seedZ.lower()

        # define the bond length
        self.bondLength = bondLength

        # Process the angle ranges that defining angles between all the bonds. Angles in degrees in the file, but always radians in the code.
        # Alpha is the allowed range of dihedral angles with respect to the last three points on the polymer when we add a new point
        # beta is the allowed range of bond angles with respect to the last two point on the polymer chain when we add a new point  
        # Alpha can range from -pi to pi.  The user specifies beta in the range  0 to pi
        # but we subtract pi/2 so we can use the same spherical polar functions used throughout the project. 
        # Thus, within the code beta = 0 gives an internal bond angle of ninety. 
        self.alpha1 = alpha1 * np.pi/180.0
        self.alpha2 = alpha2 * np.pi/180.0
        if self.alpha2 < self.alpha1:
            self.alpha2 = self.alpha1
            self.alpha1 = alpha2 * np.pi/180.0

        self.beta1 = beta1 * np.pi/180.0
        self.beta2 = beta2 * np.pi/180.0
        if self.beta2 < self.beta1:
            self.beta2 = self.beta1
            self.beta1 = beta2 * np.pi/180.0
        
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockDirector(self):
        return cart.getPrincipalAxis(self.buildingBlockXYZ)

    def generateBuildingBlockRefPoint(self):
        return self.buildingBlockXYZ[0]

    def pickFirstPoints(self):
        # pick three starting points 
        
        # choose positions that are a little bit away from the frustum boundary
        
        # choose a starting z position in bands (stay away from frustum boundary)
        if self.seedZ=='lower':
            zHeight = self.z1 + rnd.uniform( .05 * self.zLength, 0.15 *self.zLength) 
        elif self.seedZ=='centre':
            zHeight = self.z1 + rnd.uniform( .45 * self.zLength, 0.55 *self.zLength)
        elif self.seedZ=='upper':
            zHeight = self.z1 + rnd.uniform( .85 * self.zLength, 0.95 *self.zLength)
        else:
            zHeight = rnd.uniform( .05 * self.zLength, 0.95 *self.zLength)

        # choose a random point on the selected z plane
        self.nList = [ coords.pickRandomPointInFrustumAtGivenZHeightXYZ(zHeight, self.z1, self.z2, 0.9 * self.radiusInner, 0.9 * self.radiusOuter) ]
        
        # choose a second random point within the frustum but at the same plane as the first 
        self.nList.append(coords.pickRandomPointInFrustumAtGivenZHeightXYZ(self.nList[0][2], self.z1, self.z2, 0.9 * self.radiusInner, 0.9 * self.radiusOuter))
    
        # set up the variables for third point
        thirdPointInBounds = False
        numLoops = 0 
        
        # loop until we find a good third point
        while thirdPointInBounds==False and numLoops < self.maxLivesPerNewNode:

                theta, phi = coords.pickRandomPointOnUnitSphere()
                if (self.seedZ== 'lower'):
                    # get a new vector on the unit sphere
                    theta, phi = coords.pickRandomPointOnUnitSphereInAngRange(0, 90*np.pi/180, -np.pi, np.pi)

                if (self.seedZ== 'upper'):
                    # get a new vector on the unit sphere
                    theta, phi = coords.pickRandomPointOnUnitSphereInAngRange(-90*np.pi/180, 0, -np.pi, np.pi)

                # compute xyz position of new vector
                thirdPosXYZ = self.nList[1] + coords.sphericalPolar2XYZ(np.array([self.bondLength, theta, phi]))
            
                # check the third pos is inside the frustum.
                thirdPointInBounds = self.checkPointInBounds(thirdPosXYZ)

                numLoops += 1

        if thirdPointInBounds:
            # append the third point to the nList
            self.nList.append(thirdPosXYZ)
            
            # now we have three points we can use the main algorithm to find the fourth.

            # assume that the new point is going to be out of bounds
            fourthPointInBounds = False
                
            # count number of times we tried to find a good point 
            numLoops = 0
                
            # trying to find a good new point up to maxLivesPerNewNode times
            while ((fourthPointInBounds==False) and (numLoops < self.maxLivesPerNewNode)): 
                    
                # pick a new point in the defined sub space.
                newPointInSpace = self.pickRandomPointInDefinedSpace()

                # check the new Point is within the acceptable part of the defined space. 
                # (in this case that point is inside the frustum)
                fourthPointInBounds = self.checkPointInBounds(newPointInSpace)
                
                # increment counter
                numLoops += 1
                    
                # if we found a good new point then perform a global check.               
                if fourthPointInBounds==True:
                      
                    # convert from point in space to XYZ using user defined function.
                    # useful for generating chains
                    # for polymers: In fact we no longer use this check to do anything since the point is already in xyz space
                    # thanks to new forumulation
                    newPointInXYZ = self.convertPointToXYZ(newPointInSpace)
    
                    # check to see if the new point satisfies the global constraints
                    fourthPointInBounds = self.checkPointAgainstList(newPointInXYZ)
    
        if fourthPointInBounds==True:
            self.nList.append(newPointInXYZ) # add the last xyz Value
            self.nList.pop(0) # eliminate the first point, which was just a construction vector to get us going and not to be included in the polymer
            self.nAttempts = [0, 0, 0]
            return True
        else:
            print "Unable to Initialise List. Try increasing maxLivesPerNewNode."
            return False
    
    def checkPointInBounds(self, pos):
        # pos in XYZ coords - must check it is inside frustum defined in parameter file

        # assume it is inBounds
        inBounds = True

        # convert pos to cylindrical coords
        posCyl = coords.XYZ2Cyl(pos)
        
        # check the zValue is in range.
        if (posCyl[2] < self.z1) or (posCyl[2] > self.z2):
            inBounds=False
            if self.verbose==1:
                print("Frustum Z violation")             

        # compute the radius of the frustum at the height of the Z-coord.     
        radAllowed = coords.FrustumRadiusAtGivenZ(posCyl[2], self.z1, self.z2, self.radiusInner, self.radiusOuter)
        
        # if the radius of the test point exceeds the allowed radius then reject it.
        if (posCyl[0] > radAllowed):
            inBounds = False
            if self.verbose==1: 
                print("Frustum radial violation")
        
        
        if inBounds:
            inBounds = NBB.checkPointInBounds(self, pos)
        
        return inBounds
             
    def pickRandomPointInDefinedSpace(self):
        # Takes the last three points of the nList and picks a new point in XYZ which 
        # is inside the given angular ranges. No need to check the ranges, because the points are only selected inside those zones.
        
        # construct the TNB vecotrs from the last three points of the nList.
        TNB = coords.constructTNBFrame(self.nList[-3], self.nList[-2], self.nList[-1])        
        
        # compute a directional vector in XYZ coords from the last point on the list to the new point, within the given angular range constraints. 
        newCoordXYZRel = coords.pickRandomTNBDirectionInAngRangeXYZ(TNB, self.beta1, self.beta2, self.alpha1, self.alpha2)
        
        # should be a unit vector already but just in case
        newCoordXYZRelHat = newCoordXYZRel/np.linalg.norm(newCoordXYZRel) 
        
        # scale the vector by the bondlength and add it to the last point in the list to give new position vector
        return self.nList[-1] + self.bondLength * newCoordXYZRelHat
        
    def getParams(self):
        return self.params 


    def checkPointAgainstList(self, pos):
        ''' Function tests to see if new position is within minDist distance of all the other points.
            Assumes that pos is not in the list already. '''
        inBounds = True
    
        # if the point is inside the frustum then check to see if it is within minDist of any other 
        # particles.
        if inBounds:
            for zPos in self.nList:
                # if new position occurs within minDist of any of the existing positions in the list 
                # then throw a wobbly.
                if np.linalg.norm((zPos - pos)) <= self.minDist:
                    inBounds= False
                    if self.verbose==1: 
                        print("Packing violation")

        return inBounds

        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    RandomPolymerBBG = RandomPolymerPackNBB(filename)

    numPoints = 100
    centerPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0.0
    rInner = 2
    rOuter = 10
    z1 = 0
    z2 = 150    
    phi1 = -90
    phi2 = 90
    minDist = 1.5
    alpha1 = -90
    alpha2 = 90
    beta1 = 100
    beta2 = 180
    bondLength = 2
    seedZ = 'lower'
       
    # generate the XYZVals in the packed spaced
    RandomPolymerBB = RandomPolymerBBG.generateBuildingBlock(numPoints, rInner, rOuter, z1, z2, alpha1, alpha2, beta1, beta2, minDist, bondLength, seedZ)
    RandomPolymerBB.transformBBToLabFrame(director, centerPos, rotation)
    # export
    RandomPolymerBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))    