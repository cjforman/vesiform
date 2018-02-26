'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB
from Library.peptideBackbone import peptideBackboneGenerator as PBG
import sys
import numpy as np
import random as rnd
from Utilities import coordSystems as coords
from Utilities import fileIO as fIO

class ConstrainedPolymerPackNBB(NBB):
    '''
    Overloads the NPackBB object to give a building block with a packing which results in a polymeric chain of N monomers
    with random dihedrals and bond angles, whose ends are constrained by a collection of six points.
    This list of six point defines dihedrals and bond angles for the end points, and allows a 
    random coil to be inserted easily into an exterior structure.  The two end points of the returned
    chain are the positions of two of these fixed points. Thus the algorithm must define N-2 new points
    consistent with the constraints.
    
    Each new monomer is added by the NPack polymer algorithm which only chooses points yielding dihedrals and bonds
    within the given range of dihedral and bond angles already. Also ensures that the polymer does not self intersect.
    
    To accelerate the method a region of space is defined outwith the new point cannot go.
    This is a sphere centered on the fixed end point, if the new tip of the polymer strays outside this region,
    then there is no possible way the remaining links could reach the target point so the point is rejected. 
    In practice a compression value is applied to the sphere radius to account for the fact that the end 
    to end distance is a fraction of the contour length.  
    
    There is also an inner and an outer restraining sphere which allows a definition of parts of space where the
    chain cannot go. This allows emulation of steric hindrances allowing the random coil to fit in around other
    parts of a structure.  
    
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)

    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        self.compressionScaleFactor = self.getParam('compressionScaleFactor')
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for constrained polymer object"
            sys.exit()

    def generateBuildingBlock( self, 
                               numPoints, 
                               pointA,
                               pointB, 
                               alpha1, 
                               alpha2,
                               beta1,
                               beta2, 
                               minDist,
                               bondLength,
                               innerSphereR,
                               outerSphereR,
                               spherePos):

        self.originalNumPoints = numPoints
        self.pointA = pointA
        self.pointB = pointB
        self.bondLength = bondLength
        
        # use these ranges for the dihedral and bond angles to find an intial chain between
        # the end points. We will apply the bonds and dihedrals constraints once we have a starting chain.
        self.alpha1 = alpha1 * np.pi/180
        self.alpha2 = alpha2 * np.pi/180
        self.beta1 = beta1 * np.pi/180
        self.beta2 = beta2 * np.pi/180
        if self.alpha2 < self.alpha1:
            self.alpha2 = self.alpha1
            self.alpha1 =  alpha2 * np.pi/180.0
        if self.beta2 < self.beta1:
            self.beta2 = self.beta1
            self.beta1 = beta2 * np.pi/180.0

        self.bannedSphereCentre = spherePos
        self.bannedSphereRadiusInner = innerSphereR 
        self.bannedSphereRadiusOuter = outerSphereR

        return NBB.generateBuildingBlock(self, numPoints, minDist)
    
    def generateBuildingBlockDirector(self):
        vec = self.pointB - self.pointA
        return vec/np.linalg.norm(vec)
    
    def generateBuildingBlockRefPoint(self):
        return self.buildingBlockXYZ[0]

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.originalNumPoints

    def generateNConstrainedRandomPoints(self):
        # compute some derived quantities
        self.maxLength = (self.numPoints-1) * self.bondLength
        self.distAToB= np.linalg.norm(self.pointB - self.pointA)
        self.polymerCompression= self.compressionScaleFactor * self.distAToB/self.maxLength

        # perform an initial sanity check on the given starting points
        if (not self.checkPointInSphereRange(self.pointA)):
            print ("Invalid starting point: pointA is not inside spherical constraints.")
            sys.exit()
        
        if (not self.checkPointInSphereRange(self.pointB)):
            print ("Invalid starting point: pointB is not inside spherical constraints.")
            sys.exit()        
        
        # calculate the maximum allowed distance that the two end points can be apart
        maxDistance = self.polymerCompression * (self.numPoints - 1) * self.bondLength
        if np.linalg.norm(self.pointA - self.pointB) > maxDistance:
            print ("Invalid starting point: Point A and B are too far apart for given number of points, and bondLength.")
            sys.exit()        

        return NBB.generateNConstrainedRandomPoints(self)
      
    def checkBondLengths(self):
        print [ np.linalg.norm(self.buildingBlockXYZ[i] - self.buildingBlockXYZ[i-1]) for i in range(1, len(self.buildingBlockXYZ)) ]       
        
    def pickFirstPoints(self):
        
        # The first point in the chain will be pointA.
        # must calculate two additional points near point A 
        # which will form the start of the chain. These must satisfy the bond angle
        # requirements to prevent the algorithm from stalling.
        
        # reset the numPoints to original number - this function increments it so it is 
        # responsible for starting with the right value
        self.numpoints = self.originalNumPoints 
        
        # a point bondLength below A on the Z axis.
        A2 = self.pointA - np.array([0, 0, self.bondLength])
        
        # pick a random bond angle within range
        bondAngle = rnd.uniform(self.beta1, self.beta2)
        
        # A1 is a point which is in the XZ plane subtending an angle of bondAngle
        # with point A2 and point A 
        A1 = A2 - np.array([self.bondLength * np.sin(np.pi- bondAngle),
                            0.0,
                            self.bondLength * np.cos(np.pi- bondAngle)])
        
        # begin the list
        self.nList = [ A1, A2, self.pointA ]
        
        # add two to the target length of the chain to account for the fact that
        # we have manually added two points to the list.
        self.numPoints = self.numPoints + 2
        
        # subtract one from the target length of the chain to account for the fact that
        # we will manually add the last point at the end
        self.numPoints = self.numPoints - 1  
        
        return True
        
    def getNList(self):
        # Return the nlist without the first two points and adding the last point manually.
        # The first two points in the NList are dummy points to 
        # get the dihedrals/bond angles right of the first couple of points in the output
        # this should return a list of self.numPoints points which was defined by the user in the input file. 
        self.nList.append(self.pointB)
        if len(self.nList[2:]) != self.originalNumPoints:
            print("Constrained Polymer Warning: NList is the wrong length in getNList")
        
        return self.nList[2:]  
    
    def checkPointInBounds(self, pos):
        # assume point is inBounds and test for outBoundness
        outBounds = False
        
        # calculate how many links are left to compute 
        numLinksLeft =  self.originalNumPoints - (len(self.nList) - 2 )  # subtract the first two points which don't count.
        #print "numLinksLeft: ", numLinksLeft, "len nList: ", len(self.nList)
        # calculate the maximum allowed distance that the current point can be from the end point
        maxDistance = self.polymerCompression * numLinksLeft * self.bondLength
    
        # check to see if we are inserting the last point. Must do some additional checks.
        if (numLinksLeft == 2):
            # check the bond angle of the last three elements in list will be within range
            #inBounds=self.checkBondInBounds([self.nList[-1], pos, self.pointB])
            #if not inBounds:
            #    outBounds = True 

            # check that the dihedral of the last four atoms is in range            
            if (outBounds==False):
                inBounds=self.checkDihedralInBounds([self.nList[-2], self.nList[-1], pos, self.pointB])
                if not inBounds:
                    outBounds = True
            
            # check that the penUltimate ball is precisely bondLength from the last point 
            if (outBounds==False):
                if ( np.abs(np.linalg.norm(pos - self.pointB) - self.bondLength) > 1e-2): 
                    outBounds = True


        if (outBounds==False) and numLinksLeft>2:                
            # compute distance of last point on chain from end point and compare to the maximum distance
            if np.abs(np.linalg.norm(pos - self.pointB)) > maxDistance:
                outBounds = True 
                if self.verbose==1: 
                    print("dynamic radius violation")             

        # Only do this check if none of the previous checks failed      
        if (outBounds == False):
            # check for a violation of the bulk exclusion space 
            if np.linalg.norm(pos - self.bannedSphereCentre) < self.bannedSphereRadiusInner:
                outBounds = True 
                if self.verbose==1: 
                    print("steric sphere inner violation")             

        # Only do this check if none of the previous checks failed      
        if (outBounds == False):
            if np.linalg.norm(pos - self.bannedSphereCentre) > self.bannedSphereRadiusOuter:
                outBounds = True 
                if self.verbose==1: 
                    print("steric sphere outer violation")             

        # Only do this check if none of the previous checks failed      
        if (outBounds == False) and numLinksLeft>2:
            inBounds=self.checkBondInBounds([self.nList[-2], self.nList[-1], pos])
            if not inBounds:
                outBounds = True 
        
        # Only do this check if none of the previous checks failed      
        if (outBounds == False) and numLinksLeft>2:
            inBounds=self.checkDihedralInBounds([self.nList[-3], self.nList[-2], self.nList[-1], pos])
            if not inBounds:
                outBounds = True 
            
        if outBounds == True:
            inBounds = False
        else:
            inBounds = True
        return inBounds

    def checkPreformedChainParticleInBounds(self, pointList):
        # Returns true:
        #     if particle 3 in pointList yields bond angles and dihedrals in the appropriate ranges,
        #     if the particle is not inside or outside the forbidden spheres, 
        #     if the particle is not within minDist of all the other points. 
        #     
        #     only does a check if it needs to.
        inBounds = self.checkBondInBounds(pointList[1:4])
        if inBounds:
            inBounds = self.checkBondInBounds(pointList[3:6])
        if inBounds:
            inBounds = self.checkDihedralInBounds(pointList[0:4])
        if inBounds:
            inBounds = self.checkDihedralInBounds(pointList[3:7])
        if inBounds:
            inBounds = self.checkPointInSphereRange(pointList[3])
        if inBounds:
            inBounds = self.checkPointAgainstList(pointList[3])

        return inBounds
            
    def checkPointInSphereRange(self, pos):
        # assume we're not out of bounds
        outBounds = False
        
        # compute dist between pos and sphere centre.
        dist = np.linalg.norm(pos - self.bannedSphereCentre)

        # check dist against radii of two spheres
        if ( dist < self.bannedSphereRadiusInner):
            outBounds = True 
            if self.verbose==0: 
                print("Steric Sphere inner violation")             
    
        # check dist against radii of two spheres
        if ( dist > self.bannedSphereRadiusOuter):
            outBounds = True 
            if self.verbose==0: 
                print("Steric Sphere outer violation")
        
        return not outBounds        
    
    def checkDihedralInBounds(self, pointList):
        outBounds = False

        if len(pointList)==4:
                    
            dihedral = coords.Dihedral( pointList[ 0 ],
                                        pointList[ 1 ], 
                                        pointList[ 2 ], 
                                        pointList[ 3 ] )
    
            if (dihedral < self.alpha1) or (dihedral > self.alpha2):
                outBounds = True 
                if self.verbose==1: 
                    print("Dihedral violation")             
        else:        
            print ("Need list of four points to perform dihedral check.")
            outBounds = None

        return not outBounds

    def checkBondInBounds(self, pointList):
        # returns true if the bond angle of 3 points is within the given ranges
        outBounds = False

        if len(pointList)==3:
                    
            bondAngle = coords.bondAngle( pointList[ 0 ],
                                          pointList[ 1 ], 
                                          pointList[ 2 ])
    
            if (bondAngle < self.beta1) or (bondAngle > self.beta2):
                outBounds = True 
                if self.verbose==1: 
                    print("Bond Angle violation")
        else:        
            print ("Need list of three points to perform bond angle check.")
            outBounds = None

        return not outBounds

        
    # checks that the dihedrals involving a point in the centre of a list of seven points are all within the specified bounds
    def checkDihedrals(self, pointList):

        outBounds = False

        if len(pointList)==7:
                    
            dihedrals = [coords.Dihedral( pointList[ n ],
                                          pointList[ n + 1 ], 
                                          pointList[ n + 2 ], 
                                          pointList[ n + 3 ] ) for n in range(0,3)]
    
            if True in [ (dihedral < self.alpha1) or (dihedral > self.alpha2) for dihedral in dihedrals ]:
                outBounds = True 
                if self.verbose==0: 
                    print("Dihedral violation")             

        else:        
            print ("Need list of seven points to perform dihedral checks.")
            outBounds = None

        return outBounds

    # checks that the 3 bond angles involving a point in the centre of a list of five points are all within the bounds
    def checkBonds(self, pointList):

        outBounds = False
        
        if len(pointList) == 5:
                
            bondAngles = [ coords.bondAngle( pointList[ n ], 
                                             pointList[ n + 1 ], 
                                             pointList[ n + 2 ] ) for n in range(0,2) ]
            
            if True in [ (bondAngle < self.beta1 ) or (bondAngle > self.beta2) for bondAngle in bondAngles]:
                outBounds = True 
                if self.verbose==0: 
                    print("Bond angle violation")
        else:        
            print ("Need list of five points to perform bond checks.")
            outBounds = None
            
        return outBounds
             
             
    def pickRandomPointInDefinedSpace(self):
        # Takes the last three points of the nList and picks a new point in XYZ which 
        # is inside the given angular ranges. No need to check the ranges, because the points are only selected inside those zones.
        
        # construct the TNB vectors from the last three points of the nList.
        TNB = coords.constructTNBFrame(self.nList[-3], self.nList[-2], self.nList[-1])        
        
        # compute a directional vector in XYZ coords from the last point on the list to the new point, within the given angular range constraints. 
        newCoordXYZRel = coords.pickRandomTNBDirectionInAngRangeXYZ(TNB, 
                                                                    self.beta1, 
                                                                    self.beta2,
                                                                    self.alpha1,
                                                                    self.alpha2)
        
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

    # create a beta backbone generator object using static file parameters
    backboneObject = PBG(filename)

    # generate backbone realtime parameters
    numPos = 3* 8
    startPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    offset =  np.array([4.8, 0.0, 0.0])
    polarity = 'CN'
    
    backBoneBuildingBlock1 = backboneObject.generateBuildingBlock(numPos, polarity) 
    backBoneBuildingBlock1.transformFromBlockFrameToLabFrame(director, startPos, rotation)
    backBoneBuildingBlock1.exportBBK("backbone1")
    backBoneBuildingBlock2 = backboneObject.generateBuildingBlock(numPos, polarity) 
    backBoneBuildingBlock2.transformFromBlockFrameToLabFrame(director, startPos + offset, rotation)
    backBoneBuildingBlock2.exportBBK("backbone2")

    startPoints = backBoneBuildingBlock1.xyzVals[3:]
    endPoints = backBoneBuildingBlock2.xyzVals[3:]
    
    # create the NPack object.
    ConstrainedPolymerPackGBB = ConstrainedPolymerPackNBB(filename)
    
    numPoints = 10 
    pointA = startPoints[-1] 
    pointB = endPoints[-1]
    alpha1 = -180
    alpha2 = 180
    beta1 = 0
    beta2 = 180
    minDist = 1.0
    bondLength = 1.5
    innerSphereR = 0.99 * np.linalg.norm(pointA-pointB)/2
    outerSphereR = 100
    spherePos = (pointA + pointB)/2 
       
    # generate the XYZVals in the packed spaced - 
    # this generates a self-avoiding chain between two points that is contained inside a spherical shell. 
    # The bond angles and dihedrals are all random from the full range of angles and dihedrals.
    ConstrainedPolymerPackBB = ConstrainedPolymerPackGBB.generateBuildingBlock(numPoints, pointA, pointB, alpha1, alpha2, beta1, beta2, minDist, bondLength, innerSphereR, outerSphereR, spherePos)
    ConstrainedPolymerPackGBB.checkBondLengths()
    ConstrainedPolymerPackBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))