'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
import sys
import numpy as np
from scipy.optimize import minimize
import random as rnd
from Utilities import coordSystems as coords
from Utilities import fileIO as fIO
from Library.peptideBackbone import peptideBackboneGenerator as PBG

class ConstrainedPolymerPackBBG(BBG):
    '''
    Overloads the BBG object to give a building block with a packing which results in a polymeric chain of N monomers
    with random dihedrals and bond angles, whose ends are constrained by a collection of six points.
    This list of six point defines dihedrals and bond angles for the two end points, and allows a 
    random coil to be inserted easily into an larger structure.  The algorithm defines N points that 
    form a chain with a given bondlength. 
    
    
    '''
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        BBG.initialiseParameters(self) 
        
        self.compressionScaleFactor = self.getParam('compressionScaleFactor')
        self.maxNumAttempts = self.getParam('maxNumAttempts')
        self.maxLivesPerNewNode = self.getParam('maxLivesPerNewNode')
        self.maxLivesPerTwoNodes = self.getParam('maxLivesPerTwoNodes')
        self.particleName = self.getParam('particleName')
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for constrained polymer object"
            sys.exit()

    def generateBuildingBlock( self, 
                               numPoints,
                               startPoints,
                               endPoints,
                               alpha1, 
                               alpha2,
                               beta1,
                               beta2, 
                               minDist,
                               bondLength,
                               pointsToAvoid=[]):

        self.numPoints = numPoints
        self.startPoints = startPoints
        self.endPoints = endPoints
        self.bondLength = bondLength
        self.minDist = minDist
        self.pointsToAvoid = pointsToAvoid
        
        # use these ranges for the allowed dihedral and bond angles
        # when it comes to the last point or couple of points at add, 
        # we can relax these contraints for one point, providing the bond length is ok. 
        self.alphaMin = min([alpha1 * np.pi/180, alpha2 * np.pi/180])
        self.alphaMax = min([alpha1 * np.pi/180, alpha2 * np.pi/180])
        self.betaMin = min([beta1 * np.pi/180, beta2 * np.pi/180])
        self.betaMax = max([beta1 * np.pi/180, beta2 * np.pi/180])
        
        # compute the overall polymerCompression that will result
        self.MaxLength = (self.numPoints + 1) * bondLength
        self.distAToB = np.linalg.norm(startPoints[-1] - endPoints[-1])
        self.polymerCompression= self.compressionScaleFactor * self.distAToB/self.MaxLength
        
        return BBG.generateBuildingBlock(self, numPoints)

    def generateBuildingBlockDirector(self):
        vec = self.buildingBlockXYZ[-1] - self.buildingBlockXYZ[0]
        return vec/np.linalg.norm(vec) 
    
    def generateBuildingBlockRefPoint(self):
            # set the reference point to be the point at the end of the startPoints array
            # (this is the same as the start point passed to the generateBuildingBlock function).
            # this way the translation at the end of the buildingBlock function won't have any affect
            # on the position of the selected points. We do not want to rotate or move the coild
            # as it is carefully positioned "just so".
        return self.startPoints[-1]

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints

    def checkBondLengths(self):
        print [ np.linalg.norm(self.buildingBlockXYZ[i] - self.buildingBlockXYZ[i-1]) for i in range(1, len(self.buildingBlockXYZ)) ]       
        
    def generateBuildingBlockXYZ(self):

        nAttempts = 0
        listsAreGood = False # assume failure
        initialStartLength = len(self.startPoints)
        initialEndLength = len(self.endPoints)
        while not listsAreGood and nAttempts < self.maxNumAttempts:
            print "nAttempts: ", nAttempts, "of ", self.maxNumAttempts
            # call the function to recursively build up the polymer from both end points.
            startPointsFinal, endPointsFinal,  listsAreGood = self.pickTwoPoints(self.startPoints[:], self.endPoints[:], self.pointsToAvoid, self.numPoints)
            nAttempts += 1

        if not listsAreGood:
            print "Unable to compute constrained polymer"
            sys.exit()

        # if we survived then chop off the initial input points from both lists.
        # (another way to add pointsToAvoid is by sneaking them in in the startPoint and endPoints lists 
        # Reverse the endPoints list.
        # Add them together to make a single polymer
        out1 = startPointsFinal[initialStartLength:]
        out2 = endPointsFinal[initialEndLength:]
        
        return out1 + list(reversed(out2)) 

    def checkPointAgainstList(self, pos, testList):
        ''' Function returns true if the test position is greater than minDist distance from all the other test points.
            Assumes that pos is not in the list already.'''
        # assume position is good
        goodPos = True

        lenList = len(testList)
        curPoint = 0
        while goodPos and curPoint<lenList:
            # if new position occurs within minDist of any of the existing positions in the list 
            # then throw a wobbly.
            if np.linalg.norm((testList[curPoint] - pos)) <= self.minDist:
                goodPos = False
            curPoint += 1
        return goodPos

    def findOnePoint(self, searchPoints, auxPoints):
        # find a point inside the angle and dihedral bounds wrt the end of the searchPoints list.
        # Also check to ensure that the new point does not violate any of the points in searchPoints or auxPoints 
        inBounds = False
        numLoops = 0
        while not inBounds and numLoops<self.maxLivesPerNewNode:
            newPoint = self.pickRandomPointInDefinedSpace(searchPoints[-3], 
                                                          searchPoints[-2],
                                                          searchPoints[-1]) 

            # return false if there's an issue
            inBounds = self.checkPointAgainstList(newPoint, searchPoints + auxPoints)
            numLoops +=1
        
        # if we hit maxLivesPerNewNode loops without finding a suitable point then inBound is returned as False
        return newPoint, inBounds
        
    def pickTwoPoints(self, startPoints, endPoints, pointsToAvoid, numPointsLeft):
        # The core recursive function adds a new valid point to the two input lists and then calls itself 
        # with the numPoints left incremented.  The function returns a true or false flag to indicate success or failure.
        # a special function is used when numPoints left is one or two to terminate the chain correctly.  
        # 
        # To be valid the two points must have the correct bond angle, dihedral and bondLength wrt to 
        # the last three points in the given start and end points list.
        # Both points must be within certain distance of each other based on numPointsLeft.
        # Both points must not be within mindist of any other points on either list,
        # including each other or the self.pointsToAvoid member variable.

        print "Num Points Left: ", numPointsLeft

        listsAreGood = False # assume failure

        # compute max length of a polymer stretched out with the current number of points left. If
        # the distance between the two points exceeds this then the chain can never connect.
        # However the bond angles must also be taken in to account, so the real value to not
        # exceed is smaller than this.
        maxLength = (numPointsLeft - 1 ) * self.bondLength

        if numPointsLeft>2:
            pointAInBounds = False
            pointBInBounds = False
            numOuterLoops = 0
            while (not pointAInBounds or not pointBInBounds) and numOuterLoops<self.maxLivesPerTwoNodes:
                # print numOuterLoops
                # Get a candidate for point A that is within bounds, if we do not have one in hand already.
                if not pointAInBounds:
                    pointA, pointAInBounds = self.findOnePoint(startPoints, endPoints + pointsToAvoid)
        
                # if we found a pointA, and need a pointB then look for a pointB in same way.
                if pointAInBounds and not pointBInBounds:
                    pointB, pointBInBounds = self.findOnePoint(endPoints, startPoints + pointsToAvoid)
                    
                # now compare conditions dependent on both pointA and pointB (distance)
                if pointAInBounds and pointBInBounds:
                    distAToB = np.linalg.norm(pointB - pointA)
                    # check to see if new points are outside convergence window
                    if (distAToB > self.polymerCompression * maxLength) or distAToB<self.minDist:
                        #print "Distance Violation", distAToB, self.polymerCompression * maxLength
                        # too far apart or too close together.
                        # no point in looking for both a new pointA and pointB
                        # only need to look for a new value for one of them 
                        if numOuterLoops % 2 == 0:
                            pointBInBounds=False # on even loops look for a new point B 
                        if numOuterLoops % 2 == 1:
                            pointAInBounds=False # on odd loops look for a new point A 
                numOuterLoops +=1
            
            if pointAInBounds and pointBInBounds:
                # we successfully found two points that satisfy the contraints
                # call the function again to compute the next two points. 
                # after having added the two new points to the start and end point arrays
                startPoints.append(pointA)
                endPoints.append(pointB)
                newNumPointsLeft = numPointsLeft - 2
                
                # function returns lists with the new end points added on and a flag to indicate success or failure
                startPoints, endPoints, listsAreGood = self.pickTwoPoints(startPoints, endPoints, pointsToAvoid, newNumPointsLeft) 
            
        # if num points left is 2 then we have a special terminating case
        if numPointsLeft==2:
            pointA, pointB, listsAreGood = self.findLastTwoPoints(startPoints, endPoints, pointsToAvoid)
            if listsAreGood:
                startPoints.append(pointA)
                endPoints.append(pointB)
            
        # if num points let is 1 then we have a special case
        if numPointsLeft==1:
            pointA, listsAreGood = self.findLastPoint(startPoints, endPoints, pointsToAvoid)
            if listsAreGood:
                startPoints.append(pointA)
            
        return startPoints, endPoints, listsAreGood

    def whizzCheck(self, startPoints, endPoints, pointsToAvoid, newPoint):
        # check newPoint is not within minDist of any of the testPoints 
        # if it is then *whizz it about and hopefully find something that satisfies all the constraints, 
        # even though there's not many options available.
        # A whizz is a rotation of the test point about the axis between start and end point.
        # Such a rotation preserves the bond lengths between the test point and the two end points. 
        inBounds = self.checkPointAgainstList(newPoint, startPoints + endPoints + pointsToAvoid)
        rotNewPoint = newPoint
        numLoops = 0
        # only do the whizzing and checking if we need to!
        while not inBounds and numLoops<self.maxLivesPerNewNode:
            angle = rnd.uniform(0, 2*np.pi)
            rotNewPoint = coords.rotatePointInChain([startPoints[-1], newPoint, endPoints[-1]], angle)
            inBounds = self.checkPointAgainstList(newPoint, startPoints + endPoints + pointsToAvoid)
            numLoops +=1

        return rotNewPoint, inBounds

    def findLastPoint(self, startPoints, endPoints, pointsToAvoid):
        # construct a TNB frame
        TNB = coords.constructTNBFrame(startPoints[-2],
                                       startPoints[-1],
                                       endPoints[-1])
        
        # Compute position in TB plane of last point
        TComponent = np.linalg.norm(startPoints[-1] - endPoints[-1])/2
        BComponent = np.sqrt(self.bondLength**2 - TComponent**2)
        newPoint = startPoints[-1] + BComponent*TNB[2] + TComponent*TNB[0]
 
        # perform a whizz Check       
        return self.whizzCheck(startPoints, endPoints, pointsToAvoid, newPoint)
        
    def findLastTwoPoints(self, startPoints, endPoints, pointsToAvoid):
        
        TNB = coords.constructTNBFrame(startPoints[-2],
                                       startPoints[-1],
                                       endPoints[-1])
        # see page 173 in book 2
        l = np.linalg.norm(startPoints[-1] - endPoints[-1])
        tComponent = np.abs(l-self.bondLength)/2.0
        bComponent = np.sqrt(self.bondLength**2 - tComponent**2)
        
        if l>self.bondLength:
            newPoint1 = startPoints[-1] - tComponent* TNB[0] - bComponent * TNB[2]
        else:
            newPoint1 = startPoints[-1] + tComponent* TNB[0] - bComponent * TNB[2]
        
        newPoint2 = newPoint1 + self.bondLength * TNB[0]

        # whizz check point 1 as if newPoint2 was at the end of endPoints
        newPoint1, inBounds1 = self.whizzCheck(startPoints, endPoints + [newPoint2], pointsToAvoid, newPoint1)
        
        if inBounds1:
            # Now whizz check point 2 as if the new newPoint1 was at the end of startPoints.
            newPoint2, inBounds2 = self.whizzCheck(startPoints + [newPoint1], endPoints, pointsToAvoid, newPoint2)
            
        listsAreGood = False
        if inBounds1 and inBounds2: # only return true if both the new points are golden.
            listsAreGood = True
        
        return newPoint1, newPoint2, listsAreGood 
        
    def findLastTwoPointsFit(self, startPoints, endPoints, distanceBetweenPoints):
        # This function figures out the coordinates of the last 
        # one or two points in a polymer chain such that they are 
        # bondLength distance apart. It stops worrying about bond angles. 
        
        # This is a many dimensional problem and there
        # are many ways to solve it. varies all four angles simultaneously
        # and calculates the resulting distance
        # between the particles. THe function minimises this distance - bondLength.

        # construct TNB frames
        TNBA = coords.constructTNBFrame(startPoints[-3],
                                        startPoints[-2],
                                        startPoints[-1])
        TNBB = coords.constructTNBFrame(endPoints[-3],
                                        endPoints[-2],
                                        endPoints[-1])
        
        # pick bond angles at random from the acceptable range
        betaA = rnd.uniform(self.beta1, self.beta2)
        betaB = rnd.uniform(self.beta1, self.beta2)
        alphaA = rnd.uniform(self.alpha1, self.alpha2)
        alphaB = rnd.uniform(self.alpha1, self.alpha2)

        # optimise the angles for the given function.
        finalAngles = minimize(coords.computeDistBetweenPoints, 
                               [alphaA, alphaB], 
                               args=(startPoints[-1], endPoints[-1], TNBA, TNBB, betaA, betaB, self.bondLength, distanceBetweenPoints))
        
        pointA = startPoints[-1] + self.bondLength * coords.generateTNBVecXYZ(TNBA, betaA, finalAngles['x'][0])
        pointB = endPoints[-1] + self.bondLength * coords.generateTNBVecXYZ(TNBB, betaB, finalAngles['x'][1])

        # true if the optimisation exited successfully (it gave a result!)
        listsAreGood = finalAngles['success']
        
        return pointA, pointB, listsAreGood
   
    def pickRandomPointInDefinedSpace(self, p1, p2, p3):
        # Takes the given three points and picks a new point in XYZ which 
        # is bondLength distance from p3, but inside the given angular ranges
        # as defined by the bond angle p2 - p3 - newPoint and dihedral p1-p2-p3-new point.
        # No need to check againt angle ranges once point is selected because angles 
        # are picked inside those ranges already.
        
        # construct the TNB vectors from the last three points of the nList.
        TNB = coords.constructTNBFrame(p1, p2, p3)        
        
        # compute a directional vector in XYZ coords from the last point on the list to the new point, within the given angular range constraints. 
        newCoordXYZRel = coords.pickRandomTNBDirectionInAngRangeXYZ(TNB, 
                                                                    self.betaMin, 
                                                                    self.betaMax,
                                                                    self.alphaMin,
                                                                    self.alphaMax)
        
        # should be a unit vector already but just in case
        newCoordXYZRelHat = newCoordXYZRel/np.linalg.norm(newCoordXYZRel) 
        
        # scale the vector by the bondlength and add it to the last point in the list to give new position vector
        return p3 + self.bondLength * newCoordXYZRelHat
        
    def getParams(self):
        return self.params 

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
    offset =  np.array([600, 0.0, 0.0])
    connectors = [[2, 1, 0], [numPos - 3, numPos - 2, numPos - 1]]
    backBoneBuildingBlock1 = backboneObject.generateBuildingBlock(numPos, 'CN', connectors) 
    backBoneBuildingBlock1.transformFromBlockFrameToLabFrame(director, startPos, rotation)
    backBoneBuildingBlock1.exportBBK("backbone1")
    backBoneBuildingBlock2 = backboneObject.generateBuildingBlock(numPos, 'NC',  connectors) 
    backBoneBuildingBlock2.transformFromBlockFrameToLabFrame(director, startPos, rotation)
    backBoneBuildingBlock2.exportBBK("backbone2")

    startPoints = backBoneBuildingBlock1.xyzVals[-3:]
    endPoints = backBoneBuildingBlock2.xyzVals[-3:]
    
    # create the generator
    ConstrainedPolymerPackGBB = ConstrainedPolymerPackBBG(filename)
    
    numPoints = 700
    alpha1 = -50
    alpha2 = -70
    beta1 = 100   
    beta2 = 120
    minDist = 2.0
    bondLength = 2.1
       
    # generate the XYZVals in the packed spaced - 
    # this generates a self-avoiding chain between two points that 
    # also avoid a spherical excluded volume region 
    # The bond angles and dihedrals are picked at random from the 
    # a given range of angles and dihedrals.
    ConstrainedPolymerPackBB = ConstrainedPolymerPackGBB.generateBuildingBlock(numPoints, 
                                                                               startPoints, 
                                                                               endPoints, 
                                                                               alpha1, 
                                                                               alpha2, 
                                                                               beta1, 
                                                                               beta2, 
                                                                               minDist, 
                                                                               bondLength)

    ConstrainedPolymerPackGBB.transformFromBlockFrameToLabFrame(director, startPos, rotation)
    ConstrainedPolymerPackGBB.checkBondLengths()
    ConstrainedPolymerPackBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "constrainedPolymer Done"