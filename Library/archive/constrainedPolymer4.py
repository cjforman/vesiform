'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Builder.BuildingBlock import BuildingBlock
import sys
import numpy as np
import random as rnd
import itertools as it
from Utilities import coordSystems as coords
from Utilities import cartesian as cart
from Utilities import fileIO as fIO

class ConstrainedPolymerPackBBG(BBG):
    '''
    Overloads the building block generator object which, when called at run time will 
    generate a building block containing a chain of N monomers whose ends are 
    constrained by two given points, and whose pathway is random. 
    
    This algorithm begins by defining a non-crossing continuous space curve which is 
    precisely the right length to space out N equispaced particles. 
    In this case two straight lines out from the startpoint and back to the end point 
    in a triangular arrangement. (takes into account odd or even N).
    
    The complete particle chain is then manipulated using random crankshaft moves 
    which preserve bond lengths between the particles, and rejecting any moves that cause 
    self-intersections.
    
    Once randomised, the chain is then packed into the envelope conditions using the same 
    crank shaft move but only on particles outside the envelope, using axes just inside 
    the envelope. The default envelope is "None". This process also avoids a list of 
    externally supplied points that are transformed into the building block frame using 
    the transformation computed using the initial pointA and pointB. To invoke the latter
    simply supply a list of numpy points using the pointsToAvoid parameter, which defaults
    to and empty array.
    
    The construction is performed in a building block frame in which two points are
    space the apposite distance apart (baselength) along the X axis. The director of the building
    block is the Z axis.  Once the final chain is complete the building block is created and transformed
    back into the lab frame so that the end points of the chain match up with the supplied pointsA and B.
    
    The building block reference point and director are then matched to the new orientation and position
    of the building block, allowing further transformations to the polymer chain as necessary. 
    '''
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        BBG.initialiseParameters(self) 
        
        self.particleName = self.getParam('particleName')
        self.numCrankShaftMoves = self.getParam('numCrankShaftMoves')
        self.stubbornGroupLife = self.getParam('stubbornGroupLife')
        self.stubbornCranks = self.getParam('stubbornCranks')
         
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for constrained polymer object"
            sys.exit()

    def generateBuildingBlock( self, 
                               numPoints,
                               pointA,
                               pointB, 
                               minDist,
                               bondLength,
                               pointsToAvoid=[],
                               envelope="None"):

        self.numPoints = numPoints
        self.pointA = pointA
        self.pointB = pointB
        self.baseLength = np.linalg.norm(pointA - pointB)
        self.bondLength = float(bondLength)
        self.minDist = minDist
        self.pointA = np.array([-self.baseLength/2, 0.0, 0.0])
        self.pointB = np.array([self.baseLength/2, 0.0, 0.0])
        self.pointsToAvoid = pointsToAvoid
        self.envelope = envelope.split()
        
        # generate the BuildingBlock reference point earlier than usual because 
        # we need the transformation for the pointsToAvoid input.
        self.blockRefPoint = self.generateBuildingBlockRefPoint()
        
        # generate the BuildingBlock director unit vector earlier than usual because 
        # we need the transformation for the pointsToAvoid input.
        self.blockDirectorHat = self.generateBuildingBlockDirector()
        
        # generate the transformation information from building block to pointA and pointB
        self.labDirector, self.labRefPoint, self.labRotation = self.computeTransform()
        
        # convert the pointsToAvoid information from the labFrame to the block frame
        self.pointsToAvoid = self.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, pointsToAvoid)
        
        return BBG.generateBuildingBlock(self, numPoints)

    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0]) 
    
    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints

    def checkBondLengths(self):
        print [ np.linalg.norm(self.buildingBlockXYZ[i] - self.buildingBlockXYZ[i-1]) for i in range(1, len(self.buildingBlockXYZ)) ]       
        
    def generateBuildingBlockXYZ(self):
        print "Generating initial conformation in blockspace."
        xyzVals = self.generateSpaceCurve()
        
        # clean up Zeros
        xyzVals = self.cleanUpZeroes(xyzVals)
        
        print "Randomising the structure with crank shaft moves."
        xyzVals = self.crankShaftMoves(xyzVals, self.numCrankShaftMoves)

        print "Ensuring structure is inside envelope"
        xyzVals = self.foldInsideEnvelope(xyzVals)

        return xyzVals
    
    def getBuildingBlock(self):
        # overload access function that returns a building block based on the backbone
        BB = BuildingBlock( self.buildingBlockXYZ, 
                            self.blockNames, 
                            self.blockConnectors, 
                            self.blockRefPoint,
                            self.blockDirectorHat)
        # perform the transformation back to the pointA, pointB perspective
        BB.transformFromBlockFrameToLabFrame(self.labDirector, self.labRefPoint, self.labRotation)
        return BB

    def computeTransform(self):
        # computes the necessary director, refPoint and rot to successfully map from the buildingBlockRefFrame to the
        # correct pointA and pointB.
        midPoint = (self.pointA + self.pointB)/2
        baseVector = (self.pointA - midPoint)
        baseVectorHat = baseVector/np.linalg.norm(baseVector)
        NormVec = np.cross(np.array([0.0, 0.0, 1.0]), baseVectorHat)
        NormVecHat = NormVec/np.linalg.norm(NormVec)
        DirectorHat = np.cross(baseVectorHat, NormVecHat) 
        rot = np.arctan2(baseVectorHat[1], baseVectorHat[0])
        return (DirectorHat, midPoint, rot)

    def cleanUpZeroes(self, xyzVals):
        return np.around(xyzVals, 5)

    def checkEnvelope(self, XYZVals):
        # runs through the list checking to see if a given point is outside or outside the envelope bounds
        # if so then make a note of the index
        indicesOutside = []
        for n, pos in enumerate(XYZVals):
            if not self.checkPointInBounds(pos):
                indicesOutside.append(n)
        return indicesOutside

    def checkPointInBounds(self, pos):
        # no constraints in this envelope
        inBounds = True

        # check to see if pos is clashing with points to avoid 
        for pta in self.pointsToAvoid:
            if np.linalg.norm(pta - pos) < self.minDist:
                print "poinstToAvoid rejection"
                inBounds = False
        
        if inBounds:
            # in this envelope any point outside the sphere is out of bounds
            if self.envelope[0] == "outersphere":
                # sphere centered at origin with input radius
                inBounds = True
                if np.linalg.norm(pos) > self.envelope[1]:
                    inBounds=False

        if inBounds:        
            # in this envelope any point inside the inner sphere is out of bounds
            if self.envelope[0] == "innersphere":
                # sphere is centered at origin.
                inBounds = True
                if np.linalg.norm(pos) < (self.baseLength - self.minDist * 2.0)/2.0:
                    inBounds=False
        
        if inBounds:
            # Envelope is frustum with radius2 at z = 0 and radius1 at z = height.
            if self.envelope[0]=='frustum':
                Radius2 = float(self.envelope[1])
                Radius1 = float(self.envelope[2])
                Height = float(self.envelope[3])
                
                # principal axis of envelope is aligned with Z axis. 
                # Compute radial distance from z axis and z height
                zComponent = pos[2] 
                r = np.sqrt(pos[0]**2 + pos[1]**2)
        
                if np.abs(zComponent - 0) < 1e-10:
                    zComponent = 0
        
                inBounds = True
                if  zComponent < 0 or zComponent > Height:
                    inBounds = False
                
                if r > coords.FrustumRadiusAtGivenZ(zComponent, 0, Height, Radius2, Radius1):
                    inBounds = False 
        
        return inBounds
        
    def findContinuousIndexBlock(self, indexList, blocksToIgnore):
        # Function returns the first block of indices that is not identical to one in the blocksToIgnore list
        # if no such block is found returns a None 

        blockNotFound = True # assume we will not find a valid block of indices        
        currentBlock = [ indexList[0] ]
        
        # loop through the index list and return the first block of continuous indices that is not in blocksToIgnore
        curPoint = 1
        while blockNotFound:
        
            # if the current index is precisely 1 away from the last value in the output
            # then add current index to the output array and increment counter 
            try:
                continuousBlock = True
                while continuousBlock and (curPoint<len(indexList)):
                    if (indexList[curPoint] - currentBlock[-1])==1: 
                        currentBlock.append(indexList[curPoint])
                    else:
                        continuousBlock = False
                    curPoint += 1
            except IndexError:
                print "Index Error" # this placed here to catch a mystery debugging error that we encountered once.
    
            # check to see if the current Block is in the blocks to ignore list
            if currentBlock in blocksToIgnore: 
                # if so then check to see if we are at the end of the list
                if curPoint<len(indexList):
                    # if not at the end of the list start a new block
                    currentBlock = [indexList[curPoint]]
                    curPoint += 1
                else:
                    # we are at end of list and all previous blocks have been ignored.
                    # Game over.
                    break
            else:
                # we found a block
                blockNotFound = False
            
        return currentBlock, not blockNotFound
        
    def foldInsideEnvelope(self, xyzVals):
        
        workingXYZVals = xyzVals[:] 
        
        # initialise loop flow control variables
        numAttemptsToFold = 0
        continuousIndexBlock = []
        BlocksToIgnore = []
        keepGoing = True # assume there is something to look at

        # loop while there are blocks of indices outside the envelope that are not being ignored        
        while keepGoing:
            
            # check which indices are outside the working array
            indicesOutside = self.checkEnvelope(workingXYZVals)

            # if there are no indices outside then we're done so bail!
            if len(indicesOutside) == 0:
                keepGoing = False
            else:
                print "num indices outside:", len(indicesOutside), "out of", len(xyzVals) 

                # remember the old index block
                oldContinuousIndexBlock = continuousIndexBlock[:] 

                # function returns the next continuous index block. 
                # If there are no more blocks available to work on then set keepGoing to False 
                continuousIndexBlock, keepGoing = self.findContinuousIndexBlock(indicesOutside, BlocksToIgnore)
                
                # check to see if we are looking at the same block as last time
                if continuousIndexBlock == oldContinuousIndexBlock:
                    # increment the counter for how many times we've looked at this block
                    numAttemptsToFold +=1
                    
                    # check to see if we have tried to fold this more than stubbornGroupLife times in a row 
                    if numAttemptsToFold >= self.stubbornGroupLife:
                        # we have tried to fold this group stubbornGroupLife or more times in a row
                        
                        print "Found a stubborn group: Peforming stubbornCranks crankShaft Moves"
                        workingXYZVals = self.crankShaftMoves(workingXYZVals, self.stubbornCranks)
                        
                        # we have a new workingXYZVals so check the indices
                        indicesOutside = self.checkEnvelope(workingXYZVals)

                        # if there are still indices outside then process them
                        if len(indicesOutside) > 0:
                            print "num indices outside:", len(indicesOutside), "out of", len(xyzVals) 
                            
                            # find the next continuous IndexBlock
                            continuousIndexBlock, keepGoing = self.findContinuousIndexBlock(indicesOutside, BlocksToIgnore)
                                
                            # check to see if we still have the stubborn group after all these crank shafts.
                            if continuousIndexBlock == oldContinuousIndexBlock:
                                numAttemptsToFold = 0 # reset the counter
                                print "Adding stubborn block to ignore list: ", continuousIndexBlock
                                BlocksToIgnore.append(continuousIndexBlock[:])
                            else:
                                # new block so reset counter
                                numAttemptsToFold = 0
                        else:
                            keepGoing = False
                else:
                    # It's a new block so reset the counter
                    numAttemptsToFold = 0 # reset counter
                
                # attempt to fold the current continuous Index Block
                
                # find the folding axis atoms - the ones just inside the envelope.
                indexMin = continuousIndexBlock[0] - 1
                indexMax = continuousIndexBlock[-1] + 1
    
                # crank those indices outside the envelope by a random amount
                # about an axis that is just inside the envelope. (rejects automatically if it caused an overlap) 
                workingXYZVals = self.crankShaftMove(workingXYZVals, indexMin, indexMax)

        if len(indicesOutside)>0:
            print "Warning: there are points outside the envelope that cannot be moved inside."
      
        return workingXYZVals
        
    def crankShaftMoves(self, xyzValsOrig, numCrankShaftMoves):
        # picks two points at random and rotates all the points between them by a random angle
        # about the axis between the points. Checks to see if any of the balls in the new positions 
        # are closer than minDist to any of the other balls. If they are then it rejects the move.
        # this kind of move preserves the bond lengths that we went to so much trouble to get right.

        workingXYZVals = xyzValsOrig[:]

        numMoves = 0 
        while numMoves < numCrankShaftMoves:

            # first choose two indices
            index1=0
            index2=0
            # check for equal or adjacent indices
            while abs(index1-index2) < 2:
                index1 = rnd.randint(0, len(workingXYZVals))
                index2 = rnd.randint(0, len(workingXYZVals))

            indexMin = min([index1, index2])
            indexMax = max([index1, index2])

            # performs a single crank on the given indices. 
            # If move is rejected due to overlaps then the original array is returned
            workingXYZVals = self.crankShaftMove(workingXYZVals, indexMin, indexMax)

            numMoves += 1
            if numMoves % 20 == 0:
                print numMoves, " out of ", numCrankShaftMoves 
       
        return workingXYZVals

    def crankShaftMove(self, workingXYZVals, indexMin, indexMax):
        # Performs rotation of random size about the axis between the atoms
        # defined by the two indices.  
        # This is known as a crank rotation which preserves the bond lengths.
        
        # create three working sublists
        lowList = workingXYZVals[0:indexMin]
        workList = workingXYZVals[indexMin:indexMax + 1]
        highList = workingXYZVals[indexMax + 1:]
            
        # get the angle to move.
        angle = rnd.uniform(0, 2*np.pi)
            
        # compute the rotation axis of the particles
        n = workList[-1] - workList[0]
        nHat = n/np.linalg.norm(n)
        newVals = [ cart.rotPAboutAxisAtPoint(p, workList[0], nHat, angle) for p in workList]
            
        validCrank = True # assume success
        # generate iterator over the product pairs with the new vals in the lower half of the list 
        lowListPairs = it.product(newVals, lowList)
        # iterate and compare each pair to see if any match
        for pair in lowListPairs:
            if np.linalg.norm(pair[0]-pair[1]) < self.minDist:
                validCrank = False
                break 
            
        if validCrank:
            # generate iterator over the product pairs with the new vals in the upper half of the list 
            highListPairs = it.product(newVals, highList)
            # iterate and compare each pair to see if any match
            for pair in highListPairs:
                if np.linalg.norm(pair[0]-pair[1]) < self.minDist:
                    validCrank = False
                    break
                    
        if validCrank:
            # generate iterator over the product pairs with the new vals and the pointsToAvoid list 
            avoidList = it.product(newVals, self.pointsToAvoid)
            # iterate and compare each pair to see if any match
            for pair in avoidList:
                if np.linalg.norm(pair[0]-pair[1]) < self.minDist:
                    validCrank = False
                    break
        
        # reconstruct workingXYZVals
        if validCrank:
            workingXYZVals = np.concatenate((lowList, newVals, highList), 0)
               
        # if valid Crank is false then the original list is returned intact
        # if valid crank is true then the new vals are inserted into the new list
        return workingXYZVals
    
    def generateSpaceCurve(self):
        # Computing a straight conformation which is just a triangle in a plane.
        # Makes the geometry simple.
        
        b = self.bondLength
        n = self.numPoints
        
        if (n % 2)==0: 
            # n is even
            nSide = n/2

            # n is odd
            nSide = (n + 1) / 2
   
            # compute triangle sides         
            TriS = b * (float(nSide) - 1.0) 
            TriH  = np.sqrt( TriS**2 - ( ( self.baseLength - b ) / 2.0 )**2 ) 
            TriG = TriH/( (self.baseLength - b) / 2.0 )
            TriHA = TriG * self.baseLength/2.0  # intercept on vertical axis not the same yVal of heighest point in array.
            
            # compute coordinates
            xVals1 = np.linspace(-self.baseLength/2.0, -b/2, num=nSide)
            xVals2 = np.linspace( b/2, self.baseLength/2.0, num=nSide)
            xVals = np.concatenate( (xVals1, xVals2), 0)
            yVals = [ x * TriG + TriHA if x < 0 else - x * TriG + TriHA for x in xVals]
        else:
            # n is odd
            nSide = (n + 1) / 2
   
            # compute triangle sides         
            TriS = b * (float(nSide) - 1.0) 
            TriH  = np.sqrt(TriS**2 - (self.baseLength/2.0)**2) 
            TriG = TriH/(self.baseLength/2.0)
            
            # compute coordinates
            xVals = np.linspace(-self.baseLength/2.0, self.baseLength/2, num=n)
            yVals = [ x * TriG + TriH if x < 0 else - x * TriG + TriH for x in xVals] 
        
        # output values as numpy arrays in 3D
        return [ np.array([xVal, yVal, 0.0]) for xVal, yVal in zip(xVals, yVals) ]
        
    def getParams(self):
        return self.params 

if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the generator
    ConstrainedPolymerPackGBB = ConstrainedPolymerPackBBG(filename)
    
    # generate backbone realtime parameters
    numPoints = 200
    pointA = np.array([10, 1, 1])
    pointB = np.array([-10, -1, -1])
    minDist = 1.0
    bondLength = 1.5
    pointsToAvoid = [ np.array([i, 0, 0]) for i in range(-10, 10)]
    
    fIO.saveXYZ([pointA, pointB], 'S', 'endPoints.xyz')
    fIO.saveXYZ(pointsToAvoid, 'K', 'pointsToAvoid.xyz')
    
    # generate a curve between the speicifed points
    ConstrainedPolymerPackBB = ConstrainedPolymerPackGBB.generateBuildingBlock(numPoints, 
                                                                               pointA,
                                                                               pointB, 
                                                                               minDist, 
                                                                               bondLength,
                                                                               envelope="frustum 11 2 20",
                                                                               pointsToAvoid=pointsToAvoid)
    ConstrainedPolymerPackGBB.checkBondLengths()
    ConstrainedPolymerPackBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "constrainedPolymer Done"