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
from scipy.optimize import minimize
from Utilities import coordSystems as coords
from Utilities import cartesian as cart
from Utilities import fileIO as fIO

class ConstrainedPolymerPackBBG(BBG):
    '''
    Overloads the BBG object to give a building block with a packing which results in a polymeric chain
    of N monomers with random dihedrals and bond angles, whose ends are constrained by 
    a collection of six points.
    This list of six point defines dihedrals and bond angles for the two end points, and allows a 
    random coil to be inserted easily into an larger structure.  The algorithm defines N points that 
    form a chain with a given bondlength. 
    
    This algorithm defines a non-crossing continuous space curve and threads it with equispaced particles.
    It uses optimisation to ensure the curve is the right length and starts and stops at the desired place. 
    
    Once threaded with the right bondlengths the particle chain is then manipulated using many random crankshaft moves 
    which preserve bond lengths between the particles, rejecting any moves that cause crossovers or 
    violate the envelope conditions.
    
    '''
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        BBG.initialiseParameters(self) 
        
        self.particleName = self.getParam('particleName')
        self.initialModelParams = self.getParam('initialModelParams')
        self.numCrankShaftMoves = self.getParam('numCrankShaftMoves')
        self.stubbornGroupLife = self.getParam('stubbornGroupLife')
        self.stubbornCranks = self.getParam('stubbornCranks')
        #self.numPasses = self.getParam('numPasses')
        #self.maxLivesPerNewNode = self.getParam('maxLivesPerNewNode') 
        #self.scrambleWindow = self.getParam('scrambleWindow')
         
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
        
        # generate the transformation information from building block to pointA and pointB
        self.labDirector, self.labRefPoint, self.labRotation = self.computeTransform()
        
        # convert the pointsToAvoid information from the labFrame to the block frame
        self.pointsToAvoid = self.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, pointsToAvoid)
        
        fIO.saveXYZ([self.pointA, self.pointB], 'S', 'endPoints.xyz')
                     
        
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
        
        fIO.saveXYZ(xyzVals, 'O', 'CP_Stage1_initCurve.xyz')
        
        # clean up Zeros"
        xyzVals = self.cleanUpZeroes(xyzVals)
        
        #print "correcting bondLengths"
        #xyzVals = self.correctBondLengths(xyzVals)
        
        #fIO.saveXYZ(xyzVals, 'O', 'CP_Stage2_evenBondLengths.xyz')

        #print "Scrambling adjacent atoms"
        #xyzVals = self.scrambleAdjacentAtoms(xyzVals, self.numPasses)
        
        #fIO.saveXYZ(xyzVals, 'O', 'CP_Stage3_scrambleAdjacentAtoms.xyz')
        
        print "Randomising the structure with crank shaft moves."
        xyzVals = self.crankShaftMoves(xyzVals, self.numCrankShaftMoves)

        fIO.saveXYZ(xyzVals, 'O', 'CP_Stage4_crankshaft_randomization.xyz')


        print "Ensuring structure is inside envelope"
        xyzVals = self.foldInsideEnvelope(xyzVals)

        fIO.saveXYZ(xyzVals, 'O', 'CP_Stage5_enveloped.xyz')

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

    def correctBondLengths(self, xyzVals):
        
        # whip through the xyz vals and subtly adjust the positions of the atoms
        # so they are all exactly bondlength apart. 
        for i, pos in enumerate(xyzVals):
            if i > 0 and i < len(xyzVals) - 1:
                # if the bondLengths
                b1 = np.linalg.norm(pos - xyzVals[i - 1])
                b2 = np.linalg.norm(pos - xyzVals[i + 1])
                
                if (np.abs(b1 - self.bondLength) > 0.01 ) or (np.abs(b2 - self.bondLength) > 0.01): 
                    tnb = coords.constructTNBFrame(pos, xyzVals[i - 1], xyzVals[i + 1])
                    midPoint = (xyzVals[i - 1] + xyzVals[i + 1])/2
                    l = np.linalg.norm(xyzVals[i - 1] - xyzVals[i + 1])/2
                    r = np.sqrt(self.bondLength**2 - l**2)
                    xyzVals[i] = midPoint + tnb[2] * r
        
        return xyzVals

    def scrambleAdjacentAtoms(self, xyzVals, numPasses):
        # loops through the list whizzing each atom into a new place
        for curPass in range(0, numPasses): 
            print curPass + 1, "passes out of ", numPasses
            for scrambleIndex in range(1, len(xyzVals) - 2, 2):
                print scrambleIndex, "out of", len(xyzVals)-2
                numLoops = 0
                inBounds = False
                while not inBounds and numLoops<self.maxLivesPerNewNode:
                    angle = rnd.uniform(0, 2*np.pi)
                    newPoint = coords.rotatePointInChain(xyzVals[scrambleIndex-1:scrambleIndex+2], angle)
                    # iterate and compare each pair to see if any clash
                    inBounds = True
                    # since we now call this function near the beginning
                    # where the space curve is deterministic and we have a priori knowledge 
                    # about it then we can restrict the scramble check to a window on either side
                    # of the scrambleAtom. 
                    lowerScrambleIndex = scrambleIndex - self.scrambleWindow
                    if lowerScrambleIndex < 0:
                        lowerScrambleIndex = 0
                    
                    upperScrambleIndex = scrambleIndex + self.scrambleWindow
                    if upperScrambleIndex > len(xyzVals) - 1:
                        upperScrambleIndex = len(xyzVals) - 1
                    
                    i=lowerScrambleIndex
                    while inBounds and i < upperScrambleIndex:
                        if np.linalg.norm(newPoint - xyzVals[i]) < self.minDist:
                            inBounds = False
                        i += 1
                    numLoops +=1
                
                if inBounds:
                    xyzVals[scrambleIndex] = newPoint
        
        return xyzVals

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
                                               
    def generateSpaceCurve2(self, curveLength):
        # optimise the parameters to yield for the given curveLength and baseLength.
        finalResult = minimize(self.deltaLengthSpaceCurve, self.initialModelParams, args=(self.baseLength, curveLength) )

        # compute the XYZ points of the curve over uniform domain spacings in parameter. 
        spaceCurveXyzVals = self.spaceCurveDefinition(self.baseLength, self.numPoints, finalResult['x'])
 
        fIO.saveXYZ(spaceCurveXyzVals, 'O', 'spaceCurve.xyz')
        
        # compute the arc length at each point of the space curve.
        arcLengthSource = self.computeArcLengthFromSpaceXYZ(spaceCurveXyzVals)
 
        # compute the desired arc lengths for the particle positions.
        sPos = [ n * self.bondLength for n in range(0, self.numPoints)]
        
        # compute the interpolated xyz points at each desired arc length
        return self.interpXYZFromArcLength(sPos, arcLengthSource, spaceCurveXyzVals)
 
     
    def deltaLengthSpaceCurve(self, params, *args):
        # compute difference between the desired length and the actual 
        # length of a space curve between the points at either end of a 
        # line segment of length baseLength. Centre of baselength is assumed 
        # to be at the origin.
        baseLength = args[0] 
        curveLength = args[1] 
        
        # compute the xyz vals of the space curve with the given parameters
        xyzVals = self.spaceCurveDefinition(baseLength, self.numPoints, params) 
        
        # find the total length of the current space curve by summing sqrt(dx^2 + dy^2 + dz^2) along the curve.
        arcLengthSource = self.computeArcLengthFromSpaceXYZ(xyzVals)
        
        return (curveLength - arcLengthSource[-1])**2          

    def spaceCurveDefinition(self, l, numPoints, *params):
        # Compute the coords along a parametric space curve relative to an origin at the mid point 
        # between the two fixed points p1 and p2 that we're interested in, which are a distance l apart.
        # without loss of generality we can assume that the two points are in a space where the x axis runs between them
        # and each point is at (-l/2, 0, 0) and (l/2, 0, 0). 
        
        # We parameterise the space curve with t such that t=-l/2 corresponds to one fixed point and t=l/2 corresponds to the other point
        # thus t is basically the x coordinate of the space curve. where x runs from p1 to p2.
        # This makes it easier to avoid a self-intersecting space curve with judicious choice of functions.
        # 
        # For each t define a list of r(t), theta(t) and phi(t) values to yield a space curve, and then convert it to xyz coords.  
        #
        # Must ensure that whatever function we put in here:
        # r(-l/2) = l/2 and r(l/2) = l/2.
        # theta(-l/2) = 0 and theta(l/2) = 0
        # phi(-l/2) = -pi/2 and phi(l/2) = pi/2.
        t = np.linspace(-l/2, l/2, num=numPoints)
        
        RMax = params[0][0]
        b = params[0][1]
        
        ry = -4 * RMax * t**2 /(l * l) + RMax    # 0 at +/- l/2
        rList = np.sqrt(t**2 + ry**2)  + b * np.sin(5 * t * 2 * np.pi/l)
        thetaList =  b * np.sin(20 * t * 2 * np.pi/l)  # 0 at t=+/- l/2 for f an integer
        phiList = t * np.pi/l  # -pi/2 at t= -l/2 and pi/2 at t=l/2
        
        return [ np.array([r * np.cos(theta) * np.sin(phi), r * np.cos(theta) * np.cos(phi), r * np.sin(theta)]) for r, theta, phi in zip(rList, thetaList, phiList)]
        
    def computeArcLengthFromSpaceXYZ(self, xyzPoints):
        segLengths = [ np.linalg.norm(a-b) for a,b in zip(xyzPoints[0:-1], xyzPoints[1:])]
        return np.concatenate( ([0], np.cumsum(segLengths)), 0)
    
    def interpXYZFromArcLength(self, arcLengths, arcLengthSource, xyzList ):
        # returns the interpolated xyz positions of a list of points that are defined by an arclength.
        # Finds each sValue (arcLength) in sPosList, and interpolates the x, y, z values
        xyzValsRet = [] # set up output array
        highestVal = False
        indexHigher = 0
        for uniformArcLength in arcLengths:
            while uniformArcLength > arcLengthSource[indexHigher]:
                indexHigher += 1
                if indexHigher==len(arcLengthSource):
                    highestVal = True
                    break
            
            if not highestVal:
                # indexHigher is the index of arcLengthSource 
                # referencing a value of arcLengthSource just above uniformArcLength
                if indexHigher > 0: 
                    deltaArcLengthSource = arcLengthSource[indexHigher] - arcLengthSource[indexHigher - 1] 
                    deltaX = xyzList[indexHigher][0] - xyzList[indexHigher- 1][0]
                    deltaY = xyzList[indexHigher][1] - xyzList[indexHigher- 1][1]
                    deltaZ = xyzList[indexHigher][2] - xyzList[indexHigher- 1][2]
    
                    s = (arcLengthSource[indexHigher] - uniformArcLength)/deltaArcLengthSource
                    
                    newX = xyzList[indexHigher][0] - s * deltaX
                    newY = xyzList[indexHigher][1] - s * deltaY
                    newZ = xyzList[indexHigher][2] - s * deltaZ
                else:
                    newX = xyzList[0][0]
                    newY = xyzList[0][1]
                    newZ = xyzList[0][2]
            else:
                newX = xyzList[-1][0]
                newY = xyzList[-1][1]
                newZ = xyzList[-1][2]
                
            # append the interpolated value of xyz to the output array
            xyzValsRet.append(np.array([newX, newY, newZ]))        
        
        return xyzValsRet 
        
    def getParams(self):
        return self.params 

if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the generator
    ConstrainedPolymerPackGBB = ConstrainedPolymerPackBBG(filename)
    
    # generate backbone realtime parameters
    startPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    numPoints = 200
    baseLine = 10
    minDist = 1.0
    bondLength = 1.5
    
    # generate a curve between the points
    ConstrainedPolymerPackBB = ConstrainedPolymerPackGBB.generateBuildingBlock(numPoints, 
                                                                               baseLine, 
                                                                               minDist, 
                                                                               bondLength,
                                                                               envelope="frustum 11 2 20")
    ConstrainedPolymerPackBB.transformFromBlockFrameToLabFrame(director, startPos, rotation)
    ConstrainedPolymerPackGBB.checkBondLengths()
    ConstrainedPolymerPackBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "constrainedPolymer Done"