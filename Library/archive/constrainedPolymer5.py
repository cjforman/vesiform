'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Builder.BuildingBlock import BuildingBlock
from Library.Ellipsoid import EllipsoidPackNBB as EPBBG
import sys
import time
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
    constrained by two given points, and whose pathway is randomised and then packed 
    randomly into a confined space. 
    
    This algorithm begins by defining a non-crossing continuous space curve which is 
    precisely the right length to space out N equispaced particles. 
    In this case two straight lines out from the startpoint and back to the end point 
    in a triangular arrangement. Takes into account odd or even N.
    
    The complete particle chain is then manipulated using numCrankshaftmoves random 
    "crankshaft moves" which preserve bond lengths between the particles, while rejecting 
    any moves that cause self-intersections or intersections with pointsToAvoid.  
    Acceptable points to use as crankshaft axis are defined in the allowedList
    which is defined in an member function. For example, if you specify every third point,
    then the dihedral angles of every 1 and 2 points will never change, which is useful
    to preserve peptide bond rigidity in a peptide sequence.  The allow list defaults to 
    all points, but can be overridden.
    
    Once randomised, the chain is then packed into the envelope conditions using the same 
    crankshaft move.  This uses a minimisation algorithm which makes individual random 
    crankshaft moves and assesses the resulting structure. Any move which causes self-crossing,
    or intersects with pointsToAvoid is rejected.  Any move reducing the number of points
    outside the envelope is accepted.  Moves that increase the number of points are accepted
    with a probability that decreases exponentially with the increase in number of points.
    For this algorithm to be effective in all cases the allowed List must include the 
    0th and final points, otherwise some particles cannot be moved inside the envelope 
    with crankshaft moves. 
    
    The default envelope is "None". Other possibilities include frustum r1 r2 Z, innersphere
    and outersphere R.   Simply supply a string in the envelope parameter with the instructions
    for the envelope.   You may supply your own envelopes by over riding the checkPointInBounds
    function. Frustum generates a circular frustum whose top radius R1 is centred 
    in the XY plane at the origin (the midpoint of the anchor points). R2 is place at a 
    distance Z immediately below the origin.  The axis of the frustum is aligned with the
    blockDirector (z axis).  Outsphere defines a sphere of radius R, centered on the 
    building block origin into which the points must be packed. innersphere defines a 
    spherical region of diameter 0.9 * baseLength centered between the anchor points which must 
    be avoided. Halfspace prevents any points from being placed on either side of a plane. 
    
    If the numSteps taken by the packing algorithm exceeeds maxNumFoldingMoves then the 
    algorithm terminates with points left outside the envelope, and issues a warning.   
    
    The two anchoring points are supplied as any two points in the lab frame, and these are 
    transformed into (actually define!) the "building block frame". The two block anchor 
    points are placed the same distance apart as the user supplied points, but arranged 
    symmetrically about the origin on the X axis. The director of the building block 
    is then taken as the Z axis. Once the final chain is complete the building block 
    is created and transformed back into the lab frame so that the end points of the 
    chain match up with the supplied lab pointsA and B. The building block reference 
    point and director are then matched to the lab orientation and position of the 
    building block, allowing further transformations to the polymer chain as necessary.
    
    A list of externally supplied lab frame points are automatically transformed into the 
    building block frame using the internal transformation. This allows the user to supply
    a list of points that the chain must not be allowed to intersect with. Such a list increases
    the time taken to find a suitable conformation.
    
     
    '''
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        BBG.initialiseParameters(self) 
        
        self.particleName = self.getParam('particleName')
        self.foldingTemp = self.getParam('foldingTemp')
        self.maxNumFoldingMoves = self.getParam('maxNumFoldingMoves')
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for constrained polymer object"
            sys.exit()

    def generateBuildingBlock( self, 
                               numPoints,
                               pointA,
                               pointB, 
                               minDist,
                               bondLength,
                               numCrankMoves,
                               pointsToAvoid=[],
                               envelopeList=["None"]):

        self.numPoints = numPoints
        self.labPointA = pointA
        self.labPointB = pointB
        self.baseLength = np.linalg.norm(pointA - pointB)
        self.bondLength = float(bondLength)
        self.minDist = minDist
        self.blockPointA = np.array([-self.baseLength/2, 0.0, 0.0])
        self.blockPointB = np.array([self.baseLength/2, 0.0, 0.0])
        self.numCrankMoves = numCrankMoves
        self.allowedList = self.generateAllowedList()
        
        if self.dumpInterimFiles:
            fIO.saveXYZList([self.blockPointA, self.blockPointB], ['Ca', 'O'], 'blockPoints.xyz')
        
        # generate the BuildingBlock reference point earlier than usual because 
        # we need the transformation for the pointsToAvoid input.
        self.blockRefPoint = self.generateBuildingBlockRefPoint()
        
        # generate the BuildingBlock director unit vector earlier than usual because 
        # we need the transformation for the pointsToAvoid input.
        self.blockDirectorHat = self.generateBuildingBlockDirector()
        
        # generate the transformation information from building block to labPointA and labPointB
        self.labDirector, self.labRefPoint, self.labRotation = self.computeTransform()
        
        # convert the pointsToAvoid information from the labFrame to the block frame
        blockPointsToAvoid = coords.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, self.blockDirectorHat, self.blockRefPoint, pointsToAvoid)
        
        if self.dumpInterimFiles == 1:
            # these are debugging tests to make sure the transform is correct
            blockPointATrans = coords.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, self.blockDirectorHat, self.blockRefPoint, [pointA])[0]
            blockPointBTrans = coords.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, self.blockDirectorHat, self.blockRefPoint, [pointB])[0]       
        
            fIO.saveXYZList([blockPointATrans, blockPointBTrans], ['Ca', 'O'], "blockPointsTransFromLab.xyz")
            fIO.saveXYZ(blockPointsToAvoid, 'Li', "blockPointsToAvoid.xyz")
            labPointsToAvoidTrans = coords.transformFromBlockFrameToLabFrame(self.labDirector, self.labRefPoint, self.labRotation, self.blockDirectorHat, self.blockRefPoint, blockPointsToAvoid)
            fIO.saveXYZ(labPointsToAvoidTrans , 'Be', "labPointsToAvoidTrans.xyz")
        
        # parse the envelope list now to check pointA and point B
        self.parseEnvelopeList(envelopeList)
        
        # store points to avoid to check pointA and point B
        self.pointsToAvoid = blockPointsToAvoid
        
        # check starting points are legal or it's gonna be a long wait.
        if not self.checkPointInBounds(self.blockPointA):
            print "Error Warning: PointA out of bounds"
            time.sleep(3)
             
        if not self.checkPointInBounds(self.blockPointB):
            print "Error Warning: PointB out of bounds"
            time.sleep(3)
        
        return BBG.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=blockPointsToAvoid)

    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0]) 
    
    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints

    def checkBondLengths(self):
        print [ np.linalg.norm(self.buildingBlockXYZ[i] - self.buildingBlockXYZ[i-1]) for i in range(1, len(self.buildingBlockXYZ)) ]       
        
    def generateAllowedList(self):
        return [i for i in range(0, self.numPoints) ]
        
    def generateBuildingBlockXYZ(self):
        print "Generating initial conformation in blockspace."
        xyzVals = self.generateSpaceCurve()
        
        if self.dumpInterimFiles==1:
            self.blockNames = self.generateBuildingBlockNames()
            fIO.saveXYZList(xyzVals, self.blockNames, 'chainConnectedPointB.xyz')
        
        # clean up Zeros
        xyzVals = self.cleanUpZeroes(xyzVals)
        
        print "Randomising chain using crankshaft moves."
        xyzVals = self.crankShaftMoves(xyzVals, self.numCrankMoves, 1)


        if self.dumpInterimFiles==1 and self.numCrankMoves > 0:
            fIO.saveXYZList(xyzVals, self.blockNames, 'crankedChain.xyz')
            print "Visualising envelope"
            self.visualiseEnvelope(10000, 20.0, 20.0, 20.0)
        
        print "Ensuring structure is inside envelope"
        xyzVals = self.foldInsideEnvelope(xyzVals)

        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'foldedChain.xyz')

        return xyzVals
    
    def getBuildingBlock(self):
        # overload access function that returns a building block based on the backbone
        BB = BuildingBlock( self.buildingBlockXYZ, 
                            self.blockNames, 
                            self.blockConnectors, 
                            self.blockRefPoint,
                            self.blockDirectorHat)
        # perform the transformation back to the pointA, pointB perspective
        BB.transformBBToLabFrame(self.labDirector, self.labRefPoint, self.labRotation)
        return BB

    def computeTransform(self):
        # computes the necessary director, refPoint and rot to successfully map from 
        # the buildingBlockRefFrame to the correct pointA and pointB.
        midPoint = (self.labPointA + self.labPointB)/2
        baseVector = (midPoint - self.labPointA )
        baseVectorHat = baseVector/np.linalg.norm(baseVector)
        NormVec = np.cross(np.array([0.0, 0.0, 1.0]), baseVectorHat)
        NormVecHat = NormVec/np.linalg.norm(NormVec)
        DirectorHat = np.cross(baseVectorHat, NormVecHat) 
        rot = np.arctan2(baseVectorHat[1], baseVectorHat[0])
        return (DirectorHat, midPoint, rot)

    def cleanUpZeroes(self, xyzVals):
        return np.around(xyzVals, 5)
        
    def foldInsideEnvelope(self, xyzVals):
        
        # Perform a random crank shaft using the allowed list. If there are points intersecting reject it out of hand.
        # if there are less points outside the zone than previously accept it.
        # if there are more points outside the zone then accept it with a probability that depends
        # exponentially on the difference between the current minimum number of points inside.
        
        # check which indices are outside the array
        curIndicesOutside = self.checkEnvelope(xyzVals)
        curXYZVals = xyzVals
        minXYZVals = xyzVals
        curNumIndicesOutside = len(curIndicesOutside)
        minNumIndicesOutside = curNumIndicesOutside  
        minIndices = curIndicesOutside 
        maxStepScale = 1.0
        numMoves = 0
        
        while minNumIndicesOutside > 0 and numMoves < self.maxNumFoldingMoves:
            
            # Do the crank shaft move on the current working set. 
            # Already rejects crossovers with pointsToAvoid and selfcrossovers.
            # Will return the current set in those cases so no move on number of points.   
            newXYZVals = self.crankShaftMoves(curXYZVals, 1, maxStepScale)

            # find indices outside the envelope and count them.        
            newIndicesOutside = self.checkEnvelope(newXYZVals)
            newNumIndicesOutside = len(newIndicesOutside)
            
        
        
            # if a new global minimum number of indices inside the envelope then keep the move.
            if newNumIndicesOutside < minNumIndicesOutside:
                # keep a permanent log of the the new best set.
                minXYZVals = newXYZVals[:]
                minNumIndicesOutside = newNumIndicesOutside
                minIndices = newIndicesOutside
                
                # update the working copy
                curXYZVals = newXYZVals[:]
                curNumIndicesOutside = newNumIndicesOutside
                
                print "step: ", numMoves, " minNumIndicesOutside: ", minNumIndicesOutside, "curNumIndicesOutside:", curNumIndicesOutside
                 
            # if the new move means that the num indices outside has gone up above the curXYZVals then 
            # roll the die to decide whether or not to replace the curXYZVals
            if newNumIndicesOutside >= curNumIndicesOutside:
                # roll the die:
                if rnd.uniform(0.0, 1.0) < np.exp( (minNumIndicesOutside - newNumIndicesOutside) / self.foldingTemp ): 
                    # update the working copy of coords but don't replace the current global minimum
                    curXYZVals = newXYZVals
                    curNumIndicesOutside = newNumIndicesOutside 
            
            numMoves += 1

            if numMoves % 10==0:                 
                print "step: ", numMoves, " minNumIndicesOutside: ", minNumIndicesOutside, "curNumIndicesOutside:", curNumIndicesOutside
        
        if minNumIndicesOutside > 0:
            print "Warning: there are points outside the envelope that were not moved inside."
            if self.dumpInterimFiles:
                fIO.saveXYZ(minXYZVals[minIndices], 'B', 'outsideEnvelope.xyz')
                
                
        return minXYZVals
        
    def crankShaftMoves(self, xyzValsOrig, numCrankShaftMoves, maxStepScale):
        # picks two points at random and rotates all the points between them by a random angle
        # about the axis between the points.
        #
        # The angle scale is defined by maxStepScale.
        #
        # If the picked point is not on the allowed list then pick another point.
        #
        # Checks to see if any of the balls in the new positions 
        # are closer than minDist to any of the other balls. If they are then it rejects the move.
        # this kind of move preserves the bond lengths that we went to so much trouble to get right.

        workingXYZVals = xyzValsOrig[:]

        numMoves = 0 
        while numMoves < numCrankShaftMoves:

            # first choose two indices
            index1=0
            index2=0
            # check for equal or adjacent indices and indices on the allowed list
            while (abs(index1-index2) < 2) or not (index1 in self.allowedList) or not (index2 in self.allowedList):
                index1 = rnd.randint(0, len(workingXYZVals))
                index2 = rnd.randint(0, len(workingXYZVals))

            indexMin = min([index1, index2])
            indexMax = max([index1, index2])

            # performs a single crank on the given indices. 
            # If move is rejected due to overlaps then the original array is returned
            workingXYZVals = self.crankShaftMove(workingXYZVals, indexMin, indexMax, maxStepScale)

            numMoves += 1
            if numMoves % 20 == 0:
                print numMoves, " out of ", numCrankShaftMoves 
       
        return workingXYZVals

    def crankShaftMove(self, workingXYZVals, indexMin, indexMax, maxStepScale):
        # Performs rotation of random size about the axis between the atoms
        # defined by the two indices.  
        # This is known as a crank rotation which preserves the bond lengths.
        
        # create three working sublists
        lowList = workingXYZVals[0:indexMin]
        workList = workingXYZVals[indexMin:indexMax + 1]
        highList = workingXYZVals[indexMax + 1:]
            
        # get the angle to move.
        angle = rnd.uniform(-maxStepScale * np.pi, maxStepScale * np.pi)
            
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
    numPoints = 100
    pointA = 11.0 * np.array([1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 1.0/np.sqrt(3)])
    pointB = 11.0 * np.array([-1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), -1.0/np.sqrt(3)])
    minDist = 1.0
    bondLength = 1.5
    crankMoves = 20
    spherePointGenerator = EPBBG('../../Library/EllipsoidPacking.txt')
    spherePoints = spherePointGenerator.generateBuildingBlock(30, 11, 11, 11, -90, 90, -180, 180, 4)
    pointsToAvoid = spherePoints.getAtomsXYZ()
    
    fIO.saveXYZList([pointA, pointB], ['Ca', 'O'], 'labPoints.xyz')
    fIO.saveXYZ(pointsToAvoid, 'Na', 'labPointsToAvoid.xyz')
    
    # generate a curve between the speicifed points
    ConstrainedPolymerPackBB = ConstrainedPolymerPackGBB.generateBuildingBlock(numPoints, 
                                                                               pointA,
                                                                               pointB, 
                                                                               minDist, 
                                                                               bondLength,
                                                                               crankMoves,
                                                                               envelopeList=["outersphere 12.0", "innersphere 10.0"],
                                                                               pointsToAvoid=pointsToAvoid)
    ConstrainedPolymerPackGBB.checkBondLengths()
    ConstrainedPolymerPackBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "constrainedPolymer Done"