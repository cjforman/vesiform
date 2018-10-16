'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
import sys
import time
import numpy as np
import random as rnd
from Utilities import coordSystems as coords
from Utilities import cartesian as cart
from Utilities import fileIO as fIO

class RandomPolymerPackBBG(BBG):
    '''
    Overloads the building block generator class to allow construction of a polymer chain 
    from a single point, whose pathway is randomised and packed into a confined space, while
    avoiding certain points and ensuring that the bond lengths are preserved between each 
    particle in the chain. 
    
    An intial arbitrary chain is created starting at one of the points and stretching out
    # which has well defined bond angles and dihedrals.
    
    The particle chain is then manipulated by picking a random bond on the chain and 
    rotating the entire free end of the chain about that bond. Such a move preserves bond lengths 
    between the particles and any move that causes self-intersections or intersections with 
    pointsToAvoid is rejected.  
    
    A list of acceptable points to use as crankshaft axis points are defined in the allowedList
    which is defined in a member function. For example, if you specify every third point,
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
                               alpha1,
                               alpha2,
                               beta1,
                               beta2, 
                               minDist,
                               bondLength,
                               pointsToAvoid=[],
                               visualiseEnvelope=(0, 20, 'envelope.xyz'),
                               envelopeList=["None"],
                               showBlockDirector=False):

        self.numPoints = numPoints
        self.pointA = pointA
        self.bondLength = float(bondLength)
        self.minDist = minDist
        self.blockDirectorHat = np.array([0.0, 0.0, 1.0])
        self.allowedList = self.generateAllowedList()
        self.blockNames = self.generateBuildingBlockNames()
        self.parseEnvelopeList(envelopeList)
        self.alpha1 = min(alpha1*np.pi/180, alpha2*np.pi/180)
        self.alpha2 = max(alpha1*np.pi/180, alpha2*np.pi/180)
        self.beta1 = min(beta1*np.pi/180, beta2*np.pi/180)
        self.beta2 = max(beta1*np.pi/180, beta2*np.pi/180)
        self.pointsToAvoid = pointsToAvoid

        # check starting points are legal or it's gonna be a long wait.
        if not self.checkPointInBounds(self.pointA):
            print "Error Warning: PointA out of bounds"
            time.sleep(3)
        
        return BBG.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, visualiseEnvelope=visualiseEnvelope, pointsToAvoid=pointsToAvoid)

    def generateBuildingBlockDirector(self):
        director = self.buildingBlockXYZ[-1] - self.buildingBlockXYZ[0]
        return  director/np.linalg.norm(director)
    
    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints

    def checkBondLengths(self):
        print [ np.linalg.norm(self.buildingBlockXYZ[i] - self.buildingBlockXYZ[i-1]) for i in range(1, len(self.buildingBlockXYZ)) ]       
        
    def generateAllowedList(self):
        # increase probability of selecting points lower down the list
        allowed=[]
        numTimes = self.numPoints
        for i in range(0, self.numPoints):
            if i==0:
                allowed = [i] * numTimes
            else:
                allowed = allowed + [i] * numTimes
            numTimes -= 1
        
        return allowed
        
    def generateBuildingBlockXYZ(self):
        xyzVals = self.generateSpaceCurve()
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'initialChain.xyz')

        #print "Ensuring structure is inside envelope"
        xyzVals = self.foldInsideEnvelope(xyzVals)

        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'foldedChain.xyz')

        return xyzVals
    
    def generateSpaceCurve(self):
        # generates a space curve with specified bondLength but with bond angles and dihedrals 
        # chosen from a range
        
        b = self.bondLength
        n = self.numPoints
        
        # generate 3 points a distance b apart one of which is pointA
        s2 = self.pointA - b * self.blockDirectorHat
        s1Dir = np.array([rnd.uniform(0.0, 1.0), rnd.uniform(0.0, 1.0), rnd.uniform(0.0, 1.0)])
        s1 = s2 - b * s1Dir/np.linalg.norm(s1Dir)
        spaceCurve = [s1, s2, self.pointA ]
        
        # have first point so add a further n-1 points
        for _ in range(n - 1):
            # construct TNB frame from last three points of the space curve 
            tnb = coords.constructTNBFrame(spaceCurve[-3], spaceCurve[-2], spaceCurve[-1])

            beta = rnd.uniform(self.beta1, self.beta2)
            alpha = rnd.uniform(self.alpha1, self.alpha2)
            
            # compute a new direction based on beta and alpha
            dirn = coords.generateTNBVecXYZ(tnb, beta, alpha)
            
            # construct next space curve point 
            spaceCurve.append(spaceCurve[-1] + b * dirn/np.linalg.norm(dirn))
        
        currentAxis = coords.axisFromHelix(spaceCurve[2:])
        currentRefPoint = self.pointA
        
        spaceCurve = coords.transformFromBlockFrameToLabFrame(self.blockDirectorHat, self.pointA, 0.0, currentAxis, currentRefPoint, spaceCurve[2:])
        
        return spaceCurve
    
    def cleanUpZeroes(self, xyzVals):
        return np.around(xyzVals, 5)

    def dihedralTwist(self, xyzVals, maxStepRange):  
        
        # can only do this if there are sufficient atoms in the array
        if len(xyzVals) > 3:
            # pick a random point in the allowedList
            index = rnd.randint(0, len(self.allowedList) - 3)
            # get the atom referred to in the allowedList
            axisAtom1Index = self.allowedList[index]
            
            # find the relevant atom positions and rotation axis
            atom1 = xyzVals[axisAtom1Index]
            try:
                atom2 = xyzVals[axisAtom1Index + 1]
            except IndexError:
                print "Index Error"
            rotAxis = atom2 - atom1
            rotAxisHat = rotAxis/np.linalg.norm(rotAxis)
            
            # pick a step size at random from within the current allowed range
            angle = rnd.uniform(maxStepRange * -np.pi, maxStepRange * np.pi)

            # rotate all the remaining points about the atom1 axis place at atom1 by angle
            xyzVals = [ p if n < axisAtom1Index + 2 else cart.rotPAboutAxisAtPoint(p, atom1, rotAxisHat, angle) for n, p in enumerate(xyzVals) ]
        return xyzVals 
        
    def foldInsideEnvelope(self, xyzVals):
        
        # Perform a random dihedral twist using the allowed list. If there are points intersecting 
        # reject it out of hand.
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
        curMin = 0
        while minNumIndicesOutside > 0 and numMoves < self.maxNumFoldingMoves:
            
            # Do the dihedral move on the current working set. 
            newXYZVals = self.dihedralTwist(curXYZVals, maxStepScale)
            
            # if new vals self-intersect then reject it
            for n, pos in enumerate(newXYZVals[0:-1]):
                goodPos = self.checkPointAgainstList(newXYZVals[n+1:], pos)
                if goodPos==False:
                    newXYZVals = curXYZVals
                    break
            
            # find indices outside the envelope and count them.        
            newIndicesOutside = self.checkEnvelope(newXYZVals)
            newNumIndicesOutside = len(newIndicesOutside)
        
            # if a new global minimum number of indices inside the envelope then keep the move.
            if newNumIndicesOutside < minNumIndicesOutside:
                curMin+=1
                # keep a permanent log of the the new best set.
                minXYZVals = newXYZVals[:]
                minNumIndicesOutside = newNumIndicesOutside
                minIndices = newIndicesOutside
                
                # update the working copy
                curXYZVals = newXYZVals[:]
                curNumIndicesOutside = newNumIndicesOutside
                
                if self.dumpInterimFiles==1:
                    fIO.saveXYZ(minXYZVals, 'C', 'min'+str(curMin)+'.xyz')
                
                
                if self.verbose==1:
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

            if numMoves % 10==0 and self.verbose==1:                 
                print "step: ", numMoves, " minNumIndicesOutside: ", minNumIndicesOutside, "curNumIndicesOutside:", curNumIndicesOutside
        
        if minNumIndicesOutside > 0:
            print "Warning: there are points outside the envelope that were not moved inside."
            if self.dumpInterimFiles:
                xyzVals = [ minXYZVals[index] for index in minIndices]
                fIO.saveXYZ(xyzVals, 'B', 'outsideEnvelope.xyz')
                
                
        return minXYZVals

    def checkPointAgainstList(self, checkList, pos):
        ''' Function returns true if the test position is less than minDist distance 
        from points in checkList.'''
        # assume test position is bad
        goodPos = True

        for zPos in checkList:
            # if new position occurs within minDist of any of the existing positions in the list 
            # then throw a wobbly.
            if np.linalg.norm((zPos - pos)) <= self.minDist:
                goodPos = False
                if self.verbose==1: 
                    print("Packing violation")

        return goodPos        
            
    def getParams(self):
        return self.params 

if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the generator
    RandomPolymerPackBBG = RandomPolymerPackBBG(filename)
    
    # generate backbone realtime parameters
    numPoints = 70
    pointA = np.array([0.0, 0.0, 0.0])
    minDist = 1.0
    bondLength = 1.5
    
    alpha1 = 45 # 45.0
    alpha2 = 65 # 75.0
    beta1 = 110.0
    beta2 = 140.0
    
    envelopeList = ['frustum 40 15 0 5']
    
    # generate a curve between the specifed points
    RandomPolymerPackBB = RandomPolymerPackBBG.generateBuildingBlock(numPoints, 
                                                                     pointA,
                                                                     alpha1,
                                                                     alpha2,
                                                                     beta1,
                                                                     beta2,
                                                                     minDist, 
                                                                     bondLength,
                                                                     envelopeList=envelopeList,
                                                                     visualiseEnvelope=(100000, 100))
    RandomPolymerPackBB .exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "RandomPolymer Done"