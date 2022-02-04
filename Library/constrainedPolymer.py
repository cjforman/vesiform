'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
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
    Overloads the building block generator class to allow construction of a polymer chain 
    between two points, whose pathway is randomised and packed into a confined space, while
    avoiding certain points and ensuring that the bond lengths are preserved between each 
    particle in the chain.
    
    An intial arbitrary chain is created starting at one of the points and stretching out
    # which has well defined bond angles and dihedrals.
    
    The particle chain is then manipulated by picking a random bond on the chain and 
    rotating the entire free end of the chain about that bond. Such a move preserves bond lengths 
    between the particles and any move that causes self-intersections or intersections with 
    pointsToAvoid is rejected.  A PE function is defined between the end of the chain and 
    pointB.  
    
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
        self.distEpsilon= self.getParam("distEpsilon")
        self.distanceScale = self.getParam("distanceScale")
        self.maxNumConnectingMoves = self.getParam("maxNumConnectingMoves")
        self.springConstant = self.getParam("springConstant")
        self.connectingTemp = self.getParam("connectingTemp")
        self.beta = self.getParam("beta") * np.pi/180.0
        self.alpha = self.getParam("alpha") * np.pi/180.0
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for constrained polymer object")
            sys.exit()

    def generateBuildingBlock( self, 
                               numPoints,
                               pointA,
                               pointB, 
                               minDist,
                               bondLength,
                               numCrankMoves,
                               pointsToAvoid=[],
                               visualiseEnvelope=(0,20,'envelope.xyz'),
                               envelopeList=["None"],
                               angularRange=["None"],
                               startDirector=["None"]):

        self.numPoints = numPoints
        self.pointA = pointA
        self.pointB = pointB
        self.baseLength = np.linalg.norm(pointA - pointB)
        self.bondLength = float(bondLength)
        self.minDist = minDist
        self.numCrankMoves = numCrankMoves
        self.angularRange = angularRange
        self.startDirector = startDirector
        self.blockNames = self.generateBuildingBlockNames()
        self.allowedList = self.generateAllowedList()
        blockDirector = self.generateBuildingBlockDirector()
        self.blockDirectorHat = blockDirector/np.linalg.norm(blockDirector)
        self.blockRefPoint= self.generateBuildingBlockRefPoint()
        self.parseEnvelopeList(envelopeList)
        
        
        # from the points to avoid list only remember those that would pass the specified envelope test
        self.pointsToAvoid = pointsToAvoid

        # check starting points are legal or it's gonna be a long wait.
        if not self.checkPointInBounds(self.pointA):
            print("Error Warning: PointA out of bounds")
            time.sleep(3)
             
        if not self.checkPointInBounds(self.pointB):
            print("Error Warning: PointB out of bounds")
            time.sleep(3)


        return BBG.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, visualiseEnvelope=visualiseEnvelope, pointsToAvoid=pointsToAvoid)

    def generateBuildingBlockDirector(self):
        if self.startDirector[0] == 'None':
            director = self.pointB - self.pointA
        else:
            director = self.startDirector
        return  director/np.linalg.norm(director)
    
    def generateBuildingBlockRefPoint(self):
        return (self.pointA + self.pointB)/2.0

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints

    def checkBondLengths(self):
        print([ np.linalg.norm(self.buildingBlockXYZ[i] - self.buildingBlockXYZ[i-1]) for i in range(1, len(self.buildingBlockXYZ)) ])       
        
    def generateAllowedList(self, short=False):
        return [i for i in range(0, self.numPoints) ]
        
    def generateBuildingBlockXYZ(self):
        print("Generate initial conformation.")
        xyzVals = self.generateSpaceCurve()
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'initialChain.xyz')

        print("Minimising End Point with dihedral moves on allowed list")
        # perform the energy minimisation that moves the free end to pointB
        xyzVals = self.minimiseEnergy(xyzVals)
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'pointBMinimised.xyz')

        print("Randomising chain using crankshaft moves.")
        xyzVals, numValidMoves = self.crankShaftMoves(xyzVals, self.numCrankMoves, 1)

        if self.dumpInterimFiles==1 and self.numCrankMoves > 0:
            fIO.saveXYZList(xyzVals, self.blockNames, 'crankedChain.xyz')
        
        if self.maxNumFoldingMoves>0:
            print("Folding structure up.")
            xyzVals = self.foldInsideEnvelope(xyzVals)

        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'foldedChain.xyz')

        return xyzVals
    
    def generateSpaceCurve(self):
        # generates a space curve with a specified bond angle, dihedrals and bondlength
        
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
            
            if self.angularRange[0]=='None':
                # compute a new direction based on beta and alpha
                alpha = self.alpha
                beta = self.beta
            else:
                # if the angular range is set then pick a random direction. (hope it's not too madly self intersecting. Got a bit chill about that later on.
                # careful narrow choices of ranges can give good results without checking for intersections.
                alpha = rnd.uniform(self.angularRange[0], self.angularRange[1])
                beta = rnd.uniform(self.angularRange[2], self.angularRange[3])
            dirn = coords.generateTNBVecXYZ(tnb, beta, alpha)
            
            # construct next space curve point 
            spaceCurve.append(spaceCurve[-1] + b * dirn/np.linalg.norm(dirn))
        
        currentAxis = coords.axisFromHelix(spaceCurve[2:])
        currentRefPoint = self.pointA
        
        spaceCurve = coords.transformFromBlockFrameToLabFrame(self.blockDirectorHat, self.pointA, 0.0, currentAxis, currentRefPoint, spaceCurve[2:])
        
        return spaceCurve
    
    def cleanUpZeroes(self, xyzVals):
        return np.around(xyzVals, 5)

    def minimiseEnergy(self, xyzVals):
        # Returns a polymer chain which minimises a simple PE function.
        # Performs random dihedral twists on the free end of a polymer, 
        # starting at a random bond.
        # Moves resulting in lower energy arrangements are accepted.
        # Moves resulting in higher energy arrangements are accepted
        # with probability that is exponentially smaller with increasing energy difference.
        # Since the primary structure is just a straight line and we are only
        # doing sparse dihedral twists, with only a single bias towards to the pointB,
        # there is a low probability of a self-intersection especially for long chains.
        # This rapidly finds a structure where the point B is within arbitrary distance of pointB.
        # The step size scales in proportion to the distance from point B. So only tiny steps are taken 
        # near to the only minimum of the entire potential. Converges fairly rapidly even for large N.
         
        lowestEnergyMinimum = xyzVals[:]
        initPE, initDist = self.PE(xyzVals)
        curPE = initPE
        minPE = initPE
        curDist = initDist
        minDist = initDist
        maxStepRange = 1.0
        numMoves = 0
        curMin = 0
        while minDist> self.distEpsilon and numMoves < self.maxNumConnectingMoves:

            # compute new conformation based on a random dihedral twist
            newXYZ = self.dihedralTwist(xyzVals, maxStepRange)

            # compute energy and distance of new move 
            newPE, newDist = self.PE(newXYZ)   

            # compute energy difference with current minimum PE
            deltaPE = newPE - minPE

            # assume we will accept the move.   
            acceptMove = True 
            # if the currentPE is greater than the minimum then only accept
            # the move with a probability given by the difference in 
            # energy of the minimum and current states  
            if deltaPE > 0 :
                # pick a random value between 0 and 1
                prob = rnd.uniform(0,1)

                # if that value is larger than the threshold reject the move.
                # The threshold decreases with increasing deltaE, so the 
                # higher the energy of the new state relative to the old one
                # the more likely it is we reject the move 
                if prob > np.exp(-deltaPE/self.connectingTemp):
                    acceptMove = False
                
            # if we accept the move then store the new coords
            # and record the energy of the newest accepted move    
            if acceptMove:
                xyzVals = newXYZ[:]
                curPE = newPE
                curDist = newDist
                
            # check the curPE against the minimum energy
            # if we have a new min energy then retain for the future and dump to file.
            if curPE < minPE:
                lowestEnergyMinimum = xyzVals[:]
                minPE = curPE
                minDist = curDist
                
                maxStepRange = min(1.0, minDist/self.distanceScale)
          
                # if the minDist is close then regenerate the allowed list with the short parameter set.
                # In this case only make moves in the ten residues closest to the target.
                #if minDist<3.0:
                #    self.generateAllowedList(short=True)
          
                curMin += 1 
                if self.dumpInterimFiles>1:
                    self.outline(numMoves, self.maxNumConnectingMoves, minDist, minPE, maxStepRange )
                if curMin <= 20 and self.dumpInterimFiles>1:
                    fIO.saveXYZ(lowestEnergyMinimum, 'Be', 'min_' + str(curMin) + '.xyz')
                if curMin > 20 and curMin % 10 ==0 and self.dumpInterimFiles>1:
                    fIO.saveXYZ(lowestEnergyMinimum, 'Be', 'min_' + str(curMin) + '.xyz')
                if curMin > 100 and curMin % 100 ==0 and self.dumpInterimFiles>1:
                    fIO.saveXYZ(lowestEnergyMinimum, 'Be', 'min_' + str(curMin) + '.xyz')
                
            numMoves += 1
            
            if numMoves % 100 == 0: 
                self.outline(numMoves, self.maxNumConnectingMoves, minDist, minPE, maxStepRange )                 
        
        # regenerate the allowed list at full length
        self.generateAllowedList(short=False)
        
        return lowestEnergyMinimum

    def outline(self, n, M, d, E, R):
        print(n, "out of ", M, "minDist:", d, "minEnergy:", E, "maxStepRange:", R)
    
    def PE(self, xyzVals):
        PE = 0.0
        # add spring between pointB and end of xyzVals with equilibrium at pointB
        dist1 = np.linalg.norm(xyzVals[-1] - self.pointB)
        PE += 0.5 * 3 * self.springConstant * dist1**2
        
        return PE, dist1
        
    def dihedralTwist(self, xyzVals, maxStepRange):  
        
        # can only do this if there are sufficient atoms in the array
        if len(xyzVals) > 3:
            # pick a random point in the allowedList
            index = rnd.randint(0, len(self.allowedList) - 2)
            # get the atom referred to in the allowedList
            axisAtom1Index = self.allowedList[index]
            
            # find the relevant atom positions and rotation axis
            atom1 = xyzVals[axisAtom1Index]
            try:
                atom2 = xyzVals[axisAtom1Index + 1]
            except IndexError:
                print("Index Error")
            rotAxis = atom2 - atom1
            rotAxisHat = rotAxis/np.linalg.norm(rotAxis)
            
            # pick a step size at random from within the current allowed range
            angle = rnd.uniform(maxStepRange * -np.pi, maxStepRange * np.pi)

            # rotate all the remaining points about the atom1 axis place at atom1 by angle
            xyzVals = [ p if n < axisAtom1Index + 2 else cart.rotPAboutAxisAtPoint(p, atom1, rotAxisHat, angle) for n, p in enumerate(xyzVals) ]
        return xyzVals 
        
    def foldInsideEnvelope(self, xyzVals):
        
        # Perform a random crank shaft using the allowed list. If there are points intersecting reject it out of hand.
        # if there are less points outside the zone accept it.
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
        threshold = 1.0
        
        while minNumIndicesOutside > 0 and numMoves < self.maxNumFoldingMoves:
            
            # Do the crank shaft move on the current working set. 
            # Already rejects crossovers with pointsToAvoid and selfcrossovers.
            # Will return the current set in those cases so no move on number of points.   
            newXYZVals, numValidMoves = self.crankShaftMoves(curXYZVals, 1, maxStepScale)

            # find indices outside the envelope and count them.        
            newIndicesOutside = self.checkEnvelope(newXYZVals)
            newNumIndicesOutside = len(newIndicesOutside)
        
            if numValidMoves==1:
                acceptedMove = False
                # if a new global minimum number of indices inside the envelope then keep the move.
                if newNumIndicesOutside < minNumIndicesOutside:
                    acceptedMove = True
                    # keep a permanent log of the the new best set.
                    minXYZVals = newXYZVals[:]
                    minNumIndicesOutside = newNumIndicesOutside
                    minIndices = newIndicesOutside
                    
                    # update the working copy
                    curXYZVals = newXYZVals[:]
                    curNumIndicesOutside = newNumIndicesOutside
                    
                    curMin += 1
                    if curMin<10 and self.dumpInterimFiles==1:
                        fIO.saveXYZList(minXYZVals, self.blockNames, "foldMin" + str(curMin) + '.xyz')
                    if curMin>10 and self.dumpInterimFiles==1 and curMin % 10==0:
                        fIO.saveXYZList(minXYZVals, self.blockNames, "foldMin" + str(curMin) + '.xyz')
                    if curMin>100 and self.dumpInterimFiles==1 and curMin % 100==0:
                        fIO.saveXYZList(minXYZVals, self.blockNames, "foldMin" + str(curMin) + '.xyz')
                
                 
                # if the new move means that the num indices outside has gone up above the curXYZVals then 
                # roll the die to decide whether or not to replace the curXYZVals
                if newNumIndicesOutside >= minNumIndicesOutside:
                    # roll the die. The larger the gap between newNum and minNum the smaller the threshold and the less likely we will accept the move and go back to the curXyzVals. 
                    threshold = np.exp(-(float(newNumIndicesOutside - minNumIndicesOutside)) / self.foldingTemp)
                    if rnd.uniform(0.0, 1.0) < threshold: 
                        acceptedMove = True
                        # update the working copy of coords but don't replace the current global minimum
                        curXYZVals = newXYZVals
                        curNumIndicesOutside = newNumIndicesOutside 
            
                acceptedString = "rejected"
                if acceptedMove==True:
                    acceptedString = "accepted"
                print("step:", numMoves, "out of", self.maxNumFoldingMoves, "minNumIndicesOutside:", minNumIndicesOutside, "curNumIndicesOutside:", curNumIndicesOutside, "threshold", threshold, acceptedString)
            numMoves += 1

        if minNumIndicesOutside > 0:
            print("Warning: there are points outside the envelope that were not moved inside.")
            if self.dumpInterimFiles:
                fIO.saveXYZ([ xyzVals[index] for index in minIndices], 'B', 'outsideEnvelope.xyz')
                
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
        numValidMoves = 0
        while numMoves < numCrankShaftMoves:

            # select some indices at random
            index1 = 0
            index2 = 0            
            
            # check for equal or adjacent indices and indices on the allowed list
            while (abs(index1-index2) < 2):

                # pick two points in the allowed list at random
                allowedIndex1 = rnd.randint(0, len(self.allowedList) - 1)
                allowedIndex2 = rnd.randint(0, len(self.allowedList) - 1)
                
                # get the index values
                index1 = self.allowedList[allowedIndex1]
                index2 = self.allowedList[allowedIndex2]

            indexMin = min([index1, index2])
            indexMax = max([index1, index2])

            # performs a single crank on the given indices. 
            # If move is rejected due to overlaps then the original array is returned
            workingXYZVals, validCrank = self.crankShaftMove(workingXYZVals, indexMin, indexMax, maxStepScale)

            # wait until we've have numCrankShaft valid moves
            if validCrank:
                numValidMoves += 1

            numMoves += 1
            if numMoves % 20 == 0:
                print(numMoves, " out of ", numCrankShaftMoves) 
       
        return workingXYZVals, numValidMoves

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
            try:
                workingXYZVals = np.concatenate((lowList, newVals, highList), 0)
            except ValueError:
                validCrank = False
               
        # if valid Crank is false then the original list is returned intact
        # if valid crank is true then the new vals are inserted into the new list
        return workingXYZVals, validCrank

    def checkPointInBounds(self, pos, ignorePTA=False, pointsToAvoid=[]):
        inBounds = True

        # in this envelope any point outside the ellipse is out of bounds
        # The ellipse has three axes Rx, Ry and Rz which are specified
        # in the envelope list.  The ellipse is oriented such that its
        # x-axis is aligned with the blockDirector and rotated about that 
        # axis by an amount specified in the envelope list
        # 
        try:
            envelopeParams = self.envelopeSummary['outerellipse']
            # ellipse centered at refPoint with input radius rx, ry and rz
            # and oriented so it's x-axis is along the pointA to pointB line.
            
            # first transform the point into the ellipse coordinates
            blockPos = coords.transformFromLabFrameToBlockFrame(self.blockDirectorHat, self.blockRefPoint, envelopeParams[3], np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]), [pos]) 
            
            # compute spherical polar coords of pos relative to the ellipsoid body axes and centre.
            polarPos = coords.XYZ2SphericalPolar(blockPos[0])
            # compute distance from centre of ellipsoid to ellipsoids surface at the bearing as the test pos.
            r_e = coords.ellipsoidalPolarUVW(polarPos[1], polarPos[2], envelopeParams[0], envelopeParams[1], envelopeParams[2])[0] # take first parameter of returned information which is distance to ellipsoid
            if polarPos[0] > r_e:
                if self.verbose==1:
                    print("outerEllipse Violation")
                inBounds=False
        except KeyError:
            pass

        # check the underlying checkPointInBounds for other checks that may be included
        if inBounds:
            inBounds = BBG.checkPointInBounds(self, pos, ignorePTA=ignorePTA, pointsToAvoid=pointsToAvoid)
        
        return inBounds
        
        
    def getParams(self):
        return self.params 

if __name__ == '__main__':
    from Library.SurfacePackEllipsoid import SurfacePackEllipsoidBBG as SPEBBG
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the generator
    ConstrainedPolymerPackGBB = ConstrainedPolymerPackBBG(filename)
    
    # generate backbone realtime parameters
    numPoints = 50
    pointA = 11.0 * np.array([1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 1.0/np.sqrt(3)])
    pointB = 11.0 * np.array([-1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), -1.0/np.sqrt(3)])
    minDist = 1.0
    bondLength = 1.5
    crankMoves = 20
    spherePointGenerator = SPEBBG('../../Library/SurfacePackEllipsoid/SurfacePackEllipsoid.txt')
    spherePoints = spherePointGenerator.generateBuildingBlock(30, 11, 11, 11, -90, 90, -180, 180, 4)
    pointsToAvoid = spherePoints.getAtomsXYZ()
    
    fIO.saveXYZList([pointA, pointB], ['Ca', 'O'], 'labPoints.xyz')
    fIO.saveXYZ(pointsToAvoid, 'Na', 'labPointsToAvoid.xyz')
    
    # generate a curve between the specifed points
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
    print("constrainedPolymer Done")