import sys
import numpy as np
import time
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
import Utilities.cartesian as cart
import itertools as it
import random as rnd
from Library.constrainedPolymer import ConstrainedPolymerPackBBG as CPBBG
from Library.SpidroinBackbone import spidroinBackboneGenerator as SBBG
from numpy import inf

class spidroinHairpinGenerator(CPBBG):
    ''' '''
      
    def __init__(self, paramFilename):
        # initialise the parameter dictionary for the base classes
        CPBBG.__init__(self, paramFilename)
        
    def initialiseParameters(self):
        # initialise the constrained polymer parent
        CPBBG.initialiseParameters(self)
        
        self.SP1NumGUnits =self.getParam('SP1NumGUnits')
        self.SP2NumGUnits =self.getParam('SP2NumGUnits')        
        self.energyScale = self.getParam('energyScale')
        self.maxNumE2Moves = self.getParam('maxNumE2Moves')
        self.E2Temp = self.getParam('E2Temp')
        self.GEpsilon = self.getParam('GEpsilon') 
        self.GRm = self.getParam('GRm')
        self.PQEpsilon = self.getParam('PQEpsilon') 
        self.PQRm = self.getParam('PQRm')        
        self.PQGEpsilon = self.getParam('PQGEpsilon') 
        self.PQGRm = self.getParam('PQGRm')        

        # load the backbone building object
        self.SBBG = SBBG(self.paramFilename)
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for hairpin"
            sys.exit()        

    def generateBuildingBlock(  self, 
                                species, 
                                pointA,
                                pointB,
                                minDist,
                                numCrankMoves, 
                                envelopeList = ['None'],
                                visualiseEnvelope = (0,200)):        
        self.species = species
        if species=='SP1':
            self.numGUnits = self.SP1NumGUnits 
            self.numPQUnits = self.SP1NumGUnits
        else:
            self.numGUnits = self.SP2NumGUnits 
            self.numPQUnits = self.SP2NumGUnits + 1
    
        self.numPoints = self.numGUnits * 3 + self.numPQUnits * 3  
        self.minDist = minDist
        
        # set up the right reference points
        self.pointA = pointA
        self.pointB = pointB
        
        # calculate the block director
        self.blockDirectorHat = self.generateBuildingBlockDirector()
        
        # set the variables that need to be set to perform checking of pointsA and B
        self.blockRefPoint = self.generateBuildingBlockRefPoint()
        self.parseEnvelopeList(envelopeList)
        self.pointsToAvoid = []
         
        # check starting points are legal or it's gonna be a long wait.
        if not self.checkPointInBounds(self.pointA):
            print "Error Warning: pointA out of bounds"
            time.sleep(3)
         
        if not self.checkPointInBounds(self.pointB):
            print "Error Warning: pointB out of bounds"
            time.sleep(3)

        bondLength = 1.0
        return CPBBG.generateBuildingBlock(self, self.numPoints, pointA, pointB, minDist, bondLength, numCrankMoves, envelopeList=envelopeList, visualiseEnvelope=visualiseEnvelope) 

    def generateAllowedList(self, short=False):
        # ensure the names are named correctly 
        names = self.generateBuildingBlockNames()
       
        # only allowed to do crank shafts or dihedral twists at points named C or N.
        return [ i for i, name in enumerate(names) if name in ['C', 'N'] ]

    def generateSpaceCurve(self):
        # Over-rides the generate space curve function of the parent to generate a spidroin coarse grained backbone
        spidroinBackbone = self.SBBG.generateBuildingBlock(self.minDist, self.species)
        # place the first atom in the building block at point A
        spidroinBackbone.placeAtom(0, self.pointA)

        # find the current axis of the system        
        currentAxis = coords.axisFromHelix(spidroinBackbone.getAtomsXYZ())
        currentRefPoint = self.pointA
        
        # return the coords such that the system is rotated along the pointA pointB axis with the first point at pointA
        return coords.transformFromBlockFrameToLabFrame(self.blockDirectorHat, self.pointA, 0.0, currentAxis, currentRefPoint, spidroinBackbone.getAtomsXYZ())
          

    def foldInsideEnvelope(self, xyzVals):
        # fold inside the envelope conventionally with the overloaded crankshaft move
        xyzVals = CPBBG.foldInsideEnvelope(self, xyzVals)

        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'FinalFolded.xyz')

        xyzVals = self.energyMinimise(xyzVals)

        if self.dumpInterimFiles==1:
            fIO.saveXYZList(xyzVals, self.blockNames, 'EnergyMinimised.xyz')

        return xyzVals

    def energyMinimise(self, xyzVals):

        # establish the force points for the chain
        self.forcePointsG = [ i for i, name in enumerate(self.blockNames) if name=='B' ]
        self.forcePointsPQ = [ i for i, name in enumerate(self.blockNames) if name=='P' ]
        
        # Use crank shaft moves to find lower energy configurations of the force points in the spidroin backbone.
        lowestEnergyConfiguration = xyzVals[:]
        curXYZVals = xyzVals[:]
        initPE = self.PE2(xyzVals)
        curPE = initPE
        minPE = initPE
        initOOB = len(self.checkAllPointsInBounds(xyzVals))
        deltaMinEnergy = inf
        maxStepRange = 1.0
        numMoves = 0
        curMin = 0
        while numMoves < self.maxNumE2Moves:

            # compute new conformation based on a random crankshaft move. If it self intersects then reject it
            newXYZ, numValidMoves = self.crankShaftMoves(curXYZVals, 1, maxStepRange)

            # reject any move that moves things out of the envelope.
            numOOB = len(self.checkAllPointsInBounds(xyzVals))
            if  numOOB > initOOB:
                numValidMoves = 0

            # only compute a new energy if we moved successfully 
            if numValidMoves==1:
                # assume the move won't be accepted
                acceptMove = False
                
                # compute new energy and out of bounds 
                newPE = self.PE2(newXYZ)   
              
                # compute energy difference with current lowest PE
                deltaPE = newPE - minPE

                # if energy decreased then keep as new minimum 
                if (deltaPE < 0.0):
                    # we have accepted the move.
                    acceptMove = True
                    
                    # make a note of what is now the best configuration; this will have lowest energy yet encountered.
                    lowestEnergyConfiguration = newXYZ[:]

                    # set the input for the next loop to the new structure
                    curXYZVals = newXYZ[:]

                    # make a note of the difference between this minimum energy and the last minimum energy
                    deltaMinEnergy  = deltaPE
                    
                    # record the newPE as minPE and curPE
                    minPE = newPE
                    curPE = newPE
                    
                    # set the new step range based on how much the energy has changed relative to the specified energy scale. This is largely guess work.
                    maxStepRange = min(0.1, 2.0 - 2.0 / (1.0 + np.exp( -(initPE - minPE)/self.energyScale)))
                    
                    # increment the minimum number and dump some stuff to file every now and then.
                    curMin += 1 
                    if curMin <= 20 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'foldMin_' + str(curMin) + '.xyz')
                    if curMin > 20 and curMin % 10 ==0 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'foldMin_' + str(curMin) + '.xyz')
                    if curMin > 100 and curMin % 100 ==0 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'foldMin_' + str(curMin) + '.xyz')

                # If the energy has increased then keep a working copy with a certain probability      
                if (deltaPE > 0):
                    # compute energy threshold. The larger the difference in energy (always +ve in this case) between newPE and the current Minimum the lower the threshold 
                    energyThreshold = np.exp(-deltaPE/self.E2Temp)
                    
                    # pick a random value between 0 and 1
                    diceRoll = rnd.uniform(0,1)
                
                    # if the dice Roll is less than the threshold then accept the newPE as the input for the next loop.
                    if diceRoll < energyThreshold:
                        acceptMove = True
                        curXYZVals = newXYZ[:] # set the new vals as the current copy
                        curPE = newPE # set curPE and curOOB
        
                acceptString="rejected"
                if acceptMove==True:
                    acceptString="accepted"
                print numMoves,"out of", self.maxNumE2Moves, "deltaPE:", deltaPE, "curEnergy:", curPE, "minEnergy:", minPE, "deltaMinE:", deltaMinEnergy, "maxStepRange:", maxStepRange, "numOOB", numOOB, "initOOB", initOOB, acceptString 
            
            numMoves += 1

        return lowestEnergyConfiguration  


    def DualEnergyFold(self, xyzVals):

        # fold the structure inside any given envelope
        # print "Folding structure into envelope"
        # xyzVals= CPBBG.foldInsideEnvelope(self, xyzVals)
        
        # then minimise the structure for energy by performing crankmoves

        print "minimising energy of structure"

        self.forcePointsG = [ i for i, name in enumerate(self.blockNames) if name=='B' ]
        self.forcePointsPQ = [ i for i, name in enumerate(self.blockNames) if name=='P' ]
        
        # overloads the original fold function to become an energy minimiser that minimises a second potential energy function.
        # Use crank shaft moves to find lower energy configurations of the force points in the spidroin backbone.
        lowestEnergyConfiguration = xyzVals[:]
        curXYZVals = xyzVals[:]
        initPE = self.PE2(xyzVals)
        curPE = initPE
        minPE = initPE
        initOOB = len(self.checkAllPointsInBounds(xyzVals))
        curOOB = initOOB
        minOOB = initOOB
        deltaMinEnergy = 1000 * self.energyEpsilon
        maxStepRange = 1.0
        numMoves = 0
        curMin = 0
        while (np.abs(deltaMinEnergy) > self.energyEpsilon or minOOB > 0) and numMoves < self.maxNumE2Moves:

            # compute new conformation based on a random crankshaft move. If it self intersects then reject it
            newXYZ, numValidMoves = self.crankShaftMoves(curXYZVals, 1, maxStepRange)

            # only compute a new energy if we moved successfully 
            if numValidMoves==1:
                # assume the new move won't be accepted
                acceptMove = False
                
                # compute new energy and out of bounds 
                newPE = self.PE2(newXYZ)   
                newOOB = len(self.checkAllPointsInBounds(newXYZ))                

                # compute energy difference with current lowest PE
                deltaPE = newPE - minPE

                # compute OOB difference
                deltaOOB = newOOB - minOOB

                # there are nine possibilities,
                # 1) OOB decreased and energy went down. Action: Keep as new minimum 
                # 2) OOB decreased and energy stayed the same. Action: Keep as new minimum
                # 3) OOB decrease and energy went up  
                # 4) OOB Stayed the same and energy went down. Action: Keep as new minimum
                # 5) OOB Stayed the same and energy stayed the same. Keep as working copy not a new minimum                    
                # 6) OOB stayed the same and energy went up. Keep as working copy with energy threshold probability
                # 7) OOB went up and energy went down
                # 8) OOB went up and energy stayed the same. Keep as working copy with OOB threshold probability
                # 9) OOB went up and energy went up. Keep as working copy with lowest probability of OOB or energy.

                # Keep as new minimum:  Options 1, 2 and 4. 
                if ((deltaOOB < 0) and (deltaPE < 0.0)) or (( deltaOOB < 0 ) and ( np.abs(deltaPE)<1e-12)) or ((deltaOOB == 0) and (deltaPE < 0.0)):
                    # we have accepted the move
                    acceptMove = True
                    
                    # make a note of what is now the best configuration; this will have better energy and the same OOB or better OOB and the same energy.
                    lowestEnergyConfiguration = newXYZ[:]

                    # also set the working configuration to the new structure
                    curXYZVals = newXYZ[:]

                    # if it was the energy that went down then retain the difference between this and the last best minimum for energy
                    # when changes in energy start getting super tiny we're maybe close to a minimum
                    if deltaPE < 0.0: 
                        deltaMinEnergy  = deltaPE
                    
                    # record the new minPE, minOOB and curPE and curOOB.
                    minPE = newPE
                    curPE = newPE
                    curOOB = newOOB
                    minOOB = newOOB
                    
                    # set the new step range based on how close to fully folded we are in the minimum. only update this when we update the minimum
                    maxStepRange = min(1.0, float(minOOB)/float(initOOB))
                    
                    # increment the minimum number and dump some stuff to file every now and then.
                    curMin += 1 
                    # self.outline2(numMoves, self.maxNumE2Moves, curOOB, minOOB, deltaPE, minPE, maxStepRange )
                    if curMin <= 20 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'foldMin_' + str(curMin) + '.xyz')
                    if curMin > 20 and curMin % 10 ==0 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'foldMin_' + str(curMin) + '.xyz')
                    if curMin > 100 and curMin % 100 ==0 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'foldMin_' + str(curMin) + '.xyz')

                # keep as working copy - option 5 - No Change!
                if (deltaOOB==0) and (np.abs(deltaPE)<1e-12):
                    # we have accepted the move
                    acceptMove = True
                    curXYZVals = newXYZ[:] # set the new vals as the current copy
                    curPE = newPE # set curPE and curOOB
                    curOOB = newOOB
                
                # need these for the next options
                energyThreshold = np.exp(-deltaPE/self.E2Temp)
                OOBThreshold = np.exp(-float(deltaOOB)/self.foldingTemp)
                # pick a random value between 0 and 1
                diceRoll = rnd.uniform(0,1)
                   
                # keep as working copy with probability of lowest threshold - option 9 - both went up      
                if (deltaOOB > 0) and (deltaPE> 0.0):
                    if diceRoll < min(energyThreshold, OOBThreshold):
                        # we have accepted the move
                        acceptMove = True
                        curXYZVals = newXYZ[:] # set the new vals as the current copy
                        curPE = newPE # set curPE and curOOB
                        curOOB = newOOB

                # oob stayed the same and energy went up: option 6                            
                if (deltaOOB==0) and (deltaPE > 0):
                    if diceRoll < energyThreshold:
                        # we have accepted the move
                        acceptMove = True
                        curXYZVals = newXYZ[:] # set the new vals as the current copy
                        curPE = newPE # set curPE and curOOB
                        curOOB = newOOB

                # oob went up and energy stayed the same: option 8
                if (deltaOOB > 0) and ( np.abs(deltaPE) < 1e-12):
                    if diceRoll < OOBThreshold:
                        # we have accepted the move
                        acceptMove = True
                        curXYZVals = newXYZ[:] # set the new vals as the current copy
                        curPE = newPE # set curPE and curOOB
                        curOOB = newOOB

                # OOB went up but energy went down. Dilemma. What to do! Option 3.
                if (deltaOOB > 0 ) and (deltaPE < 0):
                    # accept the new minimum with a probability based on oob.
                    if diceRoll < OOBThreshold:
                        # we have accepted the move
                        acceptMove = True
                        curXYZVals = newXYZ[:] # set the new vals as the current copy
                        curPE = newPE # set curPE and curOOB
                        curOOB = newOOB

                # Energy went up but oob went down. Dilemma. What to do! option 7
                if (deltaOOB < 0 ) and (deltaPE > 0):
                    # accept the new minimum with a probability based on energy but don't make it a bi-minimum
                    if diceRoll < energyThreshold:
                        # we have accepted the move
                        acceptMove = True
                        curXYZVals = newXYZ[:] # set the new vals as the current copy
                        curPE = newPE # set curPE and curOOB
                        curOOB = newOOB

                # output the results of the current crank
                # if numMoves % 100 == 0: 
                self.outline2(numMoves, self.maxNumE2Moves, deltaOOB, curOOB, minOOB, deltaPE, curPE, minPE, maxStepRange, acceptMove, deltaMinEnergy)                 
            # else:
            #     print numMoves, " out of ", self.maxNumE2Moves, "Invalid Crank"
        
            numMoves += 1
            
        return lowestEnergyConfiguration

    def outline2(self, n, M, dOOB, cOOB, mOOB, dE, cE, mE, R, acc, dminE):
        acceptString="rejected"
        if acc==True:
            acceptString="accepted"
        print n,"out of", M, "deltaOOB:", dOOB, "curOOB: ", cOOB, "minOOB:", mOOB, "deltaPE:", dE, "curEnergy:", cE, "minEnergy:", mE, "maxStepRange:", R, acceptString, "deltaMinE:", dminE
    


    def crankShaftMoves(self, xyzValsOrig, numCrankShaftMoves, maxStepScale):
        # Picks two C points or two N points at random and rotates all the objects 
        # between them by a random angle about the axis between the selected points.
        #
        # The angle scale is defined by maxStepScale.
        #
        # Checks to see if any of the objects on the list in the new positions 
        # intersect with any of the other objects. If they do then the move is rejected.
        # This kind of crank shaft move preserves the bond lengths that we went to so much 
        # trouble to get right.

        workingXYZVals = xyzValsOrig[:]

        numMoves = 0 
        numValidMoves = 0
        while numMoves < numCrankShaftMoves:

            # select some indices at random
            index1 = 0
            index2 = 0            
            
            # check for equal or adjacent indices and indices on the allowed list
            while (abs(index1-index2) < 2):

                # pick a point in the allowed list at random
                allowedIndex1 = rnd.randint(0, len(self.allowedList) - 1)
                
                # get the main xyzVals/name list index
                index1 = self.allowedList[allowedIndex1]
                
                # get the name at the first point
                name = self.blockNames[index1]
                
                
                # choose an initial second point
                allowedIndex2 = rnd.randint(0, len(self.allowedList) - 1)
                index2 = self.allowedList[allowedIndex2]
                # keep looping until the names of the first and second points match
                while self.blockNames[index2]!=name:
                    allowedIndex2 = rnd.randint(0, len(self.allowedList) - 1)
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
                print numMoves, " out of ", numCrankShaftMoves 
       
        return workingXYZVals, numValidMoves

    def crankShaftMove(self, xyzVals, indexMin, indexMax, maxStepScale):
        # Performs rotation of a random size of all objects contained in the index range 
        # about the axis between points defined by the two indices.
        # This is known as a crank rotation which preserves the bond lengths.
        
        # select the angle to move.
        angle = rnd.uniform(-maxStepScale * np.pi, maxStepScale * np.pi)
            
        # compute the rotation axis of the particles
        n = xyzVals[indexMax] - xyzVals[indexMin]
        nHat = n/np.linalg.norm(n)
        
        # do the crank rotation of all the intervening points selected by the indices.
        newVals = [ cart.rotPAboutAxisAtPoint(p, xyzVals[indexMin], nHat, angle) for p in xyzVals[indexMin:indexMax + 1]]
        
        # extract the names of the rotated values (helps with debugging!)
        # newNames = [ name for name in self.blockNames[indexMin:indexMax+1]]
        
        # Figure out the index ranges of the low list, working list and high list. 
        # The low and high lists are the static regions before and after the moved (working) part.
        if self.blockNames[indexMin]=='N':
            lowListIndexMin = 0
            lowListIndexMax = indexMin - 1 # Chop the N of the first moved object from the low list
            workListIndexMin = 0 # note the work list indices refer to objects in the newVals list.
            workListIndexMax = len(newVals) - 2 # the end point of the new list is an N that doesn't move and belongs to the first composite object of the highlist.
            highListIndexMin = indexMax # highlist goes from the first N of the high list to the end.
            highListIndexMax = len(xyzVals) - 1 # the end point of the list.
        
        if self.blockNames[indexMin]=='C':
            lowListIndexMin = 0
            lowListIndexMax = indexMin # Low static list is from first N to the last C immediately before movedList.
            workListIndexMin = 1 # note the work list indices refer to objects in the newVals list.
            workListIndexMax = len(newVals) - 1 # the first point of the new list is a C that doesn't move and belongs to the last composite object of the low list
            highListIndexMin = indexMax + 1 # highlist goes from the first N of the high list to the end.
            highListIndexMax = len(xyzVals) - 1 # the end point of the list.

        # generate the workings list of positions and names (the names are for debugging purposes to make sure it's all good)            
        lowXYZList = xyzVals[lowListIndexMin:lowListIndexMax + 1] 
        # lowNameList = self.blockNames[lowListIndexMin:lowListIndexMax + 1]
        workXYZList = newVals[workListIndexMin:workListIndexMax + 1] 
        # workNameList = newNames[workListIndexMin:workListIndexMax + 1]
        highXYZList = xyzVals[highListIndexMin:highListIndexMax + 1]  
        # highNameList = self.blockNames[highListIndexMin:highListIndexMax + 1]
        
        # count the number of triplets in each list
        numLow = int(len(lowXYZList)/3)
        numWork = int(len(workXYZList)/3)
        numHigh = int(len(highXYZList)/3)

        # reshape all the lists into groups of three
        lowXYZListTriplets = [ [lowXYZList[3 * i], lowXYZList[3 * i + 1], lowXYZList[3 * i + 2]] for i in range(0, numLow) ]
        # lowXYZNameTriplets = [ [lowNameList[3 * i], lowNameList[3 * i + 1], lowNameList[3 * i + 2]] for i in range(0, numLow) ]
        workXYZListTriplets = [ [workXYZList[3 * i], workXYZList[3 * i + 1], workXYZList[3 * i + 2]] for i in range(0, numWork) ]
        # workXYZNameTriplets = [ [workNameList[3 * i], workNameList[3 * i + 1], workNameList[3 * i + 2]] for i in range(0, numWork) ]
        highXYZListTriplets = [ [highXYZList[3 * i], highXYZList[3 * i + 1], highXYZList[3 * i + 2]] for i in range(0, numHigh) ]
        # highXYZNameTriplets = [ [highNameList[3 * i], highNameList[3 * i + 1], highNameList[3 * i + 2]] for i in range(0, numHigh) ]
            
        validCrank = True # assume success
        # generate iterator over the object pairs with the new vals in the lower half of the list 
        lowListPairs = it.product(workXYZListTriplets, lowXYZListTriplets)
        # iterate and compare each pair to see if any match
        for pair in lowListPairs:
            if self.testTripletIntersect(pair[0], pair[1]):
                validCrank = False
                break 
            
        if validCrank:
            # generate iterator over the product pairs with the new vals in the upper half of the list 
            highListPairs = it.product(workXYZListTriplets, highXYZListTriplets)
            # iterate and compare each pair to see if any match
            for pair in highListPairs:
                if self.testTripletIntersect(pair[0], pair[1]):
                    validCrank = False
                    break
                    
        # reconstruct workingXYZVals
        if validCrank:
            try:
                retXYZVals = np.concatenate((lowXYZList, workXYZList, highXYZList), 0)
            except ValueError:
                if len(highXYZList)==0:
                    retXYZVals = np.concatenate((lowXYZList, workXYZList), 0)
                elif len(lowXYZList)==0:
                    retXYZVals = np.concatenate((workXYZList, highXYZList), 0)
                else:
                    validCrank = False

        if validCrank == False:
            retXYZVals = xyzVals[:]               
        # if valid Crank is false then the original list is returned intact
        # if valid crank is true then the new vals are inserted into the new list
        return retXYZVals, validCrank

    def testTripletIntersect(self, triplet1, triplet2):
        # tests to see if any part of two spherocylinders intersect.
        # Each spherocylinder has a radius of minDist and consists
        # of a sphere at points 0 and 2 of the triplet
        # and a cylinder connecting those two spheres. 
        
        # This is also known as a capsule-capsule sphere swept volume test 
        # for two line segments. The radius of both is the same
        # and is self.minDist. 
        # only need to find the closest point of approach of the primitives - in this case the line segments along which the spheres are swept.
        intersect = False
        dist2 = cart.closestApproachTwoLineSegmentsSquared(triplet1[0], triplet1[2], triplet2[0], triplet2[2])
        if dist2 < self.minDist * self.minDist:
            intersect= True
        return intersect

    
    def PE2(self, xyzVals):
        PE = 0.0

        GUnitPairs = it.combinations(self.forcePointsG, 2)
        PQUnitPairs = it.combinations(self.forcePointsPQ, 2)
        #GPQPairs = it.product(self.forcePointsPQ, self.forcePointsG)
        
        # add the pairwise G units
        for pair in GUnitPairs:
            PE += self.LJAttr(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.GEpsilon, self.GRm)
            PE += self.LJRep(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.GEpsilon, self.GRm)
        
        # add the pairwise PQ units
        # for pair in PQUnitPairs:
        #    PE += self.LJAttr(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.PQEpsilon, self.PQRm)
        #    PE += self.LJRep(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.PQEpsilon, self.PQRm)
            
        # add the cross over terms - only include a repulsive term
        # for pair in GPQPairs:
        #    curPE = self.LJRep(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.PQGEpsilon, self.PQGRm)
        #    PE += curPE

        return PE

    def LJ(self, r, e, rm):
        rp = rm/r
        rpsq = rp * rp
        rp4 = rpsq * rpsq
        rp6 = rp4 * rpsq
        rp12 = rp6 * rp6  
        return e * rp12 - 2 * e * rp6
    
    def LJRep(self, r, e, rm):
        rp = rm/r
        rpsq = rp * rp
        rp4 = rpsq * rpsq
        rp6 = rp4 * rpsq
        rp12 = rp6 * rp6  
        return e * rp12 
    
    def LJAttr(self, r, e, rm):
        rp = rm/r
        rpsq = rp * rp
        rp4 = rpsq * rpsq
        rp6 = rp4 * rpsq
        return - 2.0 * e * rp6 
    
    def generateBuildingBlockNames(self, ca=False):
        if self.species=='SP1':
            names = ['N', 'B', 'C', 'N', 'P', 'C'] * self.numGUnits
        else:
            names = ['N', 'P', 'C'] + ['N', 'B', 'C', 'N', 'P', 'C'] * self.numGUnits
        return names
    
    def generateBuildingBlockConnectors(self):
        # N connector first, then C Connector 
        return [ [2, 1, 0], [self.numPoints-3, self.numPoints-2, self.numPoints-1] ]
        
if __name__ == "__main__":
    
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the generator objects.
    hairPinGen = spidroinHairpinGenerator(filename)
    
    # generate starting points and move the seed to those points extracting
    # the xyzVals each time. 
    pointA = np.array([10.0, 0.0, 0.0])
    pointB = np.array([-10.0, 0.0, 0.0])
    
    minDist = 1.0
    numCrankMoves = 0
    envelopeList = ['innersphere 9.0', 'frustum 2.0 40.0 -100.0 40.0']
    # build building block and dump to file
    hairpinBuildingBlock = hairPinGen.generateBuildingBlock('SP1', pointA, pointB, minDist, numCrankMoves, envelopeList=envelopeList, visualiseEnvelope=(1000000,200))
    hairpinBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "hairpin done"