import sys
import numpy as np
import time
import Utilities.fileIO as fIO
import itertools as it
import random as rnd
from Library.constrainedPolymer import ConstrainedPolymerPackBBG as CPBBG
from Library.SpidroinBackbone import spidroinBackboneGenerator as SBBG

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
        self.energyEpsilon = self.getParam('energyEpsilon')
        self.foldingTemp = self.getParam('foldingTemp')
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
                                envelopeList = ['None']):        
        self.species = species
        if species=='SP1':
            self.numGUnits = self.SP1NumGUnits 
            self.numPQUnits = self.SP1NumGUnits
        else:
            self.numGUnits = self.SP2NumGUnits 
            self.numPQUnits = self.SP2NumGUnits + 1
    
        self.numPoints = self.numGUnits * 3 + self.numPQUnits * 3  
        self.directorHat = np.array([0.0, 0.0, 1.0])
        self.minDist = minDist
        
        # set up the right reference points
        self.pointA = pointA
        self.pointB = pointB
        
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
        return CPBBG.generateBuildingBlock(self, self.numPoints, pointA, pointB, minDist, bondLength, numCrankMoves) 

    def generateAllowedList(self, short=False):
        # ensure the names are named correctly 
        names = self.generateBuildingBlockNames()
       
        # only allowed to do crank shafts or dihedral twists at points named C or N.
        return [ i for i, name in enumerate(names) if name in ['C', 'N'] ]
        

    def generateSpaceCurve(self):
        # Over-rides the generate space curve function of the parent to generate a peptide backbone
        # and a pseudo energy landscape approach to find an initial chain with the end point fixed
        # at point B.
        # create a regular backBone using PointsA as the first residue
        spidroinBackbone = self.SBBG.generateBuildingBlock(self.minDist, self.species)
        return spidroinBackbone.getAtomsXYZ()  

    def foldInsideEnvelope(self, xyzVals):

        self.forcePointsG = [ i for i, name in enumerate(self.blockNames) if name=='B' ]
        self.forcePointsPQ = [ i for i, name in enumerate(self.blockNames) if name=='P' ]
        
        # overloads the original fold function to become an energy minimiser that minimises a second potential energy function.
        # Use crank shaft moves to find lower energy configurations of the force points in the spidroin backbone.
        lowestEnergyConfiguration = xyzVals[:]
        curXYZVals = xyzVals[:]
        initPE = self.PE2(xyzVals)
        curPE = initPE
        minPE = initPE
        deltaPE = 2 * self.energyEpsilon
        maxStepRange = 1.0
        numMoves = 0
        curMin = 0
        while np.abs(deltaPE) > self.energyEpsilon and numMoves < self.maxNumFoldingMoves:

            # compute new conformation based on a random crankshaft move
            newXYZ, numValidMoves = self.crankShaftMoves(curXYZVals, 1, maxStepRange)

            # only compute a new energy if we actually moved (other DeltaPE will be zero!)
            if numValidMoves==1:
                # compute energy 
                newPE = self.PE2(newXYZ)   
    
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
                    if prob > np.exp(-deltaPE/self.foldingTemp):
                        acceptMove = False
                    
                # if we accept the move then store the new coords
                # and record the energy of the newest accepted move    
                if acceptMove:
                    xyzVals = newXYZ[:]
                    curPE = newPE
                    
                # check the curPE against the minimum energy
                # if we have a new min energy then retain for the future and dump to file.
                if curPE < minPE:
                    lowestEnergyConfiguration = xyzVals[:]
                    minPE = curPE
                    maxStepRange = min(1.0, deltaPE/self.energyScale)
                    
                    curMin += 1 
                    self.outline2(numMoves, self.maxNumFoldingMoves, deltaPE, minPE, maxStepRange )
                    if curMin <= 20 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'min_' + str(curMin) + '.xyz')
                    if curMin > 20 and curMin % 10 ==0 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'min_' + str(curMin) + '.xyz')
                    if curMin > 100 and curMin % 100 ==0 and self.dumpInterimFiles==1:
                        fIO.saveXYZ(lowestEnergyConfiguration, 'Be', 'min_' + str(curMin) + '.xyz')
                    
            numMoves += 1
            
            if numMoves % 100 == 0: 
                self.outline2(numMoves, self.maxNumFoldingMoves, deltaPE, minPE, maxStepRange )                 
        
        return lowestEnergyConfiguration

    def outline2(self, n, dE, M, E, R):
        print n, "out of ", M, "deltaE", dE, "minEnergy:", E, "maxStepRange:", R
    
    
    def PE2(self, xyzVals):
        PE = 0.0

        GUnitPairs = it.combinations(self.forcePointsG, 2)
        PQUnitPairs = it.combinations(self.forcePointsPQ, 2)
        GPQPairs = it.product(self.forcePointsPQ, self.forcePointsG)
        
        # add the pairwise G units
        for pair in GUnitPairs:
            PE += self.LJAttr(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.GEpsilon, self.GRm)
            PE += self.LJRep(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.GEpsilon, self.GRm)
        
        # add the pairwise PQ units
        for pair in PQUnitPairs:
            PE += self.LJAttr(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.PQEpsilon, self.PQRm)
            PE += self.LJRep(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.PQEpsilon, self.PQRm)
            
        # add the cross over terms - only include a repulsive term
        for pair in GPQPairs:
            PE += self.LJRep(np.linalg.norm(xyzVals[pair[0]] - xyzVals[pair[1]]), self.PQGEpsilon, self.PQGRm)

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

    def exportEllipsoid(self, filename):
        fIO.exportEllipsoid(self.blockNames, self.buildingBlockXYZ, self.minDist, self.minDist, fIO.fileRootFromInfile(filename, 'txt'))
            
if __name__ == "__main__":
    
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the generator objects.
    hairPinGen = spidroinHairpinGenerator(filename)
    
    # generate starting points and move the seed to those points extracting
    # the xyzVals each time. 
    pointA = np.array([10.0, 10.0, 10.0])
    pointB = np.array([0.0, 0.0, 0.0])
    
    minDist = 1.0
    numCrankMoves = 10
    
    # build building block and dump to file
    hairpinBuildingBlock = hairPinGen.generateBuildingBlock('SP1', pointA, pointB, minDist, numCrankMoves)
    hairpinBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "hairpin done"