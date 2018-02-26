import sys
import numpy as np
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
import Utilities.cartesian as cart
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.peptideHairpin import peptideHairpinGenerator as PHG
#from Utilities import cartesian as cart

class betasheetGen(BBG):
    # Class generates a beta sheet building block.  Two standard peptide strand building 
    # blocks are constructed with the peptideBackbone class. One is for parallel and
    # the other is for antiparallel. 
    #
    # These building blocks are then cloned and their reference points are placed along
    # a line through the origin, that is parallel with the crossStrandDirector. 
    # Each individual strand director is aligned with the user supplied inStrand director.
    # These allows for some interesting staggered geometries.
    # the strands are all identical unless the antiparallel parameter is set and in that
    # case every second strand has the opposite polarity to the first, and is rotated and 
    # offset by a user supplied amount. 
    #
    # Once the strands are constructed they are optionally connected with the peptideHairpin5 
    # class, which maintains the sequence of NCC while providing random hairpins that do not 
    # intersect either the strands or itself. The result is a single continuous peptide that
    # represents a beta sheet.
    #
    # IF the strands are internally connected with loops then the building block is
    # generated which has an N and a C terminal connector, with the N terminal first.
    # IF the strands are not connected then the ends of each strand are provided as 
    # connectors.  This allows for multiple beta sheets/strands to be connected as part 
    # of a much bigger chain.
     
    
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        
        self.betaStrandSeparation = self.getParam('betaStrandSeparation')
        self.CCbondLength = self.getParam('CCbondLength')
        self.CNbondLength = self.getParam('CNbondLength')        
        self.phi = self.getParam('phi') * np.pi / 180.0
        self.psi = self.getParam('psi') * np.pi / 180.0
        self.omega = self.getParam('omega') * np.pi / 180.0
        self.angleN = self.getParam('angleN') * np.pi / 180.0
        self.angleCA = self.getParam('angleCA') * np.pi / 180.0
        self.angleC = self.getParam('angleC') * np.pi / 180.0
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')

    def generateBuildingBlock(self, numStrands, strandLength, numLoopResidues, minDist, inStrandDirector, crossStrandDirector, offset, parallel=True, nameCA=False):
    
        self.startPos = np.array([0.0, 0.0, 0.0])
        self.numStrands = numStrands
        self.numPoints = strandLength * 3 # number of points per strand.
        self.strandLength = strandLength
        self.numLoopResidues = numLoopResidues
        self.minDist = minDist
        self.inStrandDirectorHat = inStrandDirector/np.linalg.norm(inStrandDirector)
        self.crossStrandDirectorHat = crossStrandDirector/np.linalg.norm(crossStrandDirector)
        self.oStrandDirectorHat = np.cross(self.inStrandDirectorHat, self.crossStrandDirectorHat)
        self.offsetXYZ = offset[2] * self.inStrandDirectorHat + offset[0] * self.crossStrandDirectorHat + offset[1] * self.oStrandDirectorHat
        self.parallel = parallel
        self.nameCA = nameCA
        
        return BBG.generateBuildingBlock(self, self.numPoints, minDist)
    
    def generateBuildingBlockXYZ(self):

        # create a building block generator using the standard parameter file for this object
        strandGen = PBG(self.paramFilename)

        # construct a prototype beta strand building block
        strandBB = strandGen.generateBuildingBlock(self.strandLength)
        
        # generate an array of strand building blocks
        BBs = [ strandBB.cloneBuildingBlock() for _ in range(0, self.numStrands)]
        baseLine = [ self.startPos + n * self.betaStrandSeparation * self.crossStrandDirectorHat for n in range(self.numStrands)]
        directors = [ self.inStrandDirectorHat for n in range(self.numStrands)]
        rotation= [ 0.0 for n in range(self.numStrands)]

        # modify the strands if we are anti-parallel        
        if not self.parallel:
            baseLine = [ basePoint if n % 2==0 else basePoint + self.offsetXYZ for n, basePoint in enumerate(baseLine) ]  
            directors = [ self.inStrandDirectorHat if n % 2==0 else -1 * self.inStrandDirectorHat for n in range(self.numStrands)]
            rotation = [ rot if n % 2==0 else 0 for n, rot in enumerate(rotation) ]  
        
        # do the sheet construction:
        [ BB.transformBBToLabFrame(director, pos, rot) for director, pos, rot, BB in zip(directors, baseLine, rotation, BBs)]
        
        # if we are adding loops then construct the loops
        if self.numLoopResidues > 0:
            Loops = self.constructLoops(BBs)
        
        # compile the final building block
        retBB = BBs[0]
        if self.numLoopResidues > 0:
            for BB, Loop in zip(BBs[1:], Loops):
                retBB.append(Loop)
                retBB.append(BB)
        else:
            [retBB.append(BB) for BB in BBs[1:]]

        return retBB.blockXYZVals
    
    def constructLoops(self, BBList):
        # create a loop generator object
        loopGen = PHG(self.paramFilename)
        
        # add all the beta strands to the points to avoid.
        pointsToAvoid= [] 
        for BB in BBList:
            pointsToAvoid += BB.blockXYZVals

        # generate an array to store loop building blocks
        Loops = []
        for BB1, BB2 in zip(BBList[0:self.numStrands-1], BBList[1:self.numStrands]):
            newLoop = self.constructLoop(BB1, BB2, loopGen, pointsToAvoid)
            Loops.append(newLoop)
            # add the new loop to points To Avoid.
            pointsToAvoid += newLoop.blockXYZVals
        return Loops

    def constructLoop(self, BB1, BB2, loopGen, pointsToAvoid):
        
        connectorA = BB1.getConnectionAtoms(1) # 1 is the C terminus of interest
        connectorB = BB2.getConnectionAtoms(0) # 0 is the N terminus of interest
        
        # generate TNB frames at either connector
        TNBA = coords.constructTNBFrame(connectorA[0], 
                                        connectorA[1], 
                                        connectorA[2])
        TNBB = coords.constructTNBFrame(connectorB[0], 
                                        connectorB[1], 
                                        connectorB[2])

        # C terminus is CCN triad (new atom will be an N).  
        betaA =  self.angleC
        alphaA = self.phi
        # N Terminus is CNC triad. (new atom will be a C).
        betaB =  self.angleN
        alphaB = self.psi
        
        # compute the coil start and end points
        pointA = connectorA[2] + self.CNbondLength * coords.generateTNBVecXYZ(TNBA, betaA, alphaA) 
        pointB = connectorB[2] + self.CNbondLength * coords.generateTNBVecXYZ(TNBB, betaB, alphaB)

        if self.dumpInterimFiles:
            fIO.saveXYZList([pointA, connectorA[2], pointB, connectorB[2]], ['Ca', 'S', 'O', 'S'], 'connectionDetails.xyz')
            fIO.saveXYZ(pointsToAvoid, 'K', 'pointsToAvoid.xyz')
        
        numCrankMoves = 0
        if self.parallel:
            iSphereR = 0.4 * self.strandLength * 3.5
        else:
            iSphereR = 0.9 * self.betaStrandSeparation/2.0
       
        # create the loop building Blockss
        return loopGen.generateBuildingBlock( self.numLoopResidues,
                                              pointA,
                                              pointB, 
                                              self.minDist,
                                              numCrankMoves, 
                                              pointsToAvoid = pointsToAvoid, 
                                              envelopeList = ["innersphere " + str(iSphereR)])  

    def generateBuildingBlockNames(self):
        totalNumResidues = len(self.buildingBlockXYZ)/3
        if self.nameCA:
            names = ['N', 'CA', 'C'] * totalNumResidues
        else:
            names = ['C', 'C', 'N'] * totalNumResidues
        return names

    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)

    def generateBuildingBlockDirector(self):
        return self.inStrandDirectorHat  
    
    def generateBuildingBlockConnectors(self):
        # 0th connector is N terminal.
        # last connector is C terminal
        l = len(self.buildingBlockXYZ)
        connectors = [np.array([2, 1, 0]), np.array([l - 3, l - 2, l - 1])]
        
        # if there are no loops then add all the end points as connectors
        if self.numLoopResidues==0:
            connectors=[]
            for strand in range(0, self.numStrands):
                l = 3 * self.strandLength # num points in a strand
                n = l * strand # total number of points so far
                
                # don't insert the connector if it's already there.
                connectors.insert(-1, np.array([2 + n, 1 + n, 0 + n]) )
                connectors.insert(-1, np.array([l - 3 + n, l - 2 + n, l - 1 + n]))

        return connectors

if __name__=="__main__":
    # get the file name from the command line
    filename = sys.argv[1]

    # create the backbone generator object using static file parameters
    betasheetGenerator = betasheetGen(filename)

    # generate backbone realtime parameters
    numStrands = 4
    lengthStrand = 7
    numLoopResidues = 0
    minDist = 1.0
    inStrandDirector = np.array([0.0, 0.0, 1.0])
    crossStrandDirector = np.array([1.0, 0.0, 0.0])
    offset =  np.array([0.0, 0.0, 0.0])
    
    betasheetBB = betasheetGenerator.generateBuildingBlock(  numStrands, 
                                                             lengthStrand,
                                                             numLoopResidues,
                                                             minDist,
                                                             inStrandDirector, 
                                                             crossStrandDirector, 
                                                             offset, 
                                                             parallel = True) 

    # global positioning
    startPos = np.array([0.0, 0.0, 0.0])
    globalDirector = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    betasheetBB.transformBBToLabFrame(globalDirector, startPos, rotation)
    betasheetBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print "betasheet done"
    
    
        
    