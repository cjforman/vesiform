import sys
import numpy as np
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
import Utilities.fileIO as fIO
from Library.peptideHairpin import peptideHairpinGenerator as PHG
from Library.peptideBackbone import peptideBackboneGenerator as PBG

class chainOfChainsGenerator(BBG):
    # Class generates a peptide hairpin of a specified length between pairs of input points in a list. 
    # Each hairpin is contrained in a spherocylinder along the line connecting the points with a diameter
    # specified in the diameter array.  
    
    def __init__(self, filename):
        self.filename = filename
        BBG.__init__(self, filename)


    def initialiseParameters(self):
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')
        self.PHG = PHG(self.filename)
        self.PBG = PBG(self.filename)
        
    def generateBuildingBlock(self, numResidues, pointsA, pointsB, radii, minDist, visualiseEnvelope = (0,20,'envelope.xyz')):
    
        self.numSegments = len(numResidues)
        if len(pointsA) != self.numSegments:
            print "Num supplied points does not match num Segments in chainOfChains"
        if len(pointsB) != self.numSegments:
            print "Num supplied points does not match num Segments in chainOfChains"
        if len(radii) != self.numSegments:
            print "Num supplied diameters does not match num Segments in chainOfChains"

        self.pointsA = pointsA
        self.pointsB = pointsB
        self.numResidues = numResidues
        self.radii = radii
        self.numPoints = numResidues * self.numSegments * 3
        self.minDist = minDist
        self.envSize = visualiseEnvelope[0]
        self.envRange = visualiseEnvelope[1]
        
        return BBG.generateBuildingBlock(self, self.numPoints, self.minDist)
    
    def generateBuildingBlockXYZ(self):
        
        
        # generate 1 residue seed building block
        seedResidue = self.PBG.generateBuildingBlock(1)
    
        xyzVals = []
        namesTemp = []
        chainNum = 0
        for numResidues, pointA, pointB, radius in zip(self.numResidues, self.pointsA, self.pointsB, self.radii):
            
            seedResidue.placeAtom(2, pointA)
            seedResidue.setBlockRefPoint(pointA)
            seedResidue.orientToDirector(pointA - pointB)
            pointsA = seedResidue.blockXYZVals[:]
            
            seedResidue.placeAtom(2, pointB)
            seedResidue.setBlockRefPoint(pointB)
            seedResidue.orientToDirector(np.array([0.0, 1.0, 1.0]))
            seedResidue.orientToDirector(pointB - pointA)
            
            pointsB = seedResidue.blockXYZVals[:]
            
            if self.dumpInterimFiles==1:
                fIO.saveXYZList(pointsA+pointsB, ['B']*3 + ['P'] * 3, "seedPoints.xyz")
            
            envelopeList = [ "Spherocylinder " + str(pointA[0]) + " " + str(pointA[1]) + " " + str(pointA[2]) + " " + str(pointB[0]) + " " + str(pointB[1]) + " " + str(pointB[2]) + " " + str(radius) ] 
            PHP = self.PHG.generateBuildingBlock( numResidues, 
                                                  pointsA,
                                                  pointsB,
                                                  self.minDist,
                                                  0,
                                                  visualiseEnvelope=(self.envSize, self.envRange, 'envelope'+str(chainNum)+'.xyz'),
                                                  envelopeList=envelopeList) 
            
            if len(xyzVals)==0:
                xyzVals = PHP.blockXYZVals
                namesTemp= PHP.blockAtomNames
            else:
                xyzVals = np.concatenate((xyzVals, PHP.blockXYZVals), 0)
                namesTemp = np.concatenate((namesTemp, PHP.blockAtomNames), 0)
            self.namesTemp = namesTemp
            chainNum += 1

        return xyzVals 

    def generateBuildingBlockNames(self):
        return self.namesTemp

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])  
    
if __name__=="__main__":
    # get the file name from the command line
    filename = sys.argv[1]

    # create the backbone generator object using static file parameters
    chainOfChainsGen = chainOfChainsGenerator(filename)

    # generate backbone realtime parameters
    minDist = 1.0
    residues = [7, 8, 8, 12]
    radii = [minDist * 2, minDist * 3, minDist * 3, minDist * 4]
    pointsA = [ np.array([0.0, 0.0, 0.0]),
                np.array([25.0, 0.0, 0.0]),
                np.array([25.0, 25.0, 0.0]),
                np.array([25.0, 25.0, 25.0])]
    
    pointsB = [ np.array([23.5, 0.0, 0.0]),
                np.array([25.0, 23.5, 0.0]),
                np.array([25.0, 25.0, 23.5]),
                np.array([0.0, 0.0, 25.0]) ]
    
    
    fIO.saveXYZList(pointsA + pointsB, ['Ca'] * len(pointsA) + ['O'] * len(pointsB), "labPoints.xyz")
    
    cOfChainsBB = chainOfChainsGen.generateBuildingBlock(residues, pointsA, pointsB, radii, minDist, visualiseEnvelope=(100000, 50, 'envelope.xyz'))

    cOfChainsBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print "chainOfChains done"
    
    
        
    