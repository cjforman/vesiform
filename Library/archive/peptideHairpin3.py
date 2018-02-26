import sys
from Library.constrainedPolymer3 import ConstrainedPolymerPackBBG as CPBBG
import Utilities.fileIO as fIO
import numpy as np

class peptideHairpinGenerator(CPBBG):
    # A class for a buildingBlockGenerator object
    # that creates random peptide hairpins between two points.
    # Static values are acquired from a file
    # and dynamic values are applied at runtime.
    def __init__(self, paramFilename):
        # initialise the parameter dictionary for the base classes
        CPBBG.__init__(self, paramFilename)    
    
    def initialiseParameters(self):
        # initialise the constrained polymer parent
        CPBBG.initialiseParameters(self)
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for hairpin"
            sys.exit()        

    def generateBuildingBlock(  self, 
                                numPoints, 
                                pointA,
                                pointB, 
                                minDist, 
                                bondLength, 
                                polarity,
                                pointsToAvoid=[], 
                                envelope='innersphere'):
        
        self.polarity = polarity
        self.numResidues = int(np.ceil(numPoints /3.0))
        self.numPoints = self.numResidues * 3

        if self.numPoints != numPoints:
            print "Warning: numpoints changed to give integer number of residues."
        
        return CPBBG.generateBuildingBlock(self, numPoints, pointA, pointB, minDist, bondLength, pointsToAvoid=pointsToAvoid, envelope=envelope)

    def generateBuildingBlockNames(self):
        if (self.polarity == "NC"):
            names = ['N', 'C', 'C'] * self.numResidues
        if (self.polarity == "CN"):
            names = ['C', 'C', 'N'] * self.numResidues
        return names

        
if __name__ == "__main__":
    
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the backbone generator object.
    hairPinGen = peptideHairpinGenerator(filename)

    # generate a backbone
    numPoints = 100
    pointA = np.array([10, 1, 1])
    pointB = np.array([-10, -1, -1])
    minDist = 1.0
    bondLength = 2.0
    polarity = 'NC'

    # build building block and dump to file
    hairpinBuildingBlock = hairPinGen.generateBuildingBlock(numPoints, pointA, pointB, minDist, bondLength, polarity)
    hairpinBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "hairpin done"