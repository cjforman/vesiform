import sys
from Library.constrainedPolymer import ConstrainedPolymerPackNBB as CPNBB
import Utilities.fileIO as fIO
import numpy as np

class peptideHairpinGenerator(CPNBB):
    # A class for a buildingBlockGenerator object
    # that creates random peptide hairpins between two points.
    # Static values are acquired from a file
    # and dynamic values are applied at runtime.
    def __init__(self, paramFilename):
        # initialise the parameter dictionary for the base classes
        # and create an NPack constrained polymer object
        # use the same filename for both so that variables common
        # to both have the same value
        CPNBB.__init__(self, paramFilename)    
    
    def initialiseParameters(self):
        # initialise the constrained polymer parent
        CPNBB.initialiseParameters(self)
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for hairpin"
            sys.exit()        

    def generateBuildingBlock(self, 
        numPoints, 
        pointA, 
        pointB, 
        alpha1, 
        alpha2, 
        beta1, 
        beta2, 
        minDist, 
        bondLength, 
        innerSphereR, 
        outerSphereR, 
        spherePos,
        polarity):
        
        self.polarity = polarity
        self.numResidues = int(np.ceil(numPoints /3.0))
        self.numPoints = self.numResidues * 3
        
        if self.numPoints != numPoints:
            print "Warning: numpoints changed to give integer number of residues."
        
        return CPNBB.generateBuildingBlock(self, numPoints, pointA, pointB, rotation, alpha1, alpha2, beta1, beta2, minDist, bondLength, innerSphereR, outerSphereR, spherePos)

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
    numPoints = 18
    pointA = np.array([-5, 0, 0]) 
    pointB = np.array([ 5, 0, 0])
    rotation = 90
    alpha1 = -180
    alpha2 = 180
    beta1 = 130
    beta2 = 160
    minDist = 1.0
    bondLength = 2.0
    innerSphereR = 0
    outerSphereR = 100
    spherePos = np.array([0, 0, -5]) 
    polarity = 'NC'

    # build building block and dump to file
    hairpinBuildingBlock = hairPinGen.generateBuildingBlock(numPoints, pointA, pointB, rotation, alpha1, alpha2, beta1, beta2, minDist, bondLength, innerSphereR, outerSphereR, spherePos, polarity)
    hairpinBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "hairpin done"