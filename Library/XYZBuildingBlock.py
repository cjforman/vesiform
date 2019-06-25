import sys
from Utilities.keyProc import keyProc
import Utilities.cartesian as cart
import Utilities.fileIO as fIO
import numpy as np
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG


class XYZBuildingBlockGenerator(BBG):
    # This is a class for generating a building block from an XYZ file. 
    #
    # There is no input at object creation but at generation time one must supply
    # a valid xyz file.
    
    def __init__(self):
        # This is a keyproc object because a BuildingBlockGenerator is 
        # a keyproc object. However, object has no parameters so no need 
        # to invoke the keyproc part of the object
        # keyProc.__init__(self, paramFilename)
        self.dumpInterimFiles = 0
        self.director = []
        pass
        
    def initialiseParameters(self):
        
        keyProc.initialiseParameters(self)
    
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for XYZBuildingBlockGenerator")
            sys.exit()        

    def generateBuildingBlock(self, XYZFileName, showBlockDirector=False, visualiseEnvelope=(0,20, 'envelope.xyz'), envelopeList=['None'], pointsToAvoid=[]):
        self.atomNameList, self.xyzList = fIO.loadXYZ(XYZFileName)
        self.numPoints = len(self.atomNameList)
        minDistDummy = 0.0
        return BBG.generateBuildingBlock(self, self.numPoints, minDistDummy, showBlockDirector=showBlockDirector)

    def generateBuildingBlockXYZ(self):
        return self.xyzList   
        
    def generateBuildingBlockNames(self):
        return self.atomNameList
    
    def generateBuildingBlockDirector(self):
        # check to see if there was an externally imposed director
        director = None
        if self.director==[]:
            # calculate the eigenvectors of the inertial tensor of the xyz file and 
            # return the unit vector along the principal axis ( assumes each point is of unit mass so not quite correct).
            pAxis = cart.getPrincipalAxis(self.buildingBlockXYZ)
            director = pAxis/np.linalg.norm(pAxis)
        else:
            director = self.director

        return director

    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)
        
if __name__ == "__main__":
    

    # create the pdb object generator
    xyzObject = XYZBuildingBlockGenerator()

    # generate backbone realtime parameters
    refPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    
    xyzBuildingBlock = xyzObject.generateBuildingBlock('XYZFile.xyz', showBlockDirector=True) 
    xyzBuildingBlock.transformBBToLabFrame(director, refPos, rotation)
    xyzBuildingBlock.exportBBK('XYZFile_out') # adds the xyz extension automatically
    
    print("xyz building block done")