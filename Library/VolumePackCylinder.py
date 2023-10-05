'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import numpy as np
from Utilities import fileIO as fIO
from Utilities import coordSystems as coords

class VolumePackCylinderBBG(NBB):
    '''
    Over loads the NPackBB object to give a uniform packing throughout a 3Space 
    so that no point is less than minDist away from the other in a 
    cylindrical region of space. The centre of the space is the center of the base cylinder.
    The axis of the cylinder is aligned with the z axis.  
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for PlanePack object")
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, r, h, minDist, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0, 200, 'envelope.xyz'), showBlockDirector=False):
        self.r = float(r) 
        self.h = float(h)
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=False)
    
    def generateBuildingBlockDirector(self):
        return np.array([0,0,1])

    def generateBuildingBlockRefPoint(self):
        return np.array([0, 
                         0,
                         0])
    
    def pickFirstPoints(self):
        self.nList = [self.pickRandomPointInDefinedSpace()]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        (rho, phi, z) = coords.pickRandomPointInUnitCylinder()
        return coords.cyl2XYZ(np.array([self.r * rho, phi, self.h * z]))
        
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    CylinderPackBBG = VolumePackCylinderBBG(filename)

    numPoints = 300
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 0 
    r = 10
    h = 20
    minDist = 2
       
    # generate the building block
    SpacePackBB = CylinderPackBBG.generateBuildingBlock(numPoints, r,h, minDist)
    SpacePackBB.transformBBToLabFrame(director, centrePos, rotation)
    SpacePackBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print("Done")
    