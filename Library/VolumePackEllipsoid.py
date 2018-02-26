'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import numpy as np
import Utilities.coordSystems as coords 
from Utilities import fileIO as fIO

class VolumePackEllipsoidBBG(NBB):
    '''
    Over loads the NPack object to give a uniform packing throughout an ellipsoidal 3Space 
    so that no point is less than minDist away from the other in the 
    cubic region of space. The centre of the space is the origin. 
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for NPackBBG object"
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, xRadius, yRadius, zRadius, minDist, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):
        self.xRadius = xRadius
        self.yRadius = yRadius
        self.zRadius = zRadius
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockDirector(self):
        director = np.array([0.0, 0.0, 1.0])
        return director

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])
    
    def pickFirstPoints(self):
        self.nList = [self.pickRandomPointInDefinedSpace()]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        return coords.pickRandomPointInEllipsoid(self.xRadius, self.yRadius, self.zRadius)
       
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    EllipsoidPackBBG = VolumePackEllipsoidBBG(filename)

    numPoints = 600
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 0 
    xRadius = 10
    yRadius = 10
    zRadius = 10
    minDist = 1
    
    envelopeList = ['None']
       
    # generate the building block
    EllipsoidPackBB = EllipsoidPackBBG.generateBuildingBlock(numPoints, xRadius, yRadius, zRadius, minDist, envelopeList = envelopeList)
    EllipsoidPackBB.transformBBToLabFrame(director, centrePos, rotation)
    EllipsoidPackBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    