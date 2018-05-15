'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import random as rnd
import numpy as np
from Utilities import fileIO as fIO

class VolumePackCuboidBBG(NBB):
    '''
    Over loads the NPackBB object to give a uniform packing throughout a 3Space 
    so that no point is less than minDist away from the other in the 
    cubic region of space. The centre of the space is the origin. 
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for PlanePack object"
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, xRange1, xRange2, yRange1, yRange2, zRange1, zRange2, minDist, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):
        self.xRange1 = xRange1 
        self.xRange2 = xRange2
        self.yRange1 = yRange1
        self.yRange2 = yRange2
        self.zRange1 = zRange1
        self.zRange2 = zRange2
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=False)
    
    def generateBuildingBlockDirector(self):
        return np.array([0,0,1])

    def generateBuildingBlockRefPoint(self):
        return np.array([(self.xRange1 + self.xRange2)/2, 
                         (self.yRange1 + self.yRange2)/2,
                         (self.zRange1 + self.zRange2)/2])
    
    def pickFirstPoints(self):
        self.nList = [np.array([rnd.uniform(self.xRange1, self.xRange2), rnd.uniform(self.yRange1, self.yRange2), rnd.uniform(self.zRange1, self.zRange2)])]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        return np.array( [ rnd.uniform( self.xRange1, self.xRange2 ), 
                           rnd.uniform( self.yRange1, self.yRange2 ),
                           rnd.uniform( self.zRange1, self.zRange2 )] ) 

    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    PlanePackBBG = VolumePackCuboidBBG(filename)

    numPoints = 300
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 0 
    xRange1 = -10
    xRange2 = 10
    yRange1 = -10
    yRange2 = 10
    zRange1= -10
    zRange2= 10
    minDist = 2
       
    # generate the building block
    SpacePackBB = PlanePackBBG.generateBuildingBlock(numPoints, xRange1, xRange2, yRange1, yRange2, zRange1, zRange2, minDist)
    SpacePackBB.transformBBToLabFrame(director, centrePos, rotation)
    SpacePackBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    