'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import random as rnd
import numpy as np
from Utilities import fileIO as fIO

class SurfacePackPlaneBBG(NBB):
    '''
    Over loads the NPackBB object to give a packing on a flat plane where no point is 
    less than minDist away from the other in the specified X and Y ranges. 
    During construction the plane is parallel with 
    the XY plane and is at a distance distFromOrig from the origin at the z axis.
    The plane can then be rotated so the z-axis is parallel with director and
    then rotated about director by rotation degrees. The centre is placed at center Point
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        self.blockRefPoint = np.array([0.0, 0.0, 0.0])
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for PlanePack object")
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, xRange1, xRange2, yRange1, yRange2, distFromOrig, minDist, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):
        self.distFromOrig = distFromOrig
        self.xRange1 = xRange1 
        self.xRange2 = xRange2
        self.yRange1 = yRange1
        self.yRange2 = yRange2
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockDirector(self):
        return np.array([0,0,1])

    def generateBuildingBlockRefPoint(self):
        print("holeyMonkeys")
        return np.array([(self.xRange1 + self.xRange2)/2, 
                         (self.yRange1 + self.yRange2)/2,
                         self.distFromOrig])
    
    def pickFirstPoints(self):
        self.nList = [np.array([rnd.uniform(self.xRange1, self.xRange2), rnd.uniform(self.yRange1, self.yRange2), self.distFromOrig])]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        return np.array( [ rnd.uniform( self.xRange1, self.xRange2 ), 
                           rnd.uniform( self.yRange1, self.yRange2 ),
                           self.distFromOrig ] )
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    PlanePackBBG = SurfacePackPlaneBBG(filename)

    numPoints = 100
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 45 
    xRange1 = -10
    xRange2 = 10
    yRange1 = -10
    yRange2 = 10
    distFromOrig = 0
    minDist = 1
       
    # generate the building block
    PlanePackBB = PlanePackBBG.generateBuildingBlock(numPoints, xRange1, xRange2, yRange1, yRange2, distFromOrig, minDist)
    PlanePackBB.transformBBToLabFrame(director, centrePos, rotation)
    PlanePackBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    