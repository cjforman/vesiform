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
            print("Critical Parameters are undefined for NPackBBG object")
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, xRadius, yRadius, zRadius, theta1, theta2, phi1, phi2, minDist, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200,'envelope.xyz'), showBlockDirector=False, defaultBlockRefPoint=None):
        self.xRadius = xRadius
        self.yRadius = yRadius
        self.zRadius = zRadius
        self.theta1 = theta1 * np.pi/180.0
        self.theta2 = theta2 * np.pi/180.0
        self.phi1 = phi1 * np.pi/180.0
        self.phi2 = phi2 * np.pi/180.0
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector, defaultBlockRefPoint=defaultBlockRefPoint)
    
    def generateBuildingBlockDirector(self):
        director = np.array([0.0, 0.0, 1.0])
        return director

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])
    
    def pickFirstPoints(self):
        inbounds = False
        while inbounds==False:
            newPoint = self.pickRandomPointInDefinedSpace()
            inbounds = self.checkPointInBounds(newPoint)
        self.nList = [self.convertPointToXYZ(newPoint)]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        return coords.pickRandomPointInEllipsoidRange(self.xRadius, self.yRadius, self.zRadius, self.theta1, self.theta2, self.phi1, self.phi2)
       
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    EllipsoidPackBBG = VolumePackEllipsoidBBG(filename)

    numPoints = 300
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 0 
    xRadius = 100
    yRadius = 100
    zRadius = 100
    minDist = 1
    theta1 = -90.0
    theta2 = 90.0
    phi1 = 45.0
    phi2 = 90.0
    
    envelopeList = ['None']
    envelopeList = ['innersphere 60',
                    'outersphere 80']       
    # generate the building block
    EllipsoidPackBB = EllipsoidPackBBG.generateBuildingBlock(numPoints, xRadius, yRadius, zRadius, theta1, theta2, phi1, phi2, minDist, envelopeList = envelopeList, defaultBlockRefPoint=centrePos)
    EllipsoidPackBB.transformBBToLabFrame(director, centrePos, rotation)
    EllipsoidPackBB.blockAtomNames[0]='P'
    EllipsoidPackBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    