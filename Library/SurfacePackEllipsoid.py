'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB
import sys
import numpy as np
from Utilities import coordSystems as coords
from Utilities import fileIO as fIO

class SurfacePackEllipsoidBBG(NBB):
    '''
    Over loads the NPack object to give a packing on an ellipsoid where no point is less than minDist away from the other.
    The ellispoids semi minor axes are aligned with the principle lab axes. Can restrain quadrant in which ellipsoid is populated. 
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        # add ellipsoidal radii - sort so that c > b > a for z, y and x axes accordingly.
        # if we have to flip axes, remember and swap the vector components at the end. 
        
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Ellipsoid Pack object"
            sys.exit()        

    def generateBuildingBlock( self, 
                               numPoints, 
                               rx, 
                               ry, 
                               rz, 
                               theta1, 
                               theta2, 
                               phi1, 
                               phi2, 
                               minDist,
                               envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):
        self.radiusX = rx
        self.radiusY = ry
        self.radiusZ = rz
        self.theta1 = theta1 * np.pi/180
        self.theta2 = theta2 * np.pi/180
        if self.theta2 < self.theta1:
            self.theta2 = self.theta1
            self.theta1 = theta2 * np.pi/180.0
        self.phi1 = phi1 * np.pi/180
        self.phi2 = phi2 * np.pi/180
        if self.phi2 < self.phi1:
            self.phi2 = self.phi1
            self.phi1 = phi2 * np.pi/180.0
        
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockDirector(self):
        return np.array([0,0,1])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])
        
    def pickFirstPoints(self):
        # select the first point in the NList 
        inBounds = False
        numLoops = 0

        while not inBounds and numLoops<self.maxLivesPerNewNode:
            pos = self.pickRandomPointInDefinedSpace()
            inBounds = self.checkPointInBounds(pos)
            numLoops += 1

        if inBounds:
            posXYZ = self.convertPointToXYZ(pos)
            self.nList = [posXYZ]     
            self.nAttempts = [0]
        
        return inBounds  
    
    def checkPointInBounds(self, pos):
        # spherical coords format for pos. (r, theta, phi) 
        # Return true if point is within bounds and false if it is out of bounds.
        
        outBounds = False
        
        if pos[2] <= self.phi1:
            outBounds = True
            if self.verbose==1: 
                print("phi1 violation")             
        
        if pos[2] >= self.phi2:
            outBounds = True
            if self.verbose==1: 
                print("phi2 violation")             

        if pos[1] <= self.theta1:
            outBounds = True
            if self.verbose==1: 
                print("theta1 violation")             
        
        if pos[1] >= self.theta2:
            outBounds = True
            if self.verbose==1: 
                print("theta2 violation")             

        if outBounds == True:
            inBounds = False
        else:
            inBounds = True
        return inBounds
             
    def pickRandomPointInDefinedSpace(self):
        # returns spherical polar coords for a randomly picked point on the surface ellipsoid.
        return coords.pickRandomPointUniformlyOnEllipsoid(self.radiusX, self.radiusY, self.radiusZ)
        
    def convertPointToXYZ(self, pos):
        # this function takes a numpy array ([r, theta, phi]) and returns an XYZ numpy array.
        return coords.sphericalPolar2XYZ(pos)
        
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPackBB object.
    EllipsoidNBB = SurfacePackEllipsoidBBG(filename)
    numPoints = 150
    centerPos = np.array([-0,-0, 0])
    director = np.array([0, 0, 1])
    rotation = 00
    rx = 20
    ry = 10
    rz = 8
    theta1 = -90
    theta2 = 90
    phi1 = -90
    phi2 = 90
    minDist = 2
    
    # generate the XYZVals in the packed spaced
    EllipsoidBB = EllipsoidNBB.generateBuildingBlock(numPoints, rx, ry, rz, theta1, theta2, phi1, phi2, minDist)
    EllipsoidBB.transformBBToLabFrame(director, centerPos, rotation)
    EllipsoidBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))