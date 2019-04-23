'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB
import sys
import numpy as np
from Utilities import coordSystems as coords

class SurfacePackSphereBBG(NBB):
    '''
    Over loads the NPackBB object to give a packing of numPoints particles on a sphere where 
    no point is less than minDist away from the other (euclidean distance - chord of sphere). 
    The user can specify radius, and angular range over which the sphere is populated. 
    Thes angles are defined with using standard polar coords (theta -elevation, and phi is azimuth)
    After construction the z-axis is rotated to be parallel with the specified 
    director and the object rotated about that axis by rotation degrees. The centre of
    the sphere is placed at centrePos.  
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for SpherePack Object")
            sys.exit()        

    def generateBuildingBlock(self, numPoints, radius, theta1, theta2, phi1, phi2, minDist, envelopeList=['None'], pointsToAvoid=[], showBlockDirector=False, visualiseEnvelope=(0,200)):
    
        self.radius = radius
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
        
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, showBlockDirector=showBlockDirector, visualiseEnvelope=(0,200))
    
    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def pickFirstPoints(self):
        inbounds=False
        while inbounds==False:
            # select the first point in the NList 
            newPoint = self.pickRandomPointInDefinedSpace()
            inbounds = self.checkPointInBounds(newPoint)
        
        self.nList = [self.convertPointToXYZ(newPoint)]
        self.nAttempts = [0]             
        return True
    
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
        
        if inBounds:
            inBounds = NBB.checkPointInBounds(self, pos)
        
        return inBounds
             
    def pickRandomPointInDefinedSpace(self):
        # returns spherical polar coords for a randomly picked point on a sphere
        theta, phi = coords.pickRandomPointOnUnitSphere()
        return np.array([self.radius, theta, phi])
        
    def convertPointToXYZ(self, pos):
        # this function takes a numpy array ([r, theta, phi]) and returns an XYZ numpy array.
        return coords.sphericalPolar2XYZ(pos)
        
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPackBB object.
    SphereNBB = SurfacePackSphereBBG(filename)
    numPoints = 70
    centerPos = np.array([-0,-0, -0])
    director = np.array([1, 1, 1])
    rotation = 0
    radius = 20.0
    theta1 = -90.0
    theta2 = 90.0
    phi1 = 45.0
    phi2 = 90.0
    minDist = 2.5
    
    # generate the XYZVals in the packed spaced
    SphereBB = SphereNBB.generateBuildingBlock(numPoints, radius, theta1, theta2, phi1, phi2, minDist)
    SphereBB.transformBBToLabFrame(director, centerPos, rotation)
    SphereBB.exportBBK("sphere1")
    
   # Sphere2BB = SphereNBB.generateBuildingBlock(numPoints, radius + 2 * minDist, theta1, theta2, phi1, phi2, minDist)
   # Sphere2BB.transformBBToLabFrame(director, centerPos, rotation)
   # Sphere2BB.exportBBK("sphere2")
    
   # Sphere3BB = SphereNBB.generateBuildingBlock(numPoints, radius + 4 * minDist, theta1, theta2, phi1, phi2, minDist)
   # Sphere3BB.transformBBToLabFrame(director, centerPos, rotation)
   # Sphere3BB.exportBBK("sphere3")
    
    
    
    print("Sphere Pack Done")