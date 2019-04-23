'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB
import sys
import numpy as np
import random as rnd
from Utilities import coordSystems as coords
from Utilities import fileIO as fIO

class SurfacePackCylinderBBG(NBB):
    '''
    Over loads the NPack object to give a packing on a cylinder where no point is less than minDist away from the other.
    The cylinder is parallel with the z axis and has a defined radius and length. 
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for Cylinder Packing")
            sys.exit()        

    def generateBuildingBlock( self, 
                               numPoints, 
                               rx, 
                               ry, 
                               z1,
                               z2, 
                               phi1, 
                               phi2, 
                               minDist,
                               envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):

        self.radiusX = rx
        self.radiusY = ry
        
        # process the zs
        self.z1 = min(z1, z2)
        self.z2 = max(z1, z2)
        
        # compute the length of the cylinder
        self.zLength = abs(self.z2 - self.z1) 
        
        self.phi1 = phi1 * np.pi/180
        self.phi2 = phi2 * np.pi/180
        if self.phi2 < self.phi1:
            self.phi2 = self.phi1
            self.phi1 = phi2 * np.pi/180.0
        
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockDirector(self):
        return np.array([0,0,1])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, self.z1])
        
    def pickFirstPoints(self):
        
        # select the first point in the NList and check it against the allowed bounds
        
        # assume that the point is out of bounds
        inBounds = False
        numLoops = 0
        # loop while the point is out of bounds or we hit the max number of loops
        while(inBounds==False and numLoops< self.maxLivesPerNewNode): 
            pos = self.pickRandomPointInDefinedSpace()
            inBounds = self.checkPointInBounds(pos)
            numLoops += 1

        # if we found a spot then convert to xyz and initialise list
        if inBounds==True:
            self.nList = [self.convertPointToXYZ(pos)]                        
        else:
            # if we didn't find a spot then hack a start spot right into the centre of the allowed bounds
            self.nList = [self.convertPointToXYZ( np.array( [ coords.azToRhoCyl( ( self.phi1 + self.phi2 ) / 2.0, 
                                                                                   self.radiusX, 
                                                                                   self.radiusY), 
                                                              ( self.phi1 + self.phi2)/2.0 , 
                                                              ( self.z2 - self.z1)/2 ]) ) ]
            inBounds = True
        self.nAttempts = [0]            
    
        return inBounds
    
    def checkPointInBounds(self, pos):
        # cylindrical coords format for pos. Return true if point is within bounds and false if it is out of bounds.
        
        # assume point is inBounds and test for outBoundness
        outBounds = False
        
        if pos[2] < self.z1:
            outBounds=True
            if self.verbose==1: 
                print("z1 violation")             
        
        if pos[2] > self.z2:
            outBounds = True
            if self.verbose==1: 
                print("z2 violation")             
        
        if pos[1] < self.phi1:
            outBounds = True
            if self.verbose==1: 
                print("phi1 violation")             
        
        if pos[1] > self.phi2:
            outBounds = True
            if self.verbose==1: 
                print("phi2 violation")             

        if (pos[0] < self.radiusX) and (pos[0] < self.radiusY):
            outBounds = True
            if self.verbose==1: 
                print("radius inner violation")             
        
        if (pos[0] > self.radiusX) and (pos[0] > self.radiusY):
            outBounds = True
            if self.verbose==1: 
                print("radius outer violation")
            
        if outBounds == True:
            inBounds = False
        else:
            inBounds = True
            
        if inBounds:
            inBounds = NBB.checkPointInBounds(self, pos)
            
        return inBounds
             
    def pickRandomPointInDefinedSpace(self):
        phi = coords.pickRandomPointOnUnitCircle()
        z = rnd.uniform(self.z1, self.z2)
        # compute the rho at the chosen phi
        rho = coords.azToRhoCyl(phi, self.radiusX, self.radiusY)
        return np.array([rho, phi, z]) # pos is in cylindrical co-ordinates
        
    def convertPointToXYZ(self, pos):
        # this function takes a numpy array ([rho, az, z]) and returns an XYZ numpy array.
        return coords.cyl2XYZ(pos)
        
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    CylinderBBG = SurfacePackCylinderBBG(filename)

    numPoints = 300
    centerPos = np.array([-0,-0, 0])
    director = np.array([0, 0, 1])
    rotation = 0
    rx = 15
    ry = 10
    z1 = 0
    z2 = 40    
    phi1 = -135
    phi2 = 135
    minDist = 2
       
    # generate the XYZVals in the packed spaced
    CylinderBB = CylinderBBG.generateBuildingBlock(numPoints, rx, ry, z1, z2, phi1, phi2, minDist)
    CylinderBB.transformBBToLabFrame(director, centerPos, rotation)
    CylinderBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    