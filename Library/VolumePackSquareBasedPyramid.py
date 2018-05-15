'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import random as rnd
import numpy as np
from Utilities import fileIO as fIO

class VolumePackSquareBasedPyramidBBG(NBB):
    '''
    Over loads the NPackBB object to give a uniform packing throughout a 
    rectangular based pyramid so that no point is less than minDist away from 
    the other in the specified region of space. 
    The origin is at the centre of the square base and the pyramid central
    axis if parallel with the z-axis.
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for SpacePack object"
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, xBase, yBase, zHeight, zHeightApex, minDist, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):
        self.xBase = np.abs(xBase) 
        self.yBase = np.abs(yBase)
        self.zHeight = np.abs(zHeight) # the Z Height above which no point will appear (allow a frustum shaped object) 
        self.zHeightApex = np.abs(zHeightApex) # Sets the slope of the sides of the pyramid by specified the height of the pyramid above the centre of the base. 
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=False)
    
    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])
    
    def pickFirstPoints(self):
        self.nList = [self.pickRandomPointInDefinedSpace()]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        # pick a random point on the base of the pyramid - only using half the pyramid in y.
        x = rnd.uniform( -self.xBase/2.0, self.xBase/2.0 )
        y = rnd.uniform( 0, self.yBase/2.0 )
        
        # pick a random point in the volume above the base that is below the max allowed altitude.
        z = rnd.uniform( 0, self.zHeight)
        
        # keep choosing a random point above the base until we get one 
        # that is lower than the roof of the pyramid at that point.
        while z >= self.computeHeightOfPyramidAtPoint(x, y, self.xBase, self.yBase, self.zHeightApex):
            z = rnd.uniform( 0, self.zHeight)
        
        return np.array( [x, y, z] )

    def computeHeightOfPyramidAtPoint(self, x, y, XB, YB, ZA):
        # returns the height of a square based pyramid at a point x or y in the base of the pyramid.
        # the pyramid is defined by the height ZA of the apex above the base.
        # can use the absolute value of x and y without loss og generality due to symmetry of the pyramid/
        x = np.abs(x)
        y = np.abs(y)
         
        # first step is to decide if we are under the +/- X slope of the +/-Y slope.
        # this line in the roof projects down to the eqn y = x * yBase/xBase.
        # so compute the y value on the projected line at the given x coordinate.
        # if y is greater than this value we are under one slope, 
        # if y is less we are under the other slope.
        # on the line we can use either slope to calculate z height.  
        if y >= x * YB/XB:
            z = ZA * (1.0 - (2.0 * y / YB) ) # from simple trigonometry
        else:
            z = ZA * (1.0 - (2.0 * x / XB) )
        return z

    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    PyramidPackBBG = VolumePackSquareBasedPyramidBBG(filename)

    numPoints = 550
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 0 
    xBase = 30.0
    yBase = 60.0
    zHeight = 50.0
    zHeightApex = 70.0
    minDist = 3.0
       
    # generate the building block
    SpacePackBB = PyramidPackBBG.generateBuildingBlock(numPoints, xBase, yBase, zHeight, zHeightApex, minDist)
    SpacePackBB.transformBBToLabFrame(director, centrePos, rotation)
    SpacePackBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    