'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import numpy as np
import Utilities.coordSystems as coords 
import Utilities.cartesian as cart
from Utilities import fileIO as fIO

class VolumePackEllipsoidSCParticlesBBG(NBB):
    '''
    Over loads the NPack object to give a uniform packing throughout an ellipsoidal 3Space 
    with spherocylindrical particles, such that no particle intersects 
    the others in the elliptical region of space. The centre of the space is the origin. 
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for NPackBBG object"
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, xRadius, yRadius, zRadius, theta1, theta2, phi1, phi2, thetaDir1, thetaDir2, phiDir1, phiDir2, minDist, pLength, numSpheresPerParticle, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200,'envelope.xyz'), showBlockDirector=False):
        self.xRadius = xRadius
        self.yRadius = yRadius
        self.zRadius = zRadius
        self.theta1 = theta1 * np.pi/180.0
        self.theta2 = theta2 * np.pi/180.0
        self.phi1 = phi1 * np.pi/180.0
        self.phi2 = phi2 * np.pi/180.0
        self.thetaDir1 = thetaDir1 * np.pi/180.0
        self.thetaDir2 = thetaDir2 * np.pi/180.0
        self.phiDir1 = phiDir1 * np.pi/180.0
        self.phiDir2 = phiDir2 * np.pi/180.0
        self.minDist = minDist
        self.pLength = pLength
        self.numSpheresPerParticle = numSpheresPerParticle
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockDirector(self):
        director = np.array([0.0, 0.0, 1.0])
        return director

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])
    
    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numSpheresPerParticle * self.numPoints 
        
    def pickFirstPoints(self):
        self.nList = [self.pickRandomPointInDefinedSpace()]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        theta, phi = coords.pickRandomPointOnUnitSphereInAngRange(self.thetaDir1, self.thetaDir2, self.phiDir1, self.phiDir2)
        dirn = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi]))
        return (coords.pickRandomPointInEllipsoidRange(self.xRadius, self.yRadius, self.zRadius, self.theta1, self.theta2, self.phi1, self.phi2), dirn) 
       
    def getParams(self):
        return self.params 

    def checkPointAgainstList(self, pos):
        ''' Function returns true if the test position does not intersect the other points.
            Assumes that pos is not in the list already. '''
        # assume position is good
        goodPos = True

        for zPos in self.nList:
            # define two spherocylinders by the end points of the center of the spheres at either end
            p1 = pos[0] + self.pLength/2.0 * pos[1]
            q1 = pos[0] - self.pLength/2.0 * pos[1]
            p2 = zPos[0] + self.pLength/2.0 * zPos[1]
            q2 = zPos[0] - self.pLength/2.0 * zPos[1]
            if cart.closestApproachTwoLineSegmentsSquared(p1, q1, p2, q2) < (self.minDist/2.0)**2:
                goodPos = False
                if self.verbose==1: 
                    print("Packing violation")
                break

        return goodPos
        
    def getNList(self):
        ''' Returns the longest nList found. Converts the sphero cylinders into a double length list of xyz points.
        Each pair of xyz points defines the endpoints of a spherocylinder.'''
        nList = []
        for pos in self.longestNList:
            for p in [ pos[0] + n * self.pLength/(self.numSpheresPerParticle - 1) * pos[1] for n in range(0, self.numSpheresPerParticle)]:
                nList.append(p)
        return nList 
            
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    EllipsoidPackSCBBG = VolumePackEllipsoidSCParticlesBBG(filename)

    numPoints = 100
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 0 
    xRadius = 5
    yRadius = 5
    zRadius = 5
    theta1 = -90.0
    theta2 = 90.0
    phi1 = -180.0
    phi2 = 180.0
    thetad1 = -90.0
    thetad2 = 90.0
    phid1 = -180.0
    phid2 = 180.0
    minDist = 2.0
    pLength = 2.0
    numSpheresPerParticle = 30
    
    envelopeList = ['None']
       
    # generate the building block
    EllipsoidPackSCBB = EllipsoidPackSCBBG.generateBuildingBlock(numPoints, xRadius, yRadius, zRadius, theta1, theta2, phi1, phi2, thetad1, thetad2, phid1, phid2, minDist, pLength, numSpheresPerParticle, envelopeList = envelopeList)
    EllipsoidPackSCBB.transformBBToLabFrame(director, centrePos, rotation)
    EllipsoidPackSCBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    