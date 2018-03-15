import sys
import numpy as np
import random as rnd
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.SurfacePackEllipsoid import SurfacePackEllipsoidBBG as SPEBBG

class spidroinAggregateGenerator(BBG):
    # Takes a spidroin species 1 and 2 file and creates intermediate bundles.

    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        BBG.initialiseParameters(self) 
        
        self.spidroinSpecies1FName = self.getParam('spidroinSpecies1.xyz')
        self.spidroinSpecies2FName = self.getParam('spidroinSpecies2.xyz')
        self.ellipsoidPack = SPEBBG(self.paramFilename)
        self.species2AltitudeBoost = self.getParam('species2AltitudeBoost')
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Spidroin Object"
            sys.exit()        

    def generateBuildingBlock(self, self, numSpidroinSpecies1, numSpidroinSpecies2, numClustersInAggregate):
        self.numSpidroinSpecies1 = numSpidroinSpecies1
        self.numSpidroinSpecies2 = numSpidroinSpecies2
        self.numPointsInCluster = numSpidroinSpecies1 + numSpidroinSpecies2
        self.numClustersInAggregate = numClustersInAggregate
        return BBG.generateBuildingBlock(self, self.numPointsInCluster, self.minDist ) 

    def generateBuildingBlockXYZ(self):

        print "Computing standard CTerminus for export as PDB"
        refPosCTerm = np.array([-self.spidroinTerminusSeparation/2.0, 0.0, -2.7]) # this last 2.7 is a total hack to get pointB in the envelope
        self.CTerminusAll = self.CTerminusGen.generateBuildingBlock(backboneOnly = False, director = self.CTermDirectorHat, showBlockDirector=False)
        self.CTerminusBackbone = self.CTerminusGen.generateBuildingBlock(backboneOnly = True, director = self.CTermDirectorHat, showBlockDirector=False)
        self.CTerminusAll.transformBBToLabFrame(self.spidroinDirector, refPosCTerm, self.CTermRot)
        self.CTerminusBackbone.transformBBToLabFrame(self.spidroinDirector, refPosCTerm, self.CTermRot) 
        
        print "Computing standard NTerminus for export as PDB"
        refPosNTerm  = np.array([self.spidroinTerminusSeparation/2.0, 0.0, 0.0])
        self.NTerminusAll = self.NTerminusGen.generateBuildingBlock(backboneOnly = False, director = self.NTermDirectorHat, showBlockDirector=False)
        self.NTerminusBackbone = self.NTerminusGen.generateBuildingBlock(backboneOnly = True, director = self.NTermDirectorHat, showBlockDirector=False)
        self.NTerminusAll.transformBBToLabFrame(self.spidroinDirector, refPosNTerm, self.NTermRot)
        self.NTerminusBackbone.transformBBToLabFrame(self.spidroinDirector, refPosNTerm, self.NTermRot) 
        
        
        print "Generating spidroin sites within cluster"
        spidroinPositionBB = SPEBBG.generateBuildingBlock(self, 
                                                          self.numPointsInCluster, 
                                                          self.clusterRX, 
                                                          self.clusterRY, 
                                                          self.clusterRZ, 
                                                          -90, 90, -180, 180, 
                                                          self.SpidroinRadius)
        # compute the positions and orientations of each individual spidroin within it's cluster
        spidroinPositions = spidroinPositionBB.xyzVals
        spidroinDirectors = [ spidPos for spidPos in spidroinPositions]
        spidroinDirectorsHat = [ spidDir/np.linalg.norm(spidDir) for spidDir in spidroinDirectors ]
        spidroinRots = [ rnd.uniform(0, 360) for _ in range(0, self.numSpidroins)]


        print "Generating cluster positions with in the aggregate"
        clusterPositionBB = SPEBBG.generateBuildingBlock(self, 
                                                         self.numClustersInAggregate, 
                                                         self.aggregateRX, 
                                                         self.aggregateRY, 
                                                         self.aggregateRZ, 
                                                         -90, 90, -180, 180, 
                                                         max(self.clusterRX, self.clusterRY, self.clusterRZ) + self.species2AltitudeBoost + self.TerminalExtension)
        clusterPositions = clusterPositionBB.xyzVals
        clusterDirectors = [ clustPos for clustPos in clusterPositions]
        clusterDirectorsHat = [ clustDir/np.linalg.norm(clustDir) for clustDir in clusterDirectors ]
        clusterRots = [ rnd.uniform(0, 360) for _ in range(0, self.numClustersInAggregate)]

        print "load spidroin coil information"
        spidroinCoilSpecies1XYZ, spidroinCoilSpecies1Names = fIO.loadXYZ(self.spidroinSpecies1FName) 
        spidroinCoilSpecies2XYZ, spidroinCoilSpecies2Names = fIO.loadXYZ(self.spidroinSpecies2FName) 
        
        print "Producing Cluster of clusters"
        # sety up output arrays
        NTermPoints = []
        CTermPoints = []
        xyzValCoils = []
        names = []
        
        for clustDirector, clustPosition, clustRotation in zip(clusterDirectorsHat, clusterPositions, clusterRots): 
        
            # loop through the right number of times for the species1 data
            for rotation, director, position in zip(spidroinRots[0:self.numSpidroinSpecies1], spidroinDirectorsHat[0:self.numSpidroinSpecies1], spidroinPositions[0:self.numSpidroinSpecies1]):
                if len(NTermPoints)==0:
                    xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies1XYZ[:])
                    xyzVals = coords.transformFromBlockFrameToLabFrame(clustDirector, clustPosition, clustRotation, director, position,  xyzVals)
                    NTermPoints = xyzVals[0]
                    NTermPoints = xyzVals[-1]
                    xyzValCoils = xyzVals
                    names = spidroinCoilSpecies1Names
                else:
                    xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies1XYZ[:])
                    xyzVals = coords.transformFromBlockFrameToLabFrame(clustDirector, clustPosition, clustRotation, director, position,  xyzVals)
                    NTermPoints = np.concatenate( (NTermPoints, xyzVals[0]), 0)
                    CTermPoints = np.concatenate( (CTermPoints, xyzVals[0]), 0)
                    xyzValCoils = np.concatenate( (xyzValCoils, xyzVals), 0)                
                    names = np.concatenate( (names, spidroinCoilSpecies1Names), 0)
            
            # loop throught the species2 data and add the species two data at a slightly higher altitude.
            for rotation, director, position in zip(spidroinRots[self.numSpidroinSpecies1:], spidroinDirectorsHat[self.numSpidroinSpecies1:], spidroinPositions[self.numSpidroinSpecies1:]):
                    xyzVals = coords.transformFromBlockFrameToLabFrame(director, position + director * self.species2AltitudeBoost, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies2XYZ[:])
                    xyzVals = coords.transformFromBlockFrameToLabFrame(clustDirector, clustPosition, clustRotation, director, position,  xyzVals)
                    NTermPoints = np.concatenate( (NTermPoints, xyzVals[0]), 0)
                    CTermPoints = np.concatenate( (CTermPoints, xyzVals[0]), 0)
                    xyzValCoils = np.concatenate( (xyzValCoils, xyzVals), 0)                
                    names = np.concatenate( (names, spidroinCoilSpecies2Names), 0)
        self.namesTemp = names
        
        return xyzValCoils 
        
    def generateBuildingBlockNames(self):
        return self.namesTemp

    def generateBuildingBlockDirector(self):
        return self.spidroinAggDirector

    def generateBuildingBlockRefPoint(self):
        return self.spidroinAggCentrePos

if __name__=="__main__":
    noErrors = True
    
    sphericalFilename = sys.argv[1]
    ellipticalFilename = sys.argv[2]
    
    # generating spherical spidroin aggregate
    spidroinAggregateGeneratorSpherical = spidroinAggregateGenerator(sphericalFilename)
    
    startPoint = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    rotation = rnd.uniform(-180, 180)
    numSpids = 30
    rx = 120
    ry = 120
    rz = 120
    theta1 = -90
    theta2 = 90
    phi1 = -180
    phi2 = 180
    
    spidAggBBSpherical = spidroinAggregateGeneratorSpherical.generateBuildingBlock(numSpids,
                                                                 startPoint,
                                                                 direction, 
                                                                 rotation,
                                                                 rx,
                                                                 ry,
                                                                 rz,
                                                                 theta1,
                                                                 theta2,
                                                                 phi1,
                                                                 phi2,
                                                                 alignDirectors=True,
                                                                 showDirector=False,
                                                                 nameByBuildingBlockType=True)
    
    spidAggBBSpherical.exportBBK(fIO.fileRootFromInfile(sphericalFilename, 'txt'))


    # generating elliptical structures 
    spidroinAggregateGeneratorElliptical = spidroinAggregateGenerator(ellipticalFilename)
    
    # changed parameters
    rx = 135
    ry = 180
    rz = 70

    spidAggBBElliptical = spidroinAggregateGeneratorElliptical.generateBuildingBlock(numSpids,
                                                                 startPoint,
                                                                 direction, 
                                                                 rotation,
                                                                 rx,
                                                                 ry,
                                                                 rz,
                                                                 theta1,
                                                                 theta2,
                                                                 phi1,
                                                                 phi2,
                                                                 alignDirectors=True,
                                                                 showDirector=False,
                                                                 nameByBuildingBlockType=True)
    
    spidAggBBElliptical.exportBBK(fIO.fileRootFromInfile(ellipticalFilename, 'txt'))


    print("Done.")
        
        
    