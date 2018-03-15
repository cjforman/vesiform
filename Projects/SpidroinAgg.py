import sys
import numpy as np
import random as rnd
import Utilities.fileIO as fIO
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Projects.Spidroin import spidroinProteinGenerator as SPG
from Library.SurfacePackEllipsoid import SurfacePackEllipsoidBBG as SPEBBG

class spidroinAggregateGenerator(BBG):
    # Takes a spidroin species 1 and 2 file and creates intermediate bundles.

    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        BBG.initialiseParameters(self) 
        
        self.spidroinSpecies1 = self.getParam('spidroinSpecies1.xyz')
        self.spidroinSpecies2 = self.getParam('spidroinSpecies2.xyz')
        
        self.ellipsoidPack = SPEBBG(self.paramFilename)
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Spidroin Object"
            sys.exit()        

    def generateBuildingBlock(self, numSpidroinSpecies1, numSpidroinSpecies2, numClustersInAggregate):
        
        self.numSpidroinsSpecies1 = numSpidroinSpecies1
        self.numSpidroinsSpecies2 = numSpidroinSpecies2
        self.numPointsInCluster = numSpidroinSpecies1 + numSpidroinSpecies2
        self.numClustersInAggregate = numClustersInAggregate
           
        return BBG.generateBuildingBlock(self, self.numPointsInCluster)
    
    def generateBuildingBlockXYZ(self):

        print "Generating Positions"
        spidroinPositionBB = SPEBBG.generateBuildingBlock(self, 
                                                          self.numPointsInCluster, 
                                                          self.clusterRX, 
                                                          self.clusterRY, 
                                                          self.clusterRZ, 
                                                          -90, 90, -180, 180, 
                                                          self.SpidroinRadius)
        
        spidroinPositions = spidroinPositionBB.xyzVals
        
        # compute individual spidroin generators
        spidroinGenVecs = [ spidPos for spidPos in spidroinPositions]
        spidroinGenVecsHat = [ spidGen/np.linalg.norm(spidGen) for spidGen in spidroinGenVecs ]
        
        # randomly orient the spidroins in the surface about their generators
        spidroinRots = [ rnd.uniform(0, 360) for n in range(0, self.numSpidroins)]
        
        print "Generating BuildingBlocks with given values"
        # compute the spidroin buildingblocks
        self.spidroins = [ self.SPG.generateBuildingBlock(spidPos, spidDir, spidRot, nameByBuildingBlockType=self.nameByBuildingBlockType) for spidPos, spidDir, spidRot in zip(spidroinPositions, spidroinGenVecsHat, spidroinRots) ]
        
        
        buildingBlockXYZ =[]
        for spid in self.spidroins:
            buildingBlockXYZ += spid.xyzVals
        
        return buildingBlockXYZ 
        
    def generateBuildingBlockNames(self):
        names = []
        for spid in self.spidroins:
            names += spid.atomNames
        return names

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
        
        
    