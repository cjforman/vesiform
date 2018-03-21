import sys
import numpy as np
import random as rnd
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.SurfacePackEllipsoid import SurfacePackEllipsoidBBG as SPEBBG
from Library.VolumePackEllipsoidSCP import VolumePackEllipsoidSCParticlesBBG as VPESCPBBG


class spidroinAggregateGenerator(BBG):
    # Takes a spidroin species 1 and 2 file and creates intermediate bundles.

    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        BBG.initialiseParameters(self) 
                 
        # Over all parameters used to describe the shape and size of the spidroin packing envelope.
        self.spidroinSpecies1FName = self.getParam('spidroinSpecies1Filename')
        self.spidroinSpecies2FName = self.getParam('spidroinSpecies2Filename')
        self.terminalSeparation = self.getParam('terminalSeparation')
        self.monomerDiameter = self.getParam('monomerDiameter')
        self.clusterRX = self.getParam('clusterRX')
        self.clusterRY = self.getParam('clusterRY')
        self.clusterRZ = self.getParam('clusterRZ')
        self.aggregateRX = self.getParam('aggregateRX')
        self.aggregateRY = self.getParam('aggregateRY')
        self.aggregateRZ = self.getParam('aggregateRZ')
        self.minDist = self.getParam('minDist')
    
        self.SPEBBG = SPEBBG(self.paramFilename)
        self.VPESCPBBG = VPESCPBBG(self.paramFilename)
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Spidroin Object"
            sys.exit()        

    def generateBuildingBlock(self, numSpidroinSpecies1, numSpidroinSpecies2, numClustersInAggregate):
        self.numSpidroinSpecies1 = numSpidroinSpecies1
        self.numSpidroinSpecies2 = numSpidroinSpecies2
        self.numPointsInCluster = numSpidroinSpecies1 + numSpidroinSpecies2
        self.numClustersInAggregate = numClustersInAggregate
        return BBG.generateBuildingBlock(self, self.numPointsInCluster, self.minDist )
    
    def generateBuildingBlockXYZ(self):
        
        print "load spidroin coil information"
        spidroinCoilSpecies1Names, spidroinCoilSpecies1XYZ = fIO.loadXYZ(self.spidroinSpecies1FName) 
        spidroinCoilSpecies2Names, spidroinCoilSpecies2XYZ = fIO.loadXYZ(self.spidroinSpecies2FName) 
        
        print "Generating spidroin sites within cluster"
        spidroinPositionBB = self.SPEBBG.generateBuildingBlock(self.numPointsInCluster, 
                                                               self.clusterRX, 
                                                               self.clusterRY, 
                                                               self.clusterRZ, 
                                                               -90.0, 90.0, -180.0, 180.0, 
                                                               self.terminalSeparation)

        
        # compute the positions and orientations of each individual spidroin within it's cluster
        spidroinPositions = spidroinPositionBB.blockXYZVals
        spidDirAngles = [ rnd.uniform(0, 360) for _ in spidroinPositions]
        spidroinDirectors = [ coords.sphericalPolar2XYZ(np.array([1.0, 85.0 * np.pi/180.0, angle * np.pi/180.0])) for angle in spidDirAngles]
        spidroinDirectorsHat = [ director/np.linalg.norm(director) for director in spidroinDirectors]
        spidroinRots = [ -180.0 + 180.0 * rnd.randint(0, 1) for _ in range(0, self.numPointsInCluster)]

        print "Generating cluster positions within the larger aggregate"
        clusterPositionBB = self.VPESCPBBG.generateBuildingBlock( self.numClustersInAggregate, 
                                                                  self.aggregateRX, 
                                                                  self.aggregateRY, 
                                                                  self.aggregateRZ,
                                                                  -90.0, 90.0, -180.0, 180.0,
                                                                  self.monomerDiameter, 
                                                                  self.terminalSeparation, 
                                                                  2)
        clusterPositionInfo = clusterPositionBB.blockXYZVals
        clusterPositions = [ (pos + clusterPositionInfo[2*i + 1])/2.0 for i, pos in enumerate(clusterPositionInfo[0:-2:2])]
        clusterDirectors = [ (clusterPositionInfo[2*i + 1] - pos) for i, pos in enumerate(clusterPositionInfo[0:-2:2])]
        clusterDirectorsHat = [ clustDir/np.linalg.norm(clustDir) for clustDir in clusterDirectors ]
        clusterRots = [ 0.0 for _ in range(0, self.numClustersInAggregate)]

        fIO.saveXYZList(spidroinPositions, ['Ca', 'O'] * self.numPointsInCluster, 'aggClusterPoints.xyz')

        print "Producing Cluster"
        # sety up output arrays
        clusterXYZ = []
        clusterNames = []
        spidroinNum = 0
        # loop through the right number of times for the species1 data
        for rotation, director, position in zip(spidroinRots[0:self.numSpidroinSpecies1], spidroinDirectorsHat[0:self.numSpidroinSpecies1], spidroinPositions[0:self.numSpidroinSpecies1]):
            print spidroinNum
            if len(clusterXYZ)==0:
                xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies1XYZ[:])
                clusterXYZ = xyzVals[:]
                clusterNames = spidroinCoilSpecies1Names
            else:
                xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies1XYZ[:])
                clusterXYZ = np.concatenate( (clusterXYZ, xyzVals[:]), 0)                
                clusterNames = np.concatenate( (clusterNames[:], spidroinCoilSpecies1Names), 0)
            spidroinNum += 1
        
        # loop throught the species2 data and add the species two data at a slightly higher altitude.
        for rotation, director, position in zip(spidroinRots[self.numSpidroinSpecies1:], spidroinDirectorsHat[self.numSpidroinSpecies1:], spidroinPositions[self.numSpidroinSpecies1:]):
            print spidroinNum
            xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies2XYZ[:])
            clusterXYZ = np.concatenate( (clusterXYZ, xyzVals[:]), 0)                
            clusterNames = np.concatenate( (clusterNames[:], spidroinCoilSpecies2Names), 0)
            spidroinNum += 1
        
        fIO.saveXYZList(clusterXYZ, clusterNames, "cluster.xyz")
        
        print "Producing aggregate"
        # sety up output arrays
        aggregateXYZ = []
        aggregateNames = []
        clusterNum = 0
        # loop through the right number of times for the species1 data
        for rotation, director, position in zip(clusterRots, clusterDirectorsHat, clusterPositions):
            print clusterNum
            if len(aggregateXYZ)==0:
                xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), clusterXYZ)
                aggregateXYZ = xyzVals[:]
                aggregateNames = clusterNames
            else:
                xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), clusterXYZ)
                aggregateXYZ = np.concatenate( (aggregateXYZ, xyzVals[:]), 0)                
                aggregateNames = np.concatenate( (aggregateNames, clusterNames), 0)
            clusterNum += 1
        
        fIO.saveXYZList(aggregateXYZ, aggregateNames, "spidroinAggregate.xyz")
        
        self.namesTemp = aggregateNames
        return aggregateXYZ 
        
    def generateBuildingBlockNames(self):
        return self.namesTemp

    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

if __name__=="__main__":
    noErrors = True
    
    # generating spherical spidroin aggregate
    spidroinAggregateGenerator = spidroinAggregateGenerator('sheared.txt')
    
    numSpidsSpecies1 = 10
    numSpidsSpecies2 = 2
    numClustersInAggregate = 30
    
    SpidroinCluster = spidroinAggregateGenerator.generateBuildingBlock(numSpidsSpecies1, numSpidsSpecies2, numClustersInAggregate)
    print("Spidroin Sheared Cluster Done.")
        
        
    