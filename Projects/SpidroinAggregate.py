import sys
import numpy as np
import random as rnd
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.SurfacePackEllipsoid import SurfacePackEllipsoidBBG as SPEBBG
from Library.peptideBackbonePDB import pdbPeptideBackboneGenerator as PDBGen


class spidroinAggregateGenerator(BBG):
    # Takes a spidroin species 1 and 2 file and creates intermediate bundles.

    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        BBG.initialiseParameters(self) 
        
        # N terminal PDB information and directors
        self.NTerminalFilename = self.getParam('NTerminalFilename')
        NTerminalDirector = self.getParam('NTerminalDirector')
        NTerminalDirector = np.array([component for component in NTerminalDirector[0]])
        self.NTermDirectorHat = NTerminalDirector/np.linalg.norm(NTerminalDirector) 
        self.NTermRot = self.getParam('NTermRot') * np.pi/180
        
        # C terminal PDB information and directors
        self.CTerminalFilename= self.getParam('CTerminalFilename')
        CTerminalDirector = self.getParam('CTerminalDirector')
        CTerminalDirector = np.array([component for component in CTerminalDirector[0]])
        self.CTermDirectorHat = CTerminalDirector/np.linalg.norm(CTerminalDirector)
        self.CTermRot = self.getParam('CTermRot') * np.pi/180                 
                 
        # Over all parameters used to describe the shape and size of the spidroin packing envelope.
        self.spidroinSpecies1FName = self.getParam('spidroinSpecies1Filename')
        self.spidroinSpecies2FName = self.getParam('spidroinSpecies2Filename')
        self.spidroinTerminusSeparation = self.getParam('spidroinTerminusSeparation')
        self.SpidroinRadius = self.getParam('SpidroinRadius')
        self.species2AltitudeBoost  = self.getParam('species2AltitudeBoost')
        self.TerminalExtension = self.getParam('TerminalExtension')
        self.clusterRX = self.getParam('clusterRX')
        self.clusterRY = self.getParam('clusterRY')
        self.clusterRZ = self.getParam('clusterRZ')
        self.aggregateRX = self.getParam('aggregateRX')
        self.aggregateRY = self.getParam('aggregateRY')
        self.aggregateRZ = self.getParam('aggregateRZ')
        self.minDist = self.getParam('minDist')
    
        self.SPEBBG = SPEBBG(self.paramFilename)
        self.CTerminusGen = PDBGen(self.CTerminalFilename)
        self.NTerminusGen = PDBGen(self.NTerminalFilename)
        
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

        print "Computing standard CTerminus for export as PDB and as backbone xyz"
        refPosCTerm = np.array([-self.spidroinTerminusSeparation/2.0, 0.0, -2.7]) # this last 2.7 is a total hack to get pointB in the envelope
        spidroinDirector = np.array([0.0, 0.0, 1.0])
        self.CTerminusAll = self.CTerminusGen.generateBuildingBlock(backboneOnly = False, director = self.CTermDirectorHat, showBlockDirector=False)
        self.CTerminusBackbone = self.CTerminusGen.generateBuildingBlock(backboneOnly = True, director = self.CTermDirectorHat, showBlockDirector=False)
        self.CTerminusAll.transformBBToLabFrame(spidroinDirector, refPosCTerm, self.CTermRot)
        self.CTerminusBackbone.transformBBToLabFrame(spidroinDirector, refPosCTerm, self.CTermRot) 
        
        print "Computing standard NTerminus for export as PDB and as backbone xyz"
        refPosNTerm  = np.array([self.spidroinTerminusSeparation/2.0, 0.0, 0.0])
        self.NTerminusAll = self.NTerminusGen.generateBuildingBlock(backboneOnly = False, director = self.NTermDirectorHat, showBlockDirector=False)
        self.NTerminusBackbone = self.NTerminusGen.generateBuildingBlock(backboneOnly = True, director = self.NTermDirectorHat, showBlockDirector=False)
        self.NTerminusAll.transformBBToLabFrame(spidroinDirector, refPosNTerm, self.NTermRot)
        self.NTerminusBackbone.transformBBToLabFrame(spidroinDirector, refPosNTerm, self.NTermRot) 
        
        print "Generating spidroin sites within cluster"
        spidroinPositionBB = self.SPEBBG.generateBuildingBlock(self.numPointsInCluster, 
                                                               self.clusterRX, 
                                                               self.clusterRY, 
                                                               self.clusterRZ, 
                                                               -90, 90, -180, 180, 
                                                               self.SpidroinRadius)
        
        # compute the positions and orientations of each individual spidroin within it's cluster
        spidroinPositions = spidroinPositionBB.blockXYZVals
        spidroinDirectors = [ spidPos for spidPos in spidroinPositions]
        spidroinDirectorsHat = [ spidDir/np.linalg.norm(spidDir) for spidDir in spidroinDirectors ]
        spidroinRots = [ rnd.uniform(0, 360) for _ in range(0, self.numPointsInCluster)]


        print "Generating cluster positions with in the aggregate"
        clusterPositionBB = self.SPEBBG.generateBuildingBlock(self.numClustersInAggregate, 
                                                              self.aggregateRX, 
                                                              self.aggregateRY, 
                                                              self.aggregateRZ, 
                                                              -90, 90, -180, 180, 
                                                              max(self.clusterRX, self.clusterRY, self.clusterRZ) + self.species2AltitudeBoost + self.TerminalExtension)
        clusterPositions = clusterPositionBB.blockXYZVals
        clusterDirectors = [ clustPos for clustPos in clusterPositions]
        clusterDirectorsHat = [ clustDir/np.linalg.norm(clustDir) for clustDir in clusterDirectors ]
        clusterRots = [ rnd.uniform(0, 360) for _ in range(0, self.numClustersInAggregate)]

        print "load spidroin coil information"
        spidroinCoilSpecies1Names, spidroinCoilSpecies1XYZ = fIO.loadXYZ(self.spidroinSpecies1FName) 
        spidroinCoilSpecies2Names, spidroinCoilSpecies2XYZ = fIO.loadXYZ(self.spidroinSpecies2FName) 
        
        print "Producing Cluster of clusters"
        # sety up output arrays
        NTermPoints = []
        CTermPoints = []
        xyzValCoils = []
        names = []
        clusterNum = 0
        for clustDirector, clustPosition, clustRotation in zip(clusterDirectorsHat, clusterPositions, clusterRots): 
            
            spidroinNum = 0
            # loop through the right number of times for the species1 data
            for rotation, director, position in zip(spidroinRots[0:self.numSpidroinSpecies1], spidroinDirectorsHat[0:self.numSpidroinSpecies1], spidroinPositions[0:self.numSpidroinSpecies1]):
                print clusterNum, spidroinNum
                if len(NTermPoints)==0:
                    xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies1XYZ[:])
                    # xyzVals = coords.transformFromBlockFrameToLabFrame(clustDirector, clustPosition, clustRotation, director, position,  xyzVals)
                    NTermPoints = xyzVals[0]
                    NTermPoints = xyzVals[-1]
                    xyzValCoils = xyzVals[:]
                    names = spidroinCoilSpecies1Names
                else:
                    xyzVals = coords.transformFromBlockFrameToLabFrame(director, position, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies1XYZ[:])
                    # xyzVals = coords.transformFromBlockFrameToLabFrame(clustDirector, clustPosition, clustRotation, director, position,  xyzVals)
                    NTermPoints = np.concatenate( (NTermPoints, xyzVals[0]), 0)
                    CTermPoints = np.concatenate( (CTermPoints, xyzVals[0]), 0)
                    xyzValCoils = np.concatenate( (xyzValCoils, xyzVals[:]), 0)                
                    names = np.concatenate( (names[:], spidroinCoilSpecies1Names), 0)
                spidroinNum += 1
            
            # loop throught the species2 data and add the species two data at a slightly higher altitude.
            for rotation, director, position in zip(spidroinRots[self.numSpidroinSpecies1:], spidroinDirectorsHat[self.numSpidroinSpecies1:], spidroinPositions[self.numSpidroinSpecies1:]):
                print clusterNum, spidroinNum
                xyzVals = coords.transformFromBlockFrameToLabFrame(director, position + director * self.species2AltitudeBoost, rotation, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), spidroinCoilSpecies2XYZ[:])
                # xyzVals = coords.transformFromBlockFrameToLabFrame(clustDirector, clustPosition, clustRotation, director, position,  xyzVals)
                NTermPoints = np.concatenate( (NTermPoints, xyzVals[0]), 0)
                CTermPoints = np.concatenate( (CTermPoints, xyzVals[0]), 0)
                xyzValCoils = np.concatenate( (xyzValCoils, xyzVals[:]), 0)                
                names = np.concatenate( (names[:], spidroinCoilSpecies2Names), 0)
                spidroinNum += 1
            clusterNum += 1
            self.namesTemp = names
        
        return xyzValCoils 
        
    def generateBuildingBlockNames(self):
        return self.namesTemp

    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

if __name__=="__main__":
    noErrors = True
    
    # generating spherical spidroin aggregate
    spidroinAggregateGenerator = spidroinAggregateGenerator('unsheared.txt')
    
    numSpidsSpecies1 = 15
    numSpidsSpecies2 = 10
    numClustersInAggregate = 1
    
    SpidroinCluster = spidroinAggregateGenerator.generateBuildingBlock(numSpidsSpecies1, numSpidsSpecies2, numClustersInAggregate)
    SpidroinCluster.exportBBK("species1Cluster")
    print("Spidroin Cluster Done.")
        
        
    