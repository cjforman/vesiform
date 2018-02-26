import sys
import numpy as np
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
import copy as cp
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.peptideBackbonePDB import pdbPeptideBackboneGenerator as PDBGen


class spidroinTerminalGenerator(BBG):
    # takes a spidroin dimer and reverse engineers the connector of the two dimers.
    # takes two monomers and reconstructs the dimer exporting as an xyz
    # Does this for N and C terminus proteins.

    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        BBG.initialiseParameters(self) 
        
        # get the C Terminal information about dimers and monomers
        self.CTermDiFilename = self.getParam('CTerminalFilenameDimer')
        CTermDiCI = self.getParam('CTermDimerConnectorIndices')
        self.CTermDiConnectorIndices = [ [ CTermDiCI[0], CTermDiCI[1], CTermDiCI[2] ],
                                         [ CTermDiCI[3], CTermDiCI[4], CTermDiCI[5] ] ]
        
        self.CTermMonFilename = self.getParam('CTerminalFilenameMonomer')
        CTermMonCI = self.getParam('CTermMonomerConnectorIndices')
        self.CTermMonConnectorIndices = [ [ CTermMonCI[0], CTermMonCI[1], CTermMonCI[2] ] ]

        # get the C Terminal information about dimers and monomers
        self.NTermDiFilename = self.getParam('NTerminalFilenameDimer')
        NTermDiCI = self.getParam('NTermDimerConnectorIndices')
        self.NTermDiConnectorIndices = [ [ NTermDiCI[0], NTermDiCI[1], NTermDiCI[2] ],
                                         [ NTermDiCI[3], NTermDiCI[4], NTermDiCI[5] ] ]
        
        # get the monomer information
        self.NTermMonFilename = self.getParam('NTerminalFilenameMonomer')
        NTermMonCI = self.getParam('NTermMonomerConnectorIndices')
        self.NTermMonConnectorIndices = [ [ NTermMonCI[0], NTermMonCI[1], NTermMonCI[2] ] ]
        
        self.NTermMonGen = PDBGen(self.NTermMonFilename)
        self.NTermDiGen = PDBGen(self.NTermDiFilename)
        self.CTermMonGen = PDBGen(self.CTermMonFilename)
        self.CTermDiGen = PDBGen(self.CTermDiFilename)
        
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Spidroin Object"
            sys.exit()        

    def generateBuildingBlock(self, startPos, direction, rotation, alignDirectors=True, showDirector=False, nameByBuildingBlockType=True):
        self.numPoints = 0
        self.spidroinDirector = np.array([0.0, 0.0, 1.0])
        self.spidroindRefPoint = np.array([0.0, 0.0, 0.0])
        self.nameByBuildingBlockType = nameByBuildingBlockType
        
        return BBG.generateBuildingBlock(self, self.numPoints, startPos, direction, rotation, alignDirectors=alignDirectors, showDirector=showDirector)
    
    def generateBuildingBlockXYZ(self):

        # Create the Dimer Building Blocks
        startIndex = 0
        startPos = np.array([0.0, 0.0, 0.0])
        direction = np.array([0.0, 0.0, 1.0])
        rotation = 0
        self.NTermDiBB = self.NTermDiGen.generateBuildingBlock(startIndex, self.NTermDiConnectorIndices, startPos, direction, rotation, alignDirectors=False, showDirector=False, numPoints='all' )  
        self.CTermDiBB = self.CTermDiGen.generateBuildingBlock(startIndex, self.CTermDiConnectorIndices, startPos, direction, rotation, alignDirectors=False, showDirector=False, numPoints='all' )
        
        self.CTermDiBB.replaceName(self.CTermDiConnectorIndices[1][0], 'H')  # m0
        self.CTermDiBB.replaceName(self.CTermDiConnectorIndices[1][1], 'O')  # m1
        self.CTermDiBB.replaceName(self.CTermDiConnectorIndices[1][2], 'S')  # m2
        self.CTermDiBB.replaceName(self.CTermDiConnectorIndices[0][0], 'H')  # s0
        self.CTermDiBB.replaceName(self.CTermDiConnectorIndices[0][1], 'O')  # s1
        self.CTermDiBB.replaceName(self.CTermDiConnectorIndices[0][2], 'S')  # s2
        
        self.NTermDiBB.replaceName(self.NTermDiConnectorIndices[1][0], 'H')  # m0
        self.NTermDiBB.replaceName(self.NTermDiConnectorIndices[1][1], 'O')  # m1
        self.NTermDiBB.replaceName(self.NTermDiConnectorIndices[1][2], 'S')  # m2
        self.NTermDiBB.replaceName(self.NTermDiConnectorIndices[0][0], 'H')  # s0
        self.NTermDiBB.replaceName(self.NTermDiConnectorIndices[0][1], 'O')  # s1
        self.NTermDiBB.replaceName(self.NTermDiConnectorIndices[0][2], 'S')  # s2
        
        
        self.CTermDiBB.exportBBK('CTermDimer_names')
        self.NTermDiBB.exportBBK('NTermDimer_names')

              
        # Reverse engineer connectors and parameters from the dimer models
        # the connector indices here have been manually selected from the dimer files for N and T terminals
        # and placed in the input file for this process 
        NS = self.NTermDiBB.getConnectionAtoms(0)
        NM = self.NTermDiBB.getConnectionAtoms(1)
        self.NTermConnectionInfo = coords.measureAnglesAtConnection(NS[0], NS[1], NS[2], NM[2], NM[1], NM[0])
        CS = self.CTermDiBB.getConnectionAtoms(0)
        CM = self.CTermDiBB.getConnectionAtoms(1)
        self.CTermConnectionInfo = coords.measureAnglesAtConnection(CS[0], CS[1], CS[2], CM[2], CM[1], CM[0])

        print "NTerm Angles: ", self.NTermConnectionInfo
        print "CTerm Angles: ", self.CTermConnectionInfo
        
        fIO.writeTextFile([ str(value) + '\n' for value in self.NTermConnectionInfo], "NTermAngles.txt")
        fIO.writeTextFile([ str(value) + '\n' for value in self.CTermConnectionInfo], "CTermAngles.txt")

        self.NTermMonBB = self.NTermMonGen.generateBuildingBlock(startIndex, [self.NTermMonConnectorIndices[0]], startPos, direction, rotation, alignDirectors=False, showDirector=False, numPoints='all' )  
        self.CTermMonBB = self.CTermMonGen.generateBuildingBlock(startIndex, [self.CTermMonConnectorIndices[0]], startPos, direction, rotation, alignDirectors=False, showDirector=False, numPoints='all' )

        self.NTermDiMonBB, self.NTermStableBB = self.NTermMonBB.addBuildingBlock(cp.copy(self.NTermMonBB),
                                                                                 0,
                                                                                 0,
                                                                                 self.NTermConnectionInfo[0],
                                                                                 self.NTermConnectionInfo[1],
                                                                                 self.NTermConnectionInfo[2],
                                                                                 self.NTermConnectionInfo[3],
                                                                                 self.NTermConnectionInfo[4],
                                                                                 self.NTermConnectionInfo[5],
                                                                                 startPos,
                                                                                 direction) 

        self.CTermDiMonBB, self.CTermStableBB = self.CTermMonBB.addBuildingBlock(cp.copy(self.CTermMonBB),
                                                                                 0,
                                                                                 0,
                                                                                 self.CTermConnectionInfo[0],
                                                                                 self.CTermConnectionInfo[1],
                                                                                 self.CTermConnectionInfo[2],
                                                                                 self.CTermConnectionInfo[3],
                                                                                 self.CTermConnectionInfo[4],
                                                                                 self.CTermConnectionInfo[5],
                                                                                 startPos + np.array([10.0, 0.0, 0.0]),
                                                                                 direction) 


        self.CTermDiMonBB.replaceName(self.CTermMonConnectorIndices[0][0] + self.CTermMonBB.countAtoms(), 'S')  # m0
        self.CTermDiMonBB.replaceName(self.CTermMonConnectorIndices[0][1] + self.CTermMonBB.countAtoms(), 'S')  # m1
        self.CTermDiMonBB.replaceName(self.CTermMonConnectorIndices[0][2] + self.CTermMonBB.countAtoms(), 'S')  # m2
        self.CTermDiMonBB.replaceName(self.CTermMonConnectorIndices[0][0], 'S')  # s0
        self.CTermDiMonBB.replaceName(self.CTermMonConnectorIndices[0][1], 'S')  # s1
        self.CTermDiMonBB.replaceName(self.CTermMonConnectorIndices[0][2], 'S')  # s2
        
        self.NTermDiMonBB.replaceName(self.NTermMonConnectorIndices[0][0] + self.NTermMonBB.countAtoms(), 'S')  # m0
        self.NTermDiMonBB.replaceName(self.NTermMonConnectorIndices[0][1] + self.NTermMonBB.countAtoms(), 'S')  # m1
        self.NTermDiMonBB.replaceName(self.NTermMonConnectorIndices[0][2] + self.NTermMonBB.countAtoms(), 'S')  # m2
        self.NTermDiMonBB.replaceName(self.NTermMonConnectorIndices[0][0], 'S')  # s0
        self.NTermDiMonBB.replaceName(self.NTermMonConnectorIndices[0][1], 'S')  # s1
        self.NTermDiMonBB.replaceName(self.NTermMonConnectorIndices[0][2], 'S')  # s2

        self.NTermDiMonBB.exportBBK(self.NTermDiFilename[0:-4] + '_DiMon.xyz')
        self.CTermDiMonBB.exportBBK(self.CTermDiFilename[0:-4] + '_DiMon.xyz')


        
        # assemble the components
        spidroinXYZ = self.NTermMonBB.xyzVals
        spidroinXYZ += [ xyzVal + np.array([10.0, 0.0, 0.0]) for xyzVal in self.CTermMonBB.xyzVals ]
                                                             
        return spidroinXYZ
        
    
    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockNames(self):
        spidroinNames = self.NTermMonBB.atomNames
        spidroinNames += self.CTermMonBB.atomNames
        
        return spidroinNames

if __name__=="__main__":
    noErrors = True
    
    filename = sys.argv[1]
    
    # generating an individual spidroin
    spidroinTerminalGenerator = spidroinTerminalGenerator(filename)
    startPoint = np.array([0.0, -0.0, 0.0])
    direction = np.array([-0.0, 0.0, 1.0])
    rotation = 0
    
    SpidroinBB = spidroinTerminalGenerator.generateBuildingBlock(startPoint, direction, rotation, alignDirectors=False, showDirector=False, nameByBuildingBlockType=True)  
    SpidroinBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print("Done.")
        
        
    