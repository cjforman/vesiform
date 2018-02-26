import sys
import numpy as np
import random as rnd
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.SpacePackEllipsoid import SpacePackEllipsoidBBG as SPEBBG
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.peptideHairpin import peptideHairpinGenerator as PHG
from Library.Ellipsoid import EllipsoidPackNBB as EPNBB
from Projects.betasheet import betasheetGen as BSG
from NPack.NPackBB import NPackBB as NBB

class spidroinProteinGenerator(BBG):
    # Each spidroin protein consists of two alpha helical termini and a main 
    # body that consists of regions of beta sheet structures linked 
    # by hairpin turns.

    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        BBG.initialiseParameters(self) 
        
        self.numResiduesAlphaHelix = self.getParam('numResiduesAlphaHelix')
        self.lengthPerResidueAlpha = self.getParam('lengthPerResidueAlpha')
        self.lengthPerResidueBeta = self.getParam('lengthPerResidueBeta')
        self.terminiSeparation = self.getParam('terminiSeparation')
        self.numBetaSheets = self.getParam('numBetaSheets')
        self.numBetaStrands = self.getParam('numBetaStrands')
        self.betaStrandLength = self.getParam('betaStrandLength')
        self.betaSheetRx = self.getParam('betaSheetRx')
        self.betaSheetRy = self.getParam('betaSheetRy')
        self.betaSheetRz = self.getParam('betaSheetRz')
        self.alphaBetaSeparation = self.getParam('alphaBetaSeparation')
        self.alphaHelixAtomName = self.getParam('alphaHelixAtomName') 
        self.betaStrandAtomName = self.getParam('betaStrandAtomName')
        self.coilAtomName = self.getParam('coilAtomName')
        
        
        self.BSG = BSG(self.paramFilename)
        self.AHG = PBG(self.paramFilename)
        self.PHG = PHG(self.paramFilename)
        self.USP = NBB(self.paramFilename)
        self.SPEBBG = SPEBBG(self.paramFilename)
        
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

        # build the components
        
        print "Constructing N termninus alpha helix"
        AH_NTerm_numPoints = self.numResiduesAlphaHelix * 3
        AH_NTerm_startPos = np.array([-self.terminiSeparation/2, 
                                      0.0, 
                                      self.numResiduesAlphaHelix * self.lengthPerResidueAlpha])
        AH_NTerm_director = 1 * self.spidroinDirector
        AH_NTerm_rotation = 0
        AH_NTerm_polarity = 'NC'
         
        # Create the N terminus alpha helix bundle
        self.NTerminusAlphaHelix = self.AHG.generateBuildingBlock( AH_NTerm_numPoints,
                                                                   AH_NTerm_startPos,
                                                                   AH_NTerm_director,
                                                                   AH_NTerm_rotation,
                                                                   AH_NTerm_polarity,
                                                                   alignDirectors=True,
                                                                   showDirector=False)
        if self.nameByBuildingBlockType:
            self.NTerminusAlphaHelix.replaceNames(self.alphaHelixAtomName)
        
        print "Constructing C termninus alpha helix"
        AH_CTerm_numPoints = self.numResiduesAlphaHelix * 3
        AH_CTerm_startPos = np.array([self.terminiSeparation/2, 
                                      0.0, 
                                      0.0])
        AH_CTerm_director = -1 * self.spidroinDirector
        AH_CTerm_rotation = 0
        AH_CTerm_polarity = 'NC'

        # Create the C terminus alpha helix bundle
        self.CTerminusAlphaHelix = self.AHG.generateBuildingBlock( AH_CTerm_numPoints,
                                                                   AH_CTerm_startPos,
                                                                   AH_CTerm_director,
                                                                   AH_CTerm_rotation,
                                                                   AH_CTerm_polarity,
                                                                   alignDirectors=True,
                                                                   showDirector=False)
        
        if self.nameByBuildingBlockType:
            self.CTerminusAlphaHelix.replaceNames(self.alphaHelixAtomName) 
        
        BS_centre = np.array([0.0, 0.0, -self.betaSheetRz - self.alphaBetaSeparation])
        BS_director = self.spidroinDirector
        BS_rotation = 0
        BS_minDist = 0.8 * np.sqrt(2) * max([4.8 * self.numBetaStrands, self.betaStrandLength * self.lengthPerResidueBeta]) 
        
        print "Constructing Beta Sheet Start Points"
        
        # then create the start point for a bundle of beta sheets packed uniformly in space
        betaSheetStartPointsBB = self.SPEBBG.generateBuildingBlock( self.numBetaSheets, 
                                                                    BS_centre, 
                                                                    BS_director,
                                                                    BS_rotation, 
                                                                    self.betaSheetRx, 
                                                                    self.betaSheetRy, 
                                                                    self.betaSheetRz, 
                                                                    BS_minDist)
        betaSheetStartPoints = betaSheetStartPointsBB.xyzVals 

        betaSheetDirectors =[]
        # create beta sheet directors
        for n in range(0, self.numBetaSheets):
            theta, phi = coords.pickRandomPointOnUnitSphere()
            betaSheetDirectors.append(coords.sphericalPolar2XYZ(np.array([1.0, theta, phi])))

        inStrandDirector = np.array([0.0, 0.0, 1.0])
        crossStrandDirector = np.array([1.0, 0.0, 0.0])
        rotation = rnd.uniform(0, 2*np.pi)
        offset = np.array([0.0, 0.0, 0.0])                             
        
        print "Constructing Beta Sheets"

        # create the betaSheets
        self.betaSheetBBs = [ self.BSG.generateBuildingBlock( self.numBetaStrands, 
                                                            self.betaStrandLength, 
                                                            startPos, 
                                                            globalDirector,
                                                            inStrandDirector, 
                                                            crossStrandDirector, 
                                                            rotation, 
                                                            offset, 
                                                            polarity='NC', 
                                                            parallel=True, 
                                                            loopedEnds=False) for startPos, globalDirector in zip(betaSheetStartPoints, betaSheetDirectors)]

        if self.nameByBuildingBlockType:
            [ bsheet.replaceNames(self.betaStrandAtomName) for bsheet in self.betaSheetBBs ]

        # assemble the components
        spidroinXYZ = self.NTerminusAlphaHelix.xyzVals
        for betaSheet in self.betaSheetBBs:
            spidroinXYZ += betaSheet.xyzVals
        spidroinXYZ += self.CTerminusAlphaHelix.xyzVals
        
        return spidroinXYZ
    
    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockNames(self):
        spidroinNames = self.NTerminusAlphaHelix.atomNames
        for betaSheet in self.betaSheetBBs:
            spidroinNames += betaSheet.atomNames
        spidroinNames += self.CTerminusAlphaHelix.atomNames
        return spidroinNames

if __name__=="__main__":
    noErrors = True
    
    filename = sys.argv[1]
    
    # generating an individual spidroin
    spidroinProteinGenerator = spidroinProteinGenerator(filename)
    startPoint = np.array([0.0, -0.0, 0.0])
    direction = np.array([-0.0, 0.0, 1.0])
    rotation = 0
    
    SpidroinBB = spidroinProteinGenerator.generateBuildingBlock(startPoint, direction, rotation, alignDirectors=True, showDirector=False, nameByBuildingBlockType=True)  
    SpidroinBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print("Done.")
        
        
    