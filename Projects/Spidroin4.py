import sys
import numpy as np
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.peptideBackbonePDB import pdbPeptideBackboneGenerator as PDBGen
from Projects.spidroinHairpin import spidroinHairpinGenerator as SHGen
from Projects.chainOfChains import chainOfChainsGenerator as CCGen

class spidroinProteinGenerator(BBG):
    # Each spidroin protein consists of two alpha helical termini and a main 
    # body that consists of a coarse grain region in a constrained loop

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
        self.spidroinTerminusSeparation = self.getParam('spidroinTerminusSeparation')
        self.SpidroinFrustumMaxRadius = self.getParam('SpidroinFrustumMaxRadius')
        self.SpidroinFrustumMinRadius = self.getParam('SpidroinFrustumMinRadius')
        self.SpidroinFrustumZ1 = 0.0
        self.SpidroinFrustumZ2 = self.getParam('SpidroinFrustumZ2')
        self.CNbondLength = self.getParam('CNbondLength')
        self.CCbondLength = self.getParam('CCbondLength')
        self.omega = self.getParam('omega')
        self.angleN = self.getParam('angleN')
        self.psi = self.getParam('psi')
        self.angleC = self.getParam('angleC')
        self.phi = self.getParam('phi')
        self.angleC = self.getParam('angleC')
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')       
        self.verbose = self.getParam('verbose')
        self.SP1NumGUnits = self.getParam('SP1NumGUnits')
        self.SP2NumGUnits = self.getParam('SP2NumGUnits') 
        self.envRadiusFactor = self.getParam('envRadiusFactor')
        self.SP1LenGunitRes = self.getParam('SP1LenGunitRes')
        self.SP1LenPQunitRes = self.getParam('SP1LenPQunitRes')
        self.SP2LenGunitRes = self.getParam('SP2LenGunitRes')
        self.SP2LenPQunitRes = self.getParam('SP2LenPQunitRes')
         
         
        # for colouring regions distinct colours in blender:
        self.NTerminalAtomName = self.getParam('NTerminalAtomName') 
        self.CTerminalAtomName = self.getParam('CTerminalAtomName')
        self.spidroinHairpinAtomName = self.getParam('spidroinHairpinAtomName')
        self.spidroinAtomName = self.getParam('spidroinAtomName')
        
        # create instances of all the building block generators that we will need. 
        self.SpidroinCoilGen = SHGen(self.paramFilename)
        self.CTerminusGen = PDBGen(self.CTerminalFilename)
        self.NTerminusGen = PDBGen(self.NTerminalFilename)
        self.CCGen = CCGen(self.paramFilename)
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Spidroin Object"
            sys.exit()        

    def generateBuildingBlock(self, species, minDist, showBlockDirector=False):
        self.species = species
        dummyNumPoints = 0
        self.minDist = minDist
        self.spidroinDirector = np.array([0.0, 0.0, 1.0])
        self.spidroinRefPoint = np.array([0.0, 0.0, 0.0])
        return BBG.generateBuildingBlock(self, dummyNumPoints, minDist, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockXYZ(self):
        # The spidroin model consists of two alpha helical termini proteins which are taken from PDB files
        # and a spidroinHairpin Coil which is a coarse grain representation of the spidroin coil that connects 
        # the two points at the end of the termini.  
        # The centre of mass of the two termini are placed at +/- spidroinTerminusSeparation/2.0 on the x axis.
        
        print "CTerminus"
        refPosCTerm = np.array([-self.spidroinTerminusSeparation/2.0, 0.0, -2.7]) # this last 2.7 is a total hack to get pointB in the envelope
        self.CTerminusAll = self.CTerminusGen.generateBuildingBlock(backboneOnly = False, director = self.CTermDirectorHat, showBlockDirector=False)
        self.CTerminusBackbone = self.CTerminusGen.generateBuildingBlock(backboneOnly = True, director = self.CTermDirectorHat, showBlockDirector=False)
        self.CTerminusAll.transformBBToLabFrame(self.spidroinDirector, refPosCTerm, self.CTermRot)
        self.CTerminusBackbone.transformBBToLabFrame(self.spidroinDirector, refPosCTerm, self.CTermRot) 
        self.CTerminusGen.exportPDBWithNewCoords(self.CTerminusAll.blockXYZVals, 'CTerminus_Positioned.pdb')
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZ(self.CTerminusBackbone.blockXYZVals, self.CTerminalAtomName, "CTerminal.xyz")
        
        print "NTerminus"
        refPosNTerm  = np.array([self.spidroinTerminusSeparation/2.0, 0.0, 0.0])
        self.NTerminusAll = self.NTerminusGen.generateBuildingBlock(backboneOnly = False, director = self.NTermDirectorHat, showBlockDirector=False)
        self.NTerminusBackbone = self.NTerminusGen.generateBuildingBlock(backboneOnly = True, director = self.NTermDirectorHat, showBlockDirector=False)
        self.NTerminusAll.transformBBToLabFrame(self.spidroinDirector, refPosNTerm, self.NTermRot)
        self.NTerminusBackbone.transformBBToLabFrame(self.spidroinDirector, refPosNTerm, self.NTermRot) 
        self.NTerminusGen.exportPDBWithNewCoords(self.NTerminusAll.blockXYZVals, 'NTerminus_Positioned.pdb')
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZ(self.NTerminusBackbone.blockXYZVals, self.NTerminalAtomName, "NTerminal.xyz")
        
        print "Compute end points of spidroin coil"
        # Compute the coil start (A) and end (B) points. 
        ConnectionA = self.NTerminusBackbone.getConnectionAtoms(1)

        # N (psi) C (phi) C (omega) N (psi) C (phi) C connection
        # s0      s1        s2      m2      m1        m0
        # a C to N terminus connection
        pointsA = self.findCoilPoints(ConnectionA[0], 
                                      ConnectionA[1], 
                                      ConnectionA[2], 
                                      self.CNbondLength,
                                      self.CNbondLength,
                                      self.CCbondLength,
                                      self.phi,
                                      self.angleC,
                                      self.omega,
                                      self.angleN,
                                      self.psi,
                                      self.angleC)
        pointA = pointsA[2]
        
        # find the N-terminus end of the C terminus and calculate where the coil should start
        ConnectionB = self.CTerminusBackbone.getConnectionAtoms(0)

        # C (phi) C (omega) N (psi) C (phi) C (omega) N connection
        # s0      s1        s2      m2      m1        m0
        # an N to C terminus connection
        pointsB = self.findCoilPoints(ConnectionB[0], 
                                      ConnectionB[1], 
                                      ConnectionB[2], 
                                      self.CNbondLength,
                                      self.CCbondLength,
                                      self.CNbondLength,
                                      self.omega,
                                      self.angleN,
                                      self.psi,
                                      self.angleC,
                                      self.phi,
                                      self.angleC)

        pointB = pointsB[2]

        if self.dumpInterimFiles == 1:
            fIO.saveXYZList([pointA, pointB], ['Ca', 'O'], "labEndPoints.xyz")

        print "Generating Coil"
        # generate a spidroin coil between each terminus.
        minDist = 1.0
        numCrankMoves = 0
        envelopeList = ['innersphere ' + str(0.9 * self.spidroinTerminusSeparation), 'frustum ' + str(self.SpidroinFrustumZ1 + 2.0) + ' ' + str(self.SpidroinFrustumMaxRadius) + ' ' + str(self.SpidroinFrustumZ2) + ' ' + str(self.SpidroinFrustumMinRadius)]
        visualiseEnvelopeEnvelope = 400.0 
        
        # build building block and dump to file
        self.spidroinHairpin = self.SpidroinCoilGen.generateBuildingBlock(self.species, pointA, pointB, minDist, numCrankMoves, envelopeList=envelopeList, visualiseEnvelope=(1000000, visualiseEnvelopeEnvelope, 'envelope_' + str(self.species) + '.xyz'))

        if self.dumpInterimFiles==1:
            fIO.saveXYZ(self.spidroinHairpin.blockXYZVals, self.spidroinHairpinAtomName, "hPinCG.xyz")

        print "Coarse grained hairpin done"

        # populate each section of the coarse grained spidroin hairpin with a peptide hairpin between the start and end points. 
        numGUnits = self.SP2NumGUnits
        lenGUnit = self.SP2LenGUnitRes
        lenPQUnit = self.SP2LenPQUnitRes
        if self.species=='SP1':
            numGUnits = self.SP1NumGUnits 
            lenGUnit = self.SP1LenGunitRes
            lenPQUnit = self.SP1LenPQunitRes
            
        radii = [minDist * self.envRadiusFactor] * numGUnits * 2
        residues = [lenPQUnit, lenGUnit] * numGUnits
        if self.species=='SP2':
            residues = residues + [lenPQUnit]
            radii = radii + [minDist * self.envRadiusFactor]
        
        pointsA = [ xyzVal for xyzVal,name in zip(self.spidroinHairpin.blockXYZVals, self.spidroinHairpin.blockAtomNames) if name == 'N' ]
        pointsB = [ xyzVal for xyzVal,name in zip(self.spidroinHairpin.blockXYZVals, self.spidroinHairpin.blockAtomNames) if name == 'C' ]
            
        if self.dumpInterFiles ==1:
            fIO.saveXYZList(pointsA + pointsB, ['Ca'] * len(pointsA) + ['O'] * len(pointsB), "labPoints.xyz")
            
        self.COfChainsBB = CCGen.generateBuildingBlock(residues, pointsA, pointsB, radii, minDist, visualiseEnvelope=(0, 50, 'envelope.xyz'))
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZ(self.COfChainsBB.blockXYZVals, self.spidroinHairpinAtomName, "hPin.xyz")

        print "detailed hairpin done"

            
        # assemble the components into a single final block of xyz values
        spidroinXYZ = self.NTerminusBackbone.blockXYZVals
        spidroinXYZ = np.concatenate( (spidroinXYZ, self.COfChainsBB.blockXYZVals), 0 )
        spidroinXYZ = np.concatenate( (spidroinXYZ , self.CTerminusBackbone.blockXYZVals), 0 )

        if self.dumpInterimFiles==1:
            fIO.saveXYZ(spidroinXYZ, self.spidroinAtomName, "spidroinAsBlock.xyz")

        print "Spidroin Done"
        return spidroinXYZ

    def findCoilPoints(self, s0, s1, s2, displacement_s2_m2, 
                                         displacement_m2_m1, 
                                         displacement_m1_m0, 
                                         alpha_s0_s1_s2_m2, 
                                         beta_s1_s2_m2, 
                                         alpha_s1_s2_m2_m1, 
                                         beta_s2_m2_m1, 
                                         alpha_s2_m2_m1_m0,
                                         beta_m2_m1_m0):
        
        # function identifies the positions of the ends of a free coil in terms of the
        # end point it is connecting to and certain bond information. 
        # Our objective is to define m2, m1 and m0 in terms of this information.
        # the names s0, s1, s2 and m2, m1, m0 are defined in the building block function
        # for connecting two blocks together.
        
        # m2 atom in terms of s0, s1 and s2
        TNB1 = coords.constructTNBFrame(s0, s1, s2)
        m2 = s2 + displacement_s2_m2 * coords.generateTNBVecXYZ(TNB1, beta_s1_s2_m2, alpha_s0_s1_s2_m2)

        # m1 atom in terms of s1, s2 and m2
        TNB2 = coords.constructTNBFrame(s1, s2, m2)
        m1 = m2 +  displacement_m2_m1 * coords.generateTNBVecXYZ(TNB2, beta_s2_m2_m1, alpha_s1_s2_m2_m1)

        # m0 atom in terms of s2, m2 and m1
        TNB3 = coords.constructTNBFrame(s2, m2, m1)
        m0 = m1 +  displacement_m1_m0 * coords.generateTNBVecXYZ(TNB3, beta_m2_m1_m0, alpha_s2_m2_m1_m0)
        
        return [m0, m1, m2]
    
    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockNames(self):
        # name the components
        spidroinNames = self.NTerminusBackbone.blockAtomNames
        spidroinNames += self.COfChainsBB.blockAtomNames
        spidroinNames += self.CTerminusBackbone.blockAtomNames
        
        return spidroinNames

if __name__=="__main__":
    
    filename = sys.argv[1]
    
    # individual spidroins
    spidroinProteinGenerator = spidroinProteinGenerator(filename)
    species1 = 'SP1'
    species2 = 'SP2'
    startPoint = np.array([0.0, -0.0, 0.0])
    direction = np.array([-0.0, 0.0, 1.0])
    rotation = 0.0
    minDist = 1.0
    
    SpidroinSpecies1BB = spidroinProteinGenerator.generateBuildingBlock( species2,
                                                                         minDist,
                                                                         showBlockDirector=False)  
    SpidroinSpecies1BB.exportBBK('species2Spidroin.xyz')

    print("Spidroin Species 2 Done.")
        
        
    