import sys
import numpy as np
import random as rnd
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.peptideBackbonePDB import pdbPeptideBackboneGenerator as PDBGen
from SpacePack.SpacePack import SpacePackBBG as SPBBG
from Library.peptideHairpin import peptideHairpinGenerator as PHG
from Projects.betasheet2 import betasheetGen as BSG


class spidroinProteinGenerator(BBG):
    # Each spidroin protein consists of two alpha helical termini and a main 
    # body that consists of regions of beta sheet structures linked 
    # by hairpin turns.

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
                 
        # Chain parameters
        self.numResiduesLongChain = self.getParam('numResiduesLongChain')
        self.numResiduesShortChain = self.getParam('numResiduesShortChain')
        self.numBetaStrands = self.getParam('numBetaStrands')
        self.betaStrandLength = self.getParam('betaStrandLength')
        self.betaStrandSeparation = self.getParam('betaStrandSeparation')
        self.betaStrandLengthPerResidue = self.getParam('betaStrandLengthPerResidue')
        
        # Over all parameters used to describe the shape and size of the spidroin packing envelope.
        self.spidroinTerminusSeparation = self.getParam('spidroinTerminusSeparation')
        self.SpidroinFrustumMaxRadius = self.getParam('SpidroinFrustumMaxRadius')
        self.SpidroinFrustumMinRadius = self.getParam('SpidroinFrustumMinRadius')
        self.SpidroinFrustumZ1 = 0.0
        self.SpidroinFrustumZ2 = self.getParam('SpidroinFrustumZ2')
        self.betaSphereCenterZ = self.getParam('betaSphereCenterZ')
        self.betaSphereRadius = self.getParam('betaSphereRadius')
        
        # atomic level information about bond angles and so on 
        self.CCbondLength = self.getParam('CCbondLength')
        self.CNbondLength = self.getParam('CNbondLength')        
        self.phi = self.getParam('phi') * np.pi / 180.0
        self.psi = self.getParam('psi') * np.pi / 180.0
        self.omega = self.getParam('omega') * np.pi / 180.0
        self.angleN = self.getParam('angleN') * np.pi / 180.0
        self.angleCA = self.getParam('angleCA') * np.pi / 180.0
        self.angleC = self.getParam('angleC') * np.pi / 180.0  
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')       
        self.verbose = self.getParam('verbose')
         
        # for colouring regions distinct colours in blender:
        self.NTerminalAtomName = self.getParam('NTerminalAtomName') 
        self.CTerminalAtomName = self.getParam('CTerminalAtomName')
        self.betaStrandAtomName = self.getParam('betaStrandAtomName')
        self.hPinAtomName = self.getParam('hpinAtomName')
        self.spidroinAtomName = self.getParam('spidroinAtomName')
        
        # create instances of all the building block generators that we will need. 
        self.BetasheetG = BSG(self.paramFilename)
        self.CoilGen = PHG(self.paramFilename)
        self.SPBBG = SPBBG(self.paramFilename)
        self.CTerminusGen = PDBGen(self.CTerminalFilename)
        self.NTerminusGen = PDBGen(self.NTerminalFilename)

        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for Spidroin Object"
            sys.exit()        

    def generateBuildingBlock(self, numBetasheets, minDist, longChain=False, showBlockDirector=False, nameByBuildingBlockType=True):

        # compute number of residues needed in the chain. 
        # Full Short length including termini is 3129 
        # Full Long length including termini is 3779
        # 
        if longChain:
            self.numResiduesCoil = self.numResiduesLongChain
        else:
            self.numResiduesCoil = self.numResiduesShortChain
        
        self.numPoints = self.numResiduesCoil * 3
        self.numBetasheets = numBetasheets
        self.inStrandDirectorHat = np.array([0.0, 0.0, 1.0])
        self.crossStrandDirectorHat = np.array([1.0, 0.0, 0.0])
        self.minDist = minDist
             
        self.spidroinDirector = np.array([0.0, 0.0, 1.0])
        self.spidroinRefPoint = np.array([0.0, 0.0, 0.0])
        self.nameByBuildingBlockType = nameByBuildingBlockType
        
        return BBG.generateBuildingBlock(self, self.numPoints, minDist, showBlockDirector=showBlockDirector)
    
    def generateBuildingBlockXYZ(self):
        # The spidroin model consists of two alpha helical termini proteins which are taken from PDB files
        # and a long coil of about 3000 residues between them. The long coil is a random chain that must fit inside
        # a frustum so that it can pack together in a sphere.  
        #
        # The centre of mass of the two termini are placed at +/- spidroinTerminusSeparation/2.0 on the x axis.
        print "CTerminus"
        refPosCTerm = np.array([-self.spidroinTerminusSeparation/2.0, 0.0, 0.0])
        self.CTerminus = self.CTerminusGen.generateBuildingBlock(director = self.CTermDirectorHat, showBlockDirector=False)
        self.CTerminus.transformBBToLabFrame(self.spidroinDirector, refPosCTerm, self.CTermRot)
        if self.nameByBuildingBlockType:
            self.CTerminus.replaceNames(self.CTerminalAtomName) 

        print "NTerminus"
        refPosNTerm  = np.array([self.spidroinTerminusSeparation/2.0, 0.0, 0.0])
        self.NTerminus = self.NTerminusGen.generateBuildingBlock(director = self.NTermDirectorHat, showBlockDirector=False)
        self.NTerminus.transformBBToLabFrame(self.spidroinDirector, refPosNTerm, self.NTermRot)
        if self.nameByBuildingBlockType:
            self.NTerminus.replaceNames(self.NTerminalAtomName) 
        
        # accumulate a pointsToAvoid array as we go.
        pointsToAvoid = self.NTerminus.blockXYZVals + self.CTerminus.blockXYZVals 
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZ(self.NTerminus.blockXYZVals, self.NTerminalAtomName, "NTerminal.xyz")
            fIO.saveXYZ(self.CTerminus.blockXYZVals, self.CTerminalAtomName, "CTerminal.xyz")
        
        
        print "constructing BetaSheets"
        numLoopResidues = 0
        betaSheetOffset = np.array([0.0, 0.0, 0.0])

        # construct correct number of antiparallel betasheet building blocks without any loops
        self.betaSheetBBs = [ self.BetasheetG.generateBuildingBlock( self.numBetaStrands, 
                                                                     self.betaStrandLength,
                                                                     numLoopResidues,
                                                                     self.minDist, 
                                                                     self.inStrandDirectorHat,
                                                                     self.crossStrandDirectorHat,
                                                                     betaSheetOffset,
                                                                     parallel=False) for n in range(0, self.numBetasheets) ]
    
        # compute the longest distance available in the beta sheet. 
        # This is the mindist between the centre of mass of two beta sheets 
        # and twice the distance of the betasheet centre of mass from the geometric boundary.
        # add two to the betaStrandLength to cope with the connection points
        betasheetSeparation = np.sqrt( ( (self.betaStrandLength + 2) * self.betaStrandLengthPerResidue)**2 + (self.numBetaStrands * self.betaStrandSeparation) **2 )  
        
        # inhibit beta sheets from being closer than betaSheetSeparation/2 from the boundary where the coil cannot go.          
        betaEnvelope = ['betasphere ' + str(self.betaSphereCenterZ) + ' ' + str(self.betaSphereRadius)]
        
        # envelopeList=['None'] # useful line to have around to override envelope for debugging
        
        # compute a large box which surrounds the beta sheet zone
        XRange = [-1.5 * self.betaSphereRadius, 1.5 * self.betaSphereRadius]
        YRange = [-1.5 * self.betaSphereRadius, 1.5 * self.betaSphereRadius]
        ZRange = [self.betaSphereCenterZ - 1.5 * self.betaSphereRadius, 0]

        # calculate the positions of the centre of masses of the beta sheets within the beta envelope
        betaSheetCOMBB = self.SPBBG.generateBuildingBlock( self.numBetasheets,
                                                           XRange,
                                                           YRange,
                                                           ZRange, 
                                                           betasheetSeparation,
                                                           visualiseEnvelope=(10000, 2 * self.betaSphereRadius),
                                                           envelopeList=betaEnvelope)
        
        betaSheetDirectors =[]
        
        # create beta sheet directors
        for n in range(0, self.numBetasheets):
            theta, phi = coords.pickRandomPointOnUnitSphere()
            betaSheetDirectors.append(coords.sphericalPolar2XYZ(np.array([1.0, theta, phi])))

        print "Transform Beta Sheet locations"
        for director, com, betaSheetBB in zip(betaSheetDirectors, betaSheetCOMBB.blockXYZVals, self.betaSheetBBs):
            betaSheetBB.transformBBToLabFrame(director, com, 0) 

        # rename the atoms in the betasheets if the flag is set
        if self.nameByBuildingBlockType:
            for bsheet in self.betaSheetBBs:
                bsheet.replaceNames(self.betaStrandAtomName)

        # add the beta sheets to the pointsToAvoid List
        for betaSheetBB in self.betaSheetBBs:
            pointsToAvoid = np.concatenate( (pointsToAvoid, betaSheetBB.blockXYZVals), 0 )

        if self.dumpInterimFiles==1 and self.numBetasheets>0:
            # compile the beta sheets into a single entity and dump to file - for debug
            betaSheets = self.betaSheetBBs[0].blockXYZVals
            for betaSheet in self.betaSheetBBs[1:]:
                betaSheets = np.concatenate( (betaSheets, betaSheet.blockXYZVals), 0 )
            fIO.saveXYZ(betaSheets, self.betaStrandAtomName, "betaSheet.xyz")


        
        print "Constructing Coils"

        # Compute the coil start (A) and end (B) points. 
        # Define points m2, m1 and m0 for the s2, s1 and s0 end points of the coil.
        # ensure a realistic join. There is enough space at the end of the termini to not check for clashes.
        # The betasheetseparation is increased by 1 residue longer than the strand length to help with this.
        # find the C-terminus end of the N terminus and calculate where the coil should start
        ConnectionA = self.NTerminus.getConnectionAtoms(1)

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
        
        # find the N-terminus end of the C terminus and calculate where the coil should start
        ConnectionB = self.CTerminus.getConnectionAtoms(0)

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

        betaSheetNPoints = []
        betaSheetCPoints = []
        # Determine the connector atoms that define the seed points for each coil.
        # and make a list of the seed points for each free end of coil.
        # We make an assumption that there are no other beta strand in the vicinity of each beta stranf.
        # We can help to ensure this by defining an outer sphere when we come to make each coil link to inhibit
        # the coils from straying too far from the direct path between the end points.
        for betaSheetBB in self.betaSheetBBs:
            # get a list of all the connectors for the current beta sheet.
            # always alternating with N then C terminal connectors depending on how many strands 
            for n in range(0, len(betaSheetBB.getAllConnectionIndices())): 
                connector = betaSheetBB.getConnectionAtoms(n)
                if n % 2==1:
                    # C terminus of beta strand.
                    # compute m0, m1 and m2 in a sub list and append that to the list of coil points
                    # N (psi) C (phi) C (omega) N (psi) C (phi) C connection
                    # s0      s1        s2      m2      m1        m0
                    betaSheetCPoints.append(self.findCoilPoints(connector[0], 
                                                                connector[1], 
                                                                connector[2], 
                                                                self.CNbondLength,
                                                                self.CNbondLength,
                                                                self.CCbondLength,
                                                                self.phi,
                                                                self.angleC,
                                                                self.omega,
                                                                self.angleN,
                                                                self.psi,
                                                                self.angleC))
                else:
                    # N Terminus of beta strand
                    # C (phi) C (omega) N (psi) C (phi) C (omega) N connection
                    # s0      s1        s2      m2      m1        m0
                    betaSheetNPoints.append(self.findCoilPoints( connector[0], 
                                                                 connector[1],
                                                                 connector[2], 
                                                                 self.CNbondLength,
                                                                 self.CCbondLength,
                                                                 self.CNbondLength,
                                                                 self.omega,
                                                                 self.angleN,
                                                                 self.psi,
                                                                 self.angleC,
                                                                 self.phi,
                                                                 self.angleC))
        
        # seed the coil point list with the global N and C termini of the coil 
        coilPoints = [(pointsA, pointsB)] # coil pair goes from C terminus to N terminus 
        # so pair[0] is C, pair[1] is N

        
        if self.numBetasheets>0:
        
            # find a valid order of cyling through the connectors
            COrderIndex, NOrderIndex = self.findValidPairOrder(self.numBetasheets, self.numBetaStrands)
            
            # Re order N and C points with the chosen valid index order
            betaSheetNPoints = [ betaSheetNPoints[index] for index in NOrderIndex ] 
            betaSheetCPoints = [ betaSheetCPoints[index] for index in COrderIndex ]
    
            
            # loop through each CPoint, NPoint pair (which will be randomised throughout all the 
            # beta sheets) to make a master list of connections to connect with a randomcoil
            for CPoint, NPoint in zip(betaSheetCPoints, betaSheetNPoints):
                newCoilEntry = (coilPoints[-1][0], NPoint) # the previous C terminus going to a new N terminus.
                coilPoints[-1] = (CPoint, coilPoints[-1][1]) # a new C terminus going to the last N terminus.
                coilPoints.insert(-1, newCoilEntry) # insert the new coil points at the end of the list.

        if self.dumpInterimFiles==1:
            coilPointsXYZ = []
            coilPointNames = []
            curConnection = 0
            names = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'Ne', 'K'] # define 10 atoms for 10 colours and then cycle colours 
            
            # dump all the coilPoints. Colour each pair it's own color to help ID the connectors.
            for coilPair in coilPoints:
                for coilConnector in coilPair:
                    coilPointsXYZ += coilConnector
                    coilPointNames += ['H', 'He', 'Li']
                    # coilPointNames += [names[curConnection % 10], names[curConnection % 10], names[curConnection % 10]]
                curConnection += 1

            fIO.saveXYZList(coilPointsXYZ, coilPointNames, "CoilConnectors.xyz")

        print "Generating Coils"

        # The total length of the chain is to be broken up into chunks, each of which 
        # must be long enough to reach the next point. We take betaStrandLengthPerResidue as the minimum 
        # length of residue and compute the approximate number of residues needed to reach 
        # between each pair of points in a straight line.
        spatialDistanceBetweenPoints = [ np.linalg.norm(pair[0][2] - pair[1][2]) for pair in coilPoints ] 
        minVals = [ 2 * int(np.ceil(length/self.betaStrandLengthPerResidue)) for length in spatialDistanceBetweenPoints ]
        
        numResiduesInBetaStrands = numBetasheets * self.numBetaStrands * self.betaStrandLength
        numResiduesInTermini = int(float(self.NTerminus.countAtoms())/3.0 + float(self.CTerminus.countAtoms())/3.0) 
        numResiduesToDivideUp = self.numResiduesCoil - numResiduesInBetaStrands - numResiduesInTermini  
        
        if numResiduesToDivideUp < sum(minVals):
            print "Warning: Not enough residues for minimal coil connections."
        
        coilLengths = self.divideTotalSumEvenlyAmongListOfGroups(numResiduesToDivideUp, minVals)
        numCrankMoves = 0
        
        # generate a random coil between each terminus adding the final loop to pointsToAvoid.
        # make sure it's inside the global frustum and also inside the local sphere defined by the distance
        # between the connection points.
        self.hairpinBBs = []
        for coil, coilLength, distance in zip(coilPoints, coilLengths, spatialDistanceBetweenPoints):
            envelopeList = ['frustum ' + str(self.SpidroinFrustumZ1) +' ' + str(self.SpidroinFrustumMaxRadius) + ' ' + str(self.SpidroinFrustumZ2) + ' ' + str(self.SpidroinFrustumMinRadius) ]
            envelopeList.append('innersphere ' + str(0.9 * distance/2.0))
            if len(coilPoints)>1:
                envelopeList.append('outersphere ' + str(1.5 * distance/2.0))
           
            #envelopeList=['None'] # useful debug statement    
            # generate the hairpin connection
            hairpinBB = self.CoilGen.generateBuildingBlock(coilLength, 
                                                           coil[0],
                                                           coil[1],
                                                           self.minDist,
                                                           numCrankMoves,
                                                           visualiseEnvelope=(0, 300), 
                                                           pointsToAvoid=[], 
                                                           envelopeList=envelopeList)
            
            # add the hairpin to the pointsToAvoid array
            pointsToAvoid = np.concatenate( (pointsToAvoid, hairpinBB.blockXYZVals), 0 )
            self.hairpinBBs.append(hairpinBB)

        if self.nameByBuildingBlockType:
            [ hpin.replaceNames(self.hPinAtomName) for hpin in self.hairpinBBs]

        if self.dumpInterimFiles==1:
            hPins = self.hairpinBBs[0].blockXYZVals
            for hPin in self.hairpinBBs[1:]:
                hPins = np.concatenate( (hPins, hPin.blockXYZVals), 0 )
            fIO.saveXYZ(hPins, self.hPinAtomName, "hPin.xyz")
            
        print "hairpins done"
        
        # assemble the components into a single final block of xyz values
        spidroinXYZ = self.NTerminus.blockXYZVals
        for betaSheet in self.betaSheetBBs:
            spidroinXYZ = np.concatenate( (spidroinXYZ, betaSheet.blockXYZVals), 0 )
        for hPin in self.hairpinBBs:
            spidroinXYZ = np.concatenate( (spidroinXYZ, hPin.blockXYZVals), 0 )
        spidroinXYZ = np.concatenate( (spidroinXYZ , self.CTerminus.blockXYZVals), 0 )

        if self.dumpInterimFiles==1:
            fIO.saveXYZ(spidroinXYZ, self.spidroinAtomName, "spidroinAsBlock.xyz")

        print "Spidroin Done"
        return spidroinXYZ

    def findValidPairOrder(self, N, M):
        # function to find a valid random ordering of the beta sheet connectors
        CList = []
        NList = []
        
        # generate a list of all the end points of all the strands 
        for n in range(N):  # n is number of betaSheets
            for m in range(M): # m is number of strands per beta sheet
                CList.append((n, m))
                NList.append((n, m))
        
        CListShuffled = CList[:] 
        NListShuffled = NList[:]
                
        while not self.checkOpenLoop(CListShuffled, NListShuffled):
            rnd.shuffle(CListShuffled)
            rnd.shuffle(NListShuffled)

        # write down the indices in the original CList in the order they appear in the shuffled list.
        COrder = [ CList.index(pair) for pair in CListShuffled[1:] ]
        # add the first ClistShuffled pair to the end of the CList index order
        COrder.append(CList.index(CListShuffled[0]))
        
        # start the Nlist at the 0th entry in the NListShuffled list
        NOrder = [ NList.index(pair) for pair in NListShuffled ]

        return COrder, NOrder
    
    def checkOpenLoop(self, CList, NList):
        # This function returns False if CList and NList between them  
        # contain instructions that will form a closed loop before visiting all the CList entries 
        # except the first one.
        #
        # The entries in the CList and Nlist are corresponding (sheet, strand) termini 
        # that will be connected by a coil.  So the CList[j] = (p, q) will 
        # be connected to NList[j] = (n, m).
        # 
        # The other end of the strand referred to in NList[j] = (n, m) will appear in
        # the Clist at index k. So NList[j] = CList[k] = (n,m). 
        # Thus we can find the next index k that we have to visit in the CList by looking up (n,m)
        # in the CList, which leads to an NList[k] = (s, t) etc etc.
        # We keep following this chain and if we revisit an index in the Clist that we've 
        # been to before then we have formed a closed loop.
        #
        # We pick the first entry in the NList as the starting point. This is the N terminal
        # that the C terminal of the N Terminal domain will connect to.  The CList entry in position 0 
        # will therefore be the last beta sheet C terminus of the chain that will connect to the N Terminus of 
        # the global C terminal domain. 
        # 
        
        # pick the first index in the NList
        currentIndex = 0
        
        # create an array to store the indices of the CPoints we have visited already
        CPointsVisited = []
        
        openLoop = True # assume success
        
        # we must visit each entry in the CPoints list except one of them once and once only. 
        while len(CPointsVisited) < len(CList) - 1:
            # get the (sheet, strand) pair of the current Index and 
            # look up the next index to check in the CList.
            nextIndex = CList.index(NList[currentIndex])
            
            # if this index has not been visited before then add it to the visited list
            if not nextIndex in CPointsVisited:
                # add the index k to the CPointsVisited List
                CPointsVisited.append(nextIndex)
                # set the current index to the next index
                currentIndex = nextIndex
            else:
                # we have visted this index before and therefore this list ordering
                # contains a closed loop which we do not want
                openLoop = False
                break
            
        return openLoop


                

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
    
    def divideTotalSumEvenlyAmongListOfGroups(self, totalSum, minVals):
        # function takes a list of integers called minVals
        # and adds numbers to each item in proportion to the value in the list.
        
        totalMinVals = sum(minVals)
        
        
        output= [ int(np.floor(float(totalSum) * float(minVal)/float(totalMinVals))) for minVal in minVals ]
        
        remainder = totalSum - sum(output)
        
        # distribute the remainder evenly among the output
        curIndex = 0
        while remainder > 0:
            output[curIndex]  += 1
            remainder -= 1
            curIndex+=1
            if curIndex>=len(minVals):
                curIndex=0
        
        return output
    
    def divideTotalSumRandomlyAmongListOfGroups(self, totalSum, minVals):
        # function takes a list of integers called minVals
        # and randomly adds numbers to each item in the list until the sum equals totalSum 
        remainingSum = totalSum - sum(minVals)
        
        outputList = minVals[:]
        while remainingSum>0:
            if remainingSum==1:
                valToAdd = 1
            else:
                valToAdd = rnd.randrange(0, stop=remainingSum)
            valToAddTo = rnd.randrange(0, stop=len(minVals))
            outputList[valToAddTo] += valToAdd
            remainingSum = totalSum - sum(outputList)
            
        return outputList
    
    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0])

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockNames(self):
        # name the components
        spidroinNames = self.NTerminus.blockAtomNames
        for betaSheet in self.betaSheetBBs:
            spidroinNames += betaSheet.blockAtomNames
        for hpin in self.hairpinBBs:
            spidroinNames += hpin.blockAtomNames
        spidroinNames += self.CTerminus.blockAtomNames
        
        return spidroinNames

if __name__=="__main__":
    noErrors = True
    
    filename = sys.argv[1]
    
    # generating an individual spidroin
    spidroinProteinGenerator = spidroinProteinGenerator(filename)
    numBetasheets = 0
    minDist = 1.0
    longChain = False
    startPoint = np.array([0.0, -0.0, 0.0])
    direction = np.array([-0.0, 0.0, 1.0])
    rotation = 0
    
    SpidroinBB = spidroinProteinGenerator.generateBuildingBlock(numBetasheets, 
                                                                minDist, 
                                                                longChain=longChain, 
                                                                showBlockDirector=False, 
                                                                nameByBuildingBlockType=False)  
    SpidroinBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print("Done.")
        
        
    