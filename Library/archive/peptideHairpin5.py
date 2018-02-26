import sys
import numpy as np
import random as rnd
import time
import Utilities.fileIO as fIO
import Utilities.cartesian as cart
import Utilities.coordSystems as coords
from Library.constrainedPolymer5 import ConstrainedPolymerPackBBG as CPBBG
from Library.peptideBackbone import peptideBackboneGenerator as PBG

class peptideHairpinGenerator(CPBBG):
    # This class returns a randomly coiled peptide backbone between two sets of points 
    # A and B, which is packed into an external envelope. The backbone does not
    # self-intersect and the coil avoids a user supplied set of external points.  
    #
    # The class inherits the functionality of the constrainedPolymerClass.
    # but overides the generateSpaceCurve method of that class to provide the initial unfolded polymer.
    #
    # Instead of generating a geometric shape for the polymer, the function 
    # uses the peptide backbone class to generate a straight peptide backbone model with 
    # appropriate dihedrals, bond angles, bondLengths and number of residues, 
    # with the three points specified in pointsA forming the initial seed for the chain.
    # The first point in A is non negotiable and fixed. But the second and third points
    # are movable once the initial chain has been built.
    # PointsA is supplied in the format m0, m1 and m2 where m0 is the innermost point of a connector.
    # Thus point m2 (the third point in pointsA) is actually the start of the chain. This is important.
    # Throughouth the project connectors are defined as m0, m1 and m2.
    #
    #
    # The energy of two springs between the two atoms at the free end of the chain and 
    # the two atoms of the second anchor points B is 
    # then computed based on the distance between them. Conformations with smaller 
    # amounts of energy are found by picking C-CA or N-CA bonds at random and 
    # performing dihedral moves of the entire free end of chain from the selected 
    # bond onwards. 
    #
    # Each time a move is made the new energy of the two springs is calculated for the 
    # new position. If the move results in lower total spring energy than the lowest energy 
    # yet found it is accepted. If the spring energy is larger than the current minimum 
    # the move is accepted with a probability based on the boltzman factor between the 
    # current lowest energy and the new energy.
    #
    # As the distance between the end point and point B shrinks the allowed size 
    # range of rotations shrinks rapidly.
    #    
    # The naming functions and allowedList functions are also overridden so
    # that only the alpha carbons are allowed to be used in the dihedral twist 
    # rotations. This restriction preserves the peptide bond angles in the chain.
    # 
    # Once a double anchored space curve is generated, with anchor points at 
    # pointsA and B, the procedure for randomising the coil with crankshaft moves and 
    # folding the chain into an envelope is the same as for the coiled polymer baseclass.
    #
    # The final structure is exported as a building block class.  

      
    def __init__(self, paramFilename):
        # initialise the parameter dictionary for the base classes
        CPBBG.__init__(self, paramFilename)
        
    def initialiseParameters(self):
        # initialise the constrained polymer parent
        CPBBG.initialiseParameters(self)
        
        self.distEpsilon= self.getParam("distEpsilon")
        self.maxNumConnectingMoves = self.getParam("maxNumConnectingMoves")
        self.springConstant = self.getParam("springConstant")
        self.connectingTemp = self.getParam("connectingTemp")
        

        # load the backbone building object
        self.PBG = PBG(self.paramFilename)
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for hairpin"
            sys.exit()        

    def generateBuildingBlock(  self, 
                                numResidues, 
                                pointsA,
                                pointsB,
                                minDist,
                                numCrankMoves,
                                pointsToAvoid=[], 
                                envelopeList=['None'],
                                nameCA=False):
        self.numResidues = numResidues
        self.numPoints = self.numResidues * 3
        bondLength = 0.0 # not used in this game
        self.nameCA = nameCA
        self.minDist = minDist

        # set up the right reference points
        self.labPointA = pointsA[2]
        self.labPointB = pointsB[2]
        self.labPointsA = pointsA[:]
        self.labPointsB = pointsB[:]
        
        # generate the BuildingBlock reference point earlier than usual because 
        # we need the transformation for the pointsToAvoid input.
        self.blockRefPoint = self.generateBuildingBlockRefPoint()
        
        # generate the BuildingBlock director unit vector earlier than usual because 
        # we need the transformation for the pointsToAvoid input.
        self.blockDirectorHat = self.generateBuildingBlockDirector()        
        
        # generate the transformation information from building block to labPointA and labPointB
        self.labDirector, self.labRefPoint, self.labRotation = self.computeTransform()
        
        # convert the pointsA and pointsB information from the labFrame to the block frame
        self.blockPointsA = coords.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, self.blockDirectorHat, self.blockRefPoint, pointsA)
        self.blockPointsB = coords.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, self.blockDirectorHat, self.blockRefPoint, pointsB)        
        self.pointsToAvoid = coords.transformFromLabFrameToBlockFrame(self.labDirector, self.labRefPoint, self.labRotation, self.blockDirectorHat, self.blockRefPoint, pointsToAvoid)

        if self.dumpInterimFiles==1:
            fIO.saveXYZList(self.blockPointsA + self.blockPointsB, ['Ca', 'Ca', 'Ca', 'O', 'O', 'O'], "blockPointsAB.xyz")
        
        # parse the envelope list if we intend to use it.
        self.parseEnvelopeList(envelopeList)
        
        # check starting points are legal or it's gonna be a long wait.
        for pos in self.blockPointsA:
            if not self.checkPointInBounds(pos):
                print "Error Warning: PointA out of bounds"
                time.sleep(3)
             
        for pos in self.blockPointsB:
            if not self.checkPointInBounds(pos):
                print "Error Warning: PointB out of bounds"
                time.sleep(3)
        
        return CPBBG.generateBuildingBlock(self, self.numPoints, pointsA[2], pointsB[2], minDist, bondLength, numCrankMoves, pointsToAvoid=self.pointsToAvoid, envelopeList=envelopeList) 

    def generateAllowedList(self):
        # add the first and last points to the allowed list.
        # This is important to allow the full chain to be able to be mapped into the envelope.
        allowed = [0, self.numPoints-1]
        
        # over ride the naming of the atoms with ca=True. 
        names = self.generateBuildingBlockNames(ca=True)
        [ allowed.insert(-1,i) for i, name in enumerate(names) if name=='CA' ]
        
        return allowed

    def generateSpaceCurve(self):
        # Over-rides the generate space curve function of the parent to generate a peptide backbone
        # and a pseudo energy landscape approach to find an initial chain with the end point fixed
        # at point B.
        
        # create a regular b0a14ckBone using block Points A as the first residue
        peptideBackbone = self.PBG.generateBuildingBlock(self.numResidues, seedResidue = self.blockPointsA)
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZList(peptideBackbone.blockXYZVals, peptideBackbone.blockAtomNames, 'initialPeptideBackbone.xyz')
        
        # extract the xyzValues
        xyzVals = peptideBackbone.getAtomsXYZ()

        # perform the energy minimisation that moves the free end to blockPointsB
        xyzVals = self.minimiseEnergy(xyzVals, self.allowedList)

        if self.dumpInterimFiles==1:
            fIO.saveXYZ(xyzVals, 'K', 'chainConnectedBlockFrame.xyz')

        return xyzVals  
        

    def minimiseEnergy(self, xyzVals, allowedList):
        # Returns a polymer chain which minimises a simple PE function.
        # Performs random dihedral twists on the free end of a polymer, 
        # starting at a random bond.
        # Moves resulting in lower energy arrangements are accepted.
        # Moves resulting in higher energy arrangements are accepted
        # with probability that is exponentially smaller with increasing energy difference.
        # Since the primary structure is just a straight line and we are only
        # doing sparse dihedral twists, with only a single bias towards to the pointB,
        # there is a low probability of a self-intersection especially for long chains.
        # This rapidly finds a structure where the point B is within arbitrary distance of pointB.
        # The step size scales in proportion to the distance from point B. So only tiny steps are taken 
        # near to the only minimum of the entire potential. Converges rapidly even for large N.
         
        lowestEnergyMinimum = xyzVals[:]
        initPE, initDist = self.PE(xyzVals)
        curPE = initPE
        minPE = initPE
        curDist = initDist
        minDist = initDist
        maxStepRange = 1.0
        numMoves = 0
        curMin = 0
        while minDist> self.distEpsilon and numMoves < self.maxNumConnectingMoves:

            # compute new conformation based on a random dihedral twist
            newXYZ = self.dihedralTwist(xyzVals, maxStepRange)

            # compute energy and distance of new move 
            newPE, newDist = self.PE(newXYZ)   

            # compute energy difference with current minimum PE
            deltaPE = newPE - minPE

            # assume we will accept the move.   
            acceptMove = True 
            # if the currentPE is greater than the minimum then only accept
            # the move with a probability given by the difference in 
            # energy of the minimum and current states  
            if deltaPE > 0 :
                # pick a random value between 0 and 1
                prob = rnd.uniform(0,1)

                # if that value is larger than the threshold reject the move.
                # The threshold decreases with increasing deltaE, so the 
                # higher the energy of the new state relative to the old one
                # the more likely it is we reject the move 
                if prob > np.exp(-deltaPE/self.connectingTemp):
                    acceptMove = False
                
            # if we accept the move then store the new coords
            # and record the energy of the newest accepted move    
            if acceptMove:
                xyzVals = newXYZ[:]
                curPE = newPE
                curDist = newDist
                
            # check the curPE against the minimum energy
            # if we have a new min energy then retain for the future and dump to file.
            if curPE < minPE:
                lowestEnergyMinimum = xyzVals[:]
                minPE = curPE
                minDist = curDist
                maxStepRange = minDist/initDist
                curMin += 1 
                self.outline(numMoves, self.maxNumConnectingMoves, minDist, minPE, maxStepRange )
            numMoves += 1
            
            if numMoves % 100 == 0: 
                self.outline(numMoves, self.maxNumConnectingMoves, minDist, minPE, maxStepRange )                 
        
        return lowestEnergyMinimum

    def outline(self, n, M, d, E, R):
        print n, "out of ", M, "minDist:", d, "minEnergy:", E, "maxStepRange:", R
    
    def PE(self, xyzVals):
        PE = 0.0
        # add spring between pointB[2] and end of xyzVals with equilibrium at pointB
        dist1 = np.linalg.norm(xyzVals[-1] - self.blockPointsB[2])
        PE += 0.5 * 3 * self.springConstant * dist1**2
        
        # add spring between pointsB[1] and end of xyzVals[-2] 
        dist2 = np.linalg.norm(xyzVals[-2] - self.blockPointsB[1])
        PE += 0.5 * 2 * self.springConstant * dist2**2

        # add spring between pointsB[0] and end of xyzVals[-3] 
        #dist3 = np.linalg.norm(xyzVals[-3] - self.blockPointsB[0])
        #PE += 0.5 * self.springConstant * dist3**2
        
        # made the first spring the stiffest with second and third springs getting weaker.
        # first spring has largest influence on energy.
        
        
        # add weak spring between each CA and refPoint
        # and repulsion at close range
        #for i in allowedList:
        #    posDist = np.linalg.norm(xyzVals[i] - self.blockRefPoint)
        #    PE += +0.5 * 2 * self.springConstant * posDist **2
        #    #PE += -1.0/posDist**6
        
        # add contribution to PE based on distance between each CA with all the others
        # only look at ordered pairwise combinations without repetition.
        # ie. 1-2, 1-3, 1-4, 2-3, 2-4, 3-4  (leaving out 1-1, 2-2, and 3-1, 4-1 etc) 
        #pairWiseCombinations = it.combinations([ pos for n, pos in enumerate(xyzVals) if n in allowedList], 2)
        # iterate and compare each pair to compute contrib to PE
        #for pair in pairWiseCombinations:
        #    PE += -self.epsilon0/np.power(np.linalg.norm(pair[1] - pair[0]), 3)
     
        # add contribution to PE from every point in pointsToAvoid with 
        # every point in the list 
        #avoidList = it.product(xyzVals, self.pointsToAvoid)
        # iterate and compare each pair to see if any match
        #for pair in avoidList:
        #    PE += - self.LJRep(pair[0], pair[1])
        
        return PE, np.sqrt(dist1**2 + dist2**2)
        
    def dihedralTwist(self, xyzVals, maxStepRange):  
        
        #fIO.saveXYZList(xyzVals, self.blockNames, "preTwist.xyz")
        # can only do this if there are sufficient atoms in the array
        if len(xyzVals) > 3:
            # initialise the axisAtom1Index
            axisAtom1Index = 2
            # keep picking a random atom until we get one in the allowed list.
            while not axisAtom1Index in self.allowedList:
                axisAtom1Index = rnd.randint(0, len(xyzVals) - 3)
            
            # find the relevant points and rotation axis
            atom1 = xyzVals[axisAtom1Index]
            atom2 = xyzVals[axisAtom1Index + 1]
            rotAxis = atom2 - atom1
            rotAxisHat = rotAxis/np.linalg.norm(rotAxis)
            
            # pick a step size at random from within the current allowed range
            angle = rnd.uniform(maxStepRange * -np.pi, maxStepRange * np.pi)

            # rotate all the remaining points about the atom1 axis place at atom1 by angle
            xyzVals = [ p if n < axisAtom1Index + 2 else cart.rotPAboutAxisAtPoint(p, atom1, rotAxisHat, angle) for n, p in enumerate(xyzVals) ]

        #fIO.saveXYZ([atom1, atom2], 'O', "atomAxis.xyz")
        #fIO.saveXYZList(xyzVals, self.blockNames, "postTwist.xyz")
        return xyzVals 
        
    def LJRep(self, p1, p2):
        # add a short range repulsive term based on the distance between 
        # centre of each point
        return self.epsilon0/np.power(np.linalg.norm(p1 - p2), 12)
        
        
    def generateBuildingBlockNames(self, ca=False):
        
        if self.nameCA:
            names = ['N', 'CA', 'C'] * self.numResidues
        else:
            names = ['N', 'C', 'C'] * self.numResidues
        
        # override the externally controlled member variable nameCA with a local variable
        if ca==True:
            names = ['N', 'CA', 'C'] * self.numResidues
        
        return names

    def generateBuildingBlockConnectors(self):
        # N connector first, then C Connector 
        return [ [2, 1, 0], [self.numPoints-3, self.numPoints-2, self.numPoints-1] ]
        
if __name__ == "__main__":
    
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the generator objects.
    hairPinGen = peptideHairpinGenerator(filename)
    backboneGenerator = PBG(filename)
    
    # generate 1 residue seed building block
    seedResidue = backboneGenerator.generateBuildingBlock(1)
    
    # generate starting points and move the seed to those points extracting
    # the xyzVals each time. 
    numResidues = 20
    pointA = np.array([10.0, 10.0, 10.0])
    pointB = np.array([-10.0, 10.0, -10.0])
    
    print "Estimate min num residues: ", np.linalg.norm(pointA-pointB)/3.5
    
    seedResidue.placeAtom(2, pointA)
    seedResidue.setBlockRefPoint(pointA)
    seedResidue.orientToDirector(np.array([-1.0, 1.0, 0.0]))
    pointsA = seedResidue.blockXYZVals[:]
    
    seedResidue.placeAtom(2, pointB)
    seedResidue.setBlockRefPoint(pointB)
    seedResidue.orientToDirector(np.array([1.0, 0.0, 1.0]))
    pointsB = seedResidue.blockXYZVals[:]
    
    minDist = 1.0
    numCrankMoves = 5
    fIO.saveXYZList([pointsA[0], pointsA[1], pointsA[2], pointsB[0], pointsB[1], pointsB[2]], ['Ca', 'Ca', 'Ca', 'O', 'O', 'O'], 'labEndPoints.xyz')

    # build building block and dump to file
    hairpinBuildingBlock = hairPinGen.generateBuildingBlock(numResidues, pointsA, pointsB, minDist, numCrankMoves, envelopeList=["innersphere 8"])
    hairpinBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "hairpin done"