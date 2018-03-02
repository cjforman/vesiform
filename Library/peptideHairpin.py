import sys
import numpy as np
import time
import Utilities.fileIO as fIO
from Library.constrainedPolymer import ConstrainedPolymerPackBBG as CPBBG
from Library.peptideBackbone import peptideBackboneGenerator as PBG

class peptideHairpinGenerator(CPBBG):
    ''' This class returns a randomly coiled peptide backbone between two sets of points 
    A and B, which is packed into an external envelope. The backbone does not
    self-intersect and the coil avoids a user supplied set of external points.  
    
    The class inherits the functionality of the constrainedPolymerClass.
    The generateSpaceCurve method of that class is overridden to provide an unfolded 
    peptideBackbone.
    
    The user supplies three points for pointsA and pointsB. PointsA form the initial seed for the chain.
    The first point in A is non negotiable and fixed. But the second and third points
    are movable once the initial connected chain has been built.
    
    The main advantage of providing three points is so that derivatives as well as the end points match up,
    yielding a much smoother link between coils.
    
    PointsA and pointsB are supplied in the format m0, m1 and m2 where m0 is the innermost point of a connector.
    Thus point m2 (the third point in pointsA) is actually the start of the chain. This is important.
    Throughouth the project connectors are defined as m0, m1 and m2.
    
    The energy function is overloaded to compute the energy of two springs between the 
    two atoms at the free end of the chain and the two atoms of the second anchor; Points B is 
    then computed based on the distance between them. Conformations with smaller 
    amounts of energy are found by picking C-CA or N-CA bonds at random and 
    performing dihedral moves of the entire free end of chain from the selected 
    bond onwards. 
    
    Each time a move is made the new energy of the two springs is calculated for the 
    new position. If the move results in lower total spring energy than the lowest energy 
    yet found it is accepted. If the spring energy is larger than the current minimum 
    the move is accepted with a probability based on the boltzman factor between the 
    current lowest energy and the new energy.
    
    As the distance between the end point and point B shrinks the allowed size 
    range of rotations shrinks rapidly.
        
    The naming functions and allowedList functions are also overridden so
    that only the alpha carbons are allowed to be used in the dihedral twist 
    rotations. This restriction preserves the peptide bond angles in the chain.
     
    Once a double anchored space curve is generated, with anchor points at 
    pointsA and B, the procedure for randomising the coil with crankshaft moves and 
    folding the chain into an envelope is the same as for the coiled polymer baseclass.
    
    The final structure is exported as a building block class.'''

      
    def __init__(self, paramFilename):
        # initialise the parameter dictionary for the base classes
        CPBBG.__init__(self, paramFilename)
        
    def initialiseParameters(self):
        # initialise the constrained polymer parent
        CPBBG.initialiseParameters(self)
        
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
                                visualiseEnvelope=(0,20), 
                                envelopeList=['None'],
                                nameCA=False):
        self.numResidues = numResidues
        self.numPoints = self.numResidues * 3
        bondLength = 0.0 # not used in this game
        self.nameCA = nameCA
        self.minDist = minDist

        # set up the right reference points
        self.pointA = pointsA[2]
        self.pointB = pointsB[2]
        self.pointsA = pointsA[:]
        self.pointsB = pointsB[:]
        
        self.parseEnvelopeList(envelopeList)
        self.blockRefPoint = self.generateBuildingBlockRefPoint()
        
        # prematurely add pointsToAvoid and parse the envelope due to the pointsA and pointsB checks below
        self.pointsToAvoid = pointsToAvoid 
        
        # check starting points are legal or it's gonna be a long wait.
        for pos in self.pointsA:
            if not self.checkPointInBounds(pos):
                print "Error Warning: One of pointsA out of bounds"
                time.sleep(3)
             
        for pos in self.pointsB:
            if not self.checkPointInBounds(pos):
                print "Error Warning: One of pointsB out of bounds"
                time.sleep(3)
        
        return CPBBG.generateBuildingBlock(self, self.numPoints, pointsA[2], pointsB[2], minDist, bondLength, numCrankMoves, visualiseEnvelope=visualiseEnvelope, pointsToAvoid=pointsToAvoid, envelopeList=envelopeList) 

    def generateAllowedList(self, short=False):
        # add the first and last points to the allowed list.
        # This is important to allow the full chain to be able to be mapped into the envelope.
        allowed = [0, self.numPoints-1]
        
        # over ride the naming of the atoms with ca=True. 
        names = self.generateBuildingBlockNames(ca=True)
        
        # if the minimisation is already close then only make moves in the last ten CAs closest to the chain 
        if short:
            [ allowed.insert(-1,i) for i, name in enumerate(names) if name=='CA' and i>(len(names) - 10) ]
        else:
            [ allowed.insert(-1,i) for i, name in enumerate(names) if name=='CA' ]
        
        return allowed

    def generateSpaceCurve(self):
        # Over-rides the generate space curve function of the parent to generate a peptide backbone
        # and a pseudo energy landscape approach to find an initial chain with the end point fixed
        # at point B.
        # create a regular backBone using PointsA as the first residue
        peptideBackbone = self.PBG.generateBuildingBlock(self.numResidues, seedResidue = [self.pointsA[2], self.pointsA[1], self.pointsA[0]])
        return peptideBackbone.getAtomsXYZ()  
        
    def PE(self, xyzVals):
        # Over load PE function to include two springs - between the last two points
        # ensures a smooth connection to other polymers
        PE = 0.0
        # add spring between pointB[2] and end of xyzVals with equilibrium at pointB
        dist1 = np.linalg.norm(xyzVals[-1] - self.pointsB[2])
        PE += 0.5 * 3 * self.springConstant * dist1**2
        
        # add spring between pointsB[1] and end of xyzVals[-2] 
        dist2 = 0.0
        # dist2 = np.linalg.norm(xyzVals[-2] - self.pointsB[1])
        # PE += 0.5 * 2 * self.springConstant * dist2**2
        
        return PE, np.sqrt(dist1**2 + dist2**2)
        
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
    pointB = np.array([0.0, 0.0, 0.0])
    
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
    numCrankMoves = 0
    fIO.saveXYZList([pointsA[0], pointsA[1], pointsA[2], pointsB[0], pointsB[1], pointsB[2]], ['Ca', 'Ca', 'Ca', 'O', 'O', 'O'], 'labEndPoints.xyz')

    envelopeList=["innersphere 4", "frustum 11.0 20.0 -5.0 2.0"]
    
    pointsToAvoid = backboneGenerator.generateBuildingBlock(10)
    pointsToAvoid.translateRefPointToTarget( (pointsA[2] + pointsB[2])/2.0 )
    newDirector = (pointsA[2] - pointsB[2] - np.array([3.0, 0.0, 3.0]))
    newDirectorHat = newDirector/np.linalg.norm(newDirector)
    pointsToAvoid.orientToDirector(newDirectorHat)
    fIO.saveXYZ(pointsToAvoid.blockXYZVals, 'Li', "pointsToAvoid.xyz")
    
    # build building block and dump to file
    hairpinBuildingBlock = hairPinGen.generateBuildingBlock(numResidues, pointsA, pointsB, minDist, numCrankMoves, pointsToAvoid = pointsToAvoid.blockXYZVals, envelopeList=envelopeList)
    hairpinBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "hairpin done"