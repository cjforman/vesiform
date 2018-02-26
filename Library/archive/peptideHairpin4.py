import sys
from Library.constrainedPolymer4 import ConstrainedPolymerPackBBG as CPBBG
from Builder.BuildingBlock import BuildingBlock
from scipy.optimize import minimize
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO
import numpy as np

class peptideHairpinGenerator(CPBBG):
    # This class uses the constrainedPolymer routine to generate
    # a sequence of particles that represent individual residues.  
    # The residues are then packed with an NCC triplet which connects
    # between the residue particles such that the peptide bond is planar and
    # the "stress" of the externally applied chain is accommodated
    # in the psi and phi dihedral angles. Where possible the 
    # the bond angles are correct.  The coordinates of the coarse residue correspond 
    # to the location of the N in the NCC triplet for 'NC' polarity and C for 'CN' 
    # polarity.
    
    def __init__(self, paramFilename):
        # initialise the parameter dictionary for the base classes
        CPBBG.__init__(self, paramFilename)    
    
    def initialiseParameters(self):
        # initialise the constrained polymer parent
        CPBBG.initialiseParameters(self)
        self.CCbondLength = self.getParam('CCbondLength')
        self.CNbondLength = self.getParam('CNbondLength')        
        self.phi = self.getParam('phi') * np.pi / 180.0
        self.psi = self.getParam('psi') * np.pi / 180.0
        self.omega = self.getParam('omega') * np.pi / 180.0
        self.angleN = self.getParam('angleN') * np.pi / 180.0
        self.angleCA = self.getParam('angleCA') * np.pi / 180.0
        self.angleC = self.getParam('angleC') * np.pi / 180.0
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for hairpin"
            sys.exit()        

    def generateBuildingBlock(  self, 
                                numResidues, 
                                pointA,
                                pointB,
                                pointC,
                                pointD,
                                pointE, 
                                minDist, 
                                residueLength, 
                                polarity,
                                pointsToAvoid=[], 
                                envelope='None'):
        
        self.polarity = polarity
        self.numResidues = numResidues
        self.numPoints = self.numResidues * 3
        self.pointC = pointC
        self.pointD = pointD
        self.pointE = pointE

        # first generate the residue chain 
        self.residueChainBB = CPBBG.generateBuildingBlock(self, numResidues, pointA, pointB, minDist, residueLength, pointsToAvoid=pointsToAvoid, envelope=envelope)
        self.residueChainBB.exportBBK("residueList.xyz")
        return self.generateAtomicBuildingBlock() 
        
    def generateAtomicBuildingBlock(self):
        xyzVals = self.convertResidueHairpinToAtomicHairpin(self.residueChainBB.blockXYZVals, self.polarity)
        names = self.generateAtomNames() 
        xyzValsLen = len(xyzVals)
        
        # always put the N terminal connector first for peptide building blocks
        if self.polarity == 'NC':
            connectors = [ [2, 1, 0], [xyzValsLen-3, xyzValsLen-2, xyzValsLen-1] ]
        else:
            connectors = [ [xyzValsLen-3, xyzValsLen-2, xyzValsLen-1], [2, 1, 0] ]
        refPoint = self.residueChainBB.blockRefPoint
        directorHat = self.residueChainBB.blockDirectorHat

        return BuildingBlock(xyzVals, names, connectors, refPoint, directorHat)  

    def convertResidueHairpinToAtomicHairpin(self, residueXYZList, polarity):
        
        # set the required bond angles and lengths based on polarity
        alpha0 = self.omega # (this is technically -180 or 180 based on polarity but its the same angle)
        if polarity=='NC':
            b1 = self.CNbondLength
            b2 = self.CCbondLength
            b3 = self.CNbondLength
            beta1 = self.angleN
            beta2 = self.angleCA
            beta3 = self.angleC
        else:
            # polarity=='CN':
            b1 = self.CCbondLength
            b2 = self.CNbondLength
            b3 = self.CNbondLength
            beta1 = self.angleC
            beta2 = self.angleCA
            beta3 = self.angleN

        # set up output array
        xyzVals = []

        # loop through the residue chain giving names to the surrounding positions        
        for curIndex, residueXYZ in enumerate(residueXYZList):
            # figure out the named positions surrounding the current residue
            # that will determine the psi and phi values necessary.- 
            # see page 201 in notebook 2 for diagram.
            # points C, D and E must be supplied externally.
            if curIndex==0:
                    S0 = self.pointC
                    S1 = self.pointD
                    S2 = residueXYZ
                    S3 = residueXYZList[curIndex + 1]
            elif curIndex==1:
                    S0 = self.pointD
                    S1 = xyzVals[-1]
                    S2 = residueXYZ
                    S3 = residueXYZList[curIndex + 1]
            elif curIndex==len(residueXYZList) - 1:
                S0 = xyzVals[-2]
                S1 = xyzVals[-1]
                S2 = residueXYZ
                S3 = self.pointE
            else:
                S0 = xyzVals[-2]
                S1 = xyzVals[-1]
                S2 = residueXYZ
                S3 = residueXYZList[curIndex + 1]

            # calculate the positions of m0 and m1                
            M0, M1 = self.residueToAtoms(S0, S1, S2, S3, alpha0, b1, beta1, b2, beta2, b3, beta3)
            
            # add the three relavent atom positions to the output array
            xyzVals.append(S2)
            xyzVals.append(M0)
            xyzVals.append(M1)
        
        return xyzVals

    def residueToAtoms(self, S0, S1, S2, S3, alpha0, b1, beta1, b2, beta2, b3, beta3):
        # Positioning M0 correctly in the S0, S1, S2 TNB frame
        # allows for alpha0, b1 and beta1 to be set correctly.
        TNB_S0_S1_S2 = coords.constructTNBFrame(S0, S1, S2)
        M0 = S2 + b1 * coords.generateTNBVecXYZ(TNB_S0_S1_S2, beta1, alpha0)

        TNB_S1_S2_M0 = coords.constructTNBFrame(S1, S2, M0)        
        # Can position M1 such that b2, and beta2 are correct, but that leaves b3 and beta3 and only one
        # parameter - alpha2 - which is dihedral about s2 to m0 axis. So find the alpha which
        # minimises the difference between b3 and linalg.norm(M1 - S3).
        results = minimize(coords.computeDistBetweenM1AndS3, [-60 *np.pi/180], (TNB_S1_S2_M0, beta2, S3, b3))
        alpha2 = results['x']

        # compute position of M1 with b2, beta2 and alpha2. We don't care what the last angle is. 
        M1 = M0 + b2 * coords.generateTNBVecXYZ(TNB_S1_S2_M0, beta2, alpha2)
        
        return M0, M1

    def generateAtomNames(self):
        if (self.polarity == "NC"):
            names = ['N', 'C', 'C'] * self.numResidues
        if (self.polarity == "CN"):
            names = ['C', 'C', 'N'] * self.numResidues
        return names

        
if __name__ == "__main__":
    
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the backbone generator object.
    hairPinGen = peptideHairpinGenerator(filename)

    # generate a backbone
    numResidues = 5
    bondLength = 3.5
    pointA = np.array([3, 0, 0])
    pointB = np.array([-3, 0, 0])
    pointC = pointA + bondLength * np.array([0.77, 0.77, 0.0])
    pointD = pointC + bondLength * np.array([0.0, 1.0, 0.0])
    pointE = pointD + bondLength * np.array([0.0, 1.0, 1.0])
    minDist = 3
    polarity = 'NC'
    fIO.saveXYZ([pointA, pointB, pointC, pointD, pointE], 'Ca', 'externalPoint.xyz')

    # build building block and dump to file
    hairpinBuildingBlock = hairPinGen.generateBuildingBlock(numResidues, pointA, pointB, pointC, pointD, pointE, minDist, bondLength, polarity)
    hairpinBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    print "hairpin done"