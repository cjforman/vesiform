import sys
from Utilities.keyProc import keyProc
import Utilities.cartesian as cart
import Utilities.coordSystems as coords 
import Utilities.fileIO as fIO
import numpy as np
import random as rnd
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG

class peptideBackboneGenerator(BBG):
    # This is a class for generating a building block with a peptide backbone
    # The input filename contains the backbone parameters, e.g betastrand and alphahelix.
    # standard parameters for these are saved in the library directory.
    # The strand is always from N to C terminus. If you want it C to N terminus then just
    # rotate it. 
    
    def __init__(self, paramFilename):
        # initialise the parameter dictionary
        keyProc.__init__(self, paramFilename)
    
    def initialiseParameters(self):

        keyProc.initialiseParameters(self)
    
        try:
            self.CCbondLength = self.getParam('CCbondLength')
            self.CNbondLength = self.getParam('CNbondLength')        
            self.phi = self.getParam('phi') * np.pi / 180.0
            self.psi = self.getParam('psi') * np.pi / 180.0
            self.omega = self.getParam('omega') * np.pi / 180.0
            self.angleN = self.getParam('angleN') * np.pi / 180.0
            self.angleCA = self.getParam('angleCA') * np.pi / 180.0
            self.angleC = self.getParam('angleC') * np.pi / 180.0
            self.dumpInterimFiles = self.getParam('dumpInterimFiles')
        except KeyError as e:
            print(e)
            sys.exit()
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for peptide BackboneGenerator")
            sys.exit()        

    def generateBuildingBlock(self, numResidues, minDist=1.0, seedResidue=None, showBlockDirector=False, nameCA=False, warp=False, **kwds):
    
        self.numResidues = numResidues    
        self.numPoints = numResidues * 3 
        self.directorHat = np.array([0.0, 0.0, 1.0])
        self.nameCA = nameCA
        self.seedResidue = seedResidue
        if not seedResidue==None:
            if not len(seedResidue)==3:
                print("Warning: Seed residue in peptide backbone generator is not 3 atoms long") 
        self.minDist = minDist
        self.kwds = kwds
        
        return BBG.generateBuildingBlock(self, self.numPoints, minDist, showBlockDirector=showBlockDirector, warp=warp, **kwds)

    def generateBuildingBlockXYZ(self):

        # create the first residue in the BB
        if self.seedResidue==None:
            strand = self.startResidue()
        else:
            strand = self.seedResidue

        # tally the number of residues we added to strand        
        numResidues = 1

        checkClashes=False
        try:
            checkClashes = bool(self.kwds['checkClashes'])
        except KeyError:
            pass

        try:
            if checkClashes and bool(self.kwds['warp']):
                print("Warning:  You have selected check clashes and warp. The input into the warp function is checked for clashes, but not the output. You have been warned. We are not responsible for chain crossings caused by reckless warping.")
        except KeyError:
            pass


        consecutiveClashes = 0

        # loop until we've added the right number of residues        
        while numResidues < self.numResidues:
            
            # create new residue from the last three points that obeys bond angles and dihedral sampling rules etc
            residue = self.generateResidue(strand[-3:], numResidues)
            
            # check to see if we are checking for clashes
            if checkClashes:
                # only add residue if it doesn't come too close to another residue
                if self.checkResidue(residue, strand, self.minDist): 
                    # add the new Residues to the strandVec list
                    [strand.append(atom) for atom in residue]

                    numResidues += 1
                    
                    # output update. 
                    print("current length: ", numResidues, " out of ", self.numResidues)
                    # set the consecutive Clashes Counter to zero
                    consecutiveClashes=0
                else:
                    consecutiveClashes += 1
                    if consecutiveClashes > self.kwds['maxClashes']:
                        print("Max consecutive clashes reached. Removing ", self.kwds['numResiduesToPrune'], " residues (3 atoms each) from strand.")
                        consecutiveClashes = 0
                        print("original Strand Length: ", len(strand))
                        strand = strand[0:-1 * (3 * int(self.kwds['numResiduesToPrune']))]
                        print("new Strand Length: ", len(strand))
                        if len(strand)%3==0:
                            # reset the number of residues count
                            numResidues = int(len(strand)/3)
                        else:
                            print("Warning: Strand Length is no longer multiple of 3. Somethings gone wrong.")
                        print("new residue Length: ", numResidues)
                        # if we accidentally prune the entire strand restart it
                        if len(strand)<3:
                            # create the first residue in the BB
                            if self.seedResidue==None:
                                strand = self.startResidue()
                            else:
                                strand = self.seedResidue

            else:
                # add residue without checking
                [strand.append(atom) for atom in residue]
                numResidues += 1
            
        
        # perform final orientation depending on whether or not a seed residue was provided.
        # If seed was provided then don't fiddle with orientation at all.
        # if seed wasn't provided then return the strand with it's centre of mass
        # at the origin and it's helical axis pointing along the z-axis
        if self.seedResidue == None and self.numResidues > 3:
            # This structure was created in such a way that the first two atoms are
            # aligned with Z vector. We want the structure to be returned so it is at it's
            # centre of mass and also with the helical axis pointing up the z axis.
            # We will use the transformation from Lab to Block to achieve this,
            # with a little bit of hacking of the member variables.
            
            # find the helix axis and the centre of mass of the construct
            # as it is right now.
            currentDirector = coords.axisFromHelix(strand)
            currentDirectorHat = currentDirector/np.linalg.norm(currentDirector) 
            currentCOM = cart.getCentreOfMass(strand)
            
            # Now specify the orientation and position that we want to map to. 
            # store the refPoint as a member variable.
            targetRefPoint = np.array([0.0, 0.0, 0.0])
            targetDirectorHat = np.array([0.0, 0.0, 1.0])
            
            # modify the strand accordingly with the universal transformation function
            strand = coords.transformFromBlockFrameToLabFrame(targetDirectorHat, 
                                                              targetRefPoint, 
                                                              0.0, 
                                                              currentDirectorHat, 
                                                              currentCOM, 
                                                              strand)
     
        return list(np.around(strand, 12)) # round the positions to 12.dp - cleans up machine precision noise on zeroes 


    def checkResidue(self, residue, strand, minDist):

        retVal = True
        # Compare minDist against the distance between each point in residue and each point in strand-except last three points of strand.
        # There are three points in residue, but strand is just a list of points.
        # If any points are < minDist return False
        if False in [ np.linalg.norm(resPos - strandPos) > minDist for strandPos in strand[0:-3] for resPos in residue]:
            retVal = False
        # compare non-adjacent points in residue and the last point on the strand. if they are less than min Dist then reject 
        elif False in [ np.linalg.norm(resPos - strandPos) > minDist for strandPos in strand[-1] for resPos in residue[1:]]:
            retVal = False
            
        return retVal
    
    def startResidue(self):
        # generates the first residue of the backbone
       
        # set first NPos to be at origin.
        NPos = np.array([0.0, 0.0, 0.0])

        # Assume CAPos is along director axis from NPos - bad assumption - we'll correct it later
        CAPos = NPos + self.CNbondLength * self.directorHat

        # generate a third dummy vector 
        DummyPos = NPos + np.array([1.0, 0.0, 0.0]) 
        
        # determine phi and psi
        phi, _ = self.pickPhiPsi(0)

        # compute the final CPos for the initial residue
        TNB3 = coords.constructTNBFrame(DummyPos, NPos, CAPos)
        CPos = CAPos + self.CCbondLength *  coords.generateTNBVecXYZ(TNB3, self.angleCA, phi)
    
        return [ NPos, CAPos, CPos]
                
    # can set a dictionary containing several ranges of allowed phi psi angles
    # we always pick randomly from one of these ranges. 
    # can limit to one very narrow range if you want a particular secondary structure  
    def pickPhiPsi(self, resIndex):
        
        # figure out which phiPsi Range is of interest
        try:
            if self.kwds['blockLength']==0:
                phiPsiRange = rnd.choice( list(self.kwds['phiPsiRanges'].values()) )
            else:
                cycleLength = len(self.kwds['phiPsiRanges']) * self.kwds['blockLength']
                
                curPhiPsiRangeIndex  = int(np.floor( (resIndex % cycleLength ) / self.kwds['blockLength'] ) )
                
                phiPsiRange = list(self.kwds['phiPsiRanges'].values())[curPhiPsiRangeIndex]
    
            # pick a value at random from the chosen phi psi range    
            psi = rnd.uniform( float( phiPsiRange['minPsi'])* np.pi / 180.0, float(phiPsiRange['maxPsi'] )* np.pi / 180.0 )
            phi = rnd.uniform( float( phiPsiRange['minPhi'])* np.pi / 180.0, float(phiPsiRange['maxPhi'] )* np.pi / 180.0 )
        except KeyError as e:
            print("Phi/Psi Range not specified in input json. Defaulting to .txt file input via keyProc system.")
            phi = self.phi
            psi = self.psi

        return phi, psi
    
    def generateResidue(self, prevRes, resIndex):
        
        # given a list of 3 xyz positions for the previous residue
        # construct the next set of three positions with appropriate bond angles
        # and dihedral angles.

        phi, psi = self.pickPhiPsi(resIndex)
                            
        # construct the first TNB from the previous residue
        TNB1 = coords.constructTNBFrame(prevRes[0], prevRes[1], prevRes[2])
        NPos = prevRes[2] + self.CNbondLength * coords.generateTNBVecXYZ(TNB1, self.angleC, psi)
        
        TNB2 = coords.constructTNBFrame(prevRes[-2], prevRes[-1], NPos)
        CAPos = NPos + self.CNbondLength * coords.generateTNBVecXYZ(TNB2, self.angleN, self.omega)

        TNB3 = coords.constructTNBFrame(prevRes[-1], NPos, CAPos)
        CPos = CAPos + self.CCbondLength *  coords.generateTNBVecXYZ(TNB3, self.angleCA, phi)
    
        return [NPos, CAPos, CPos]

    def generateBuildingBlockNames(self):
        if self.nameCA:
            names = ['N', 'CA', 'C'] * self.numResidues
        else:
            names = ['N', 'C', 'C'] * self.numResidues
        
        return names[0:self.numPoints]
    
    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)
    
    def generateBuildingBlockDirector(self):
        if self.numResidues > 3:
            director = coords.axisFromHelix(self.buildingBlockXYZ)
        else:
            director = self.buildingBlockXYZ[-1] - self.buildingBlockXYZ[0]
        directorHat = director/np.linalg.norm(director)   
        return directorHat 

    def generateBuildingBlockConnectors(self):
        return [np.array([2, 1, 0]), np.array([self.numPoints - 3, self.numPoints - 2, self.numPoints - 1])]
                
if __name__ == "__main__":
    
    
    # get the file name from the command line
    filename = "alphahelix.txt"

    # create the backbone generator object using static file parameters
    backboneObject = peptideBackboneGenerator(filename)

    # generate backbone realtime parameters
    numResidues = 16
    startPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    
    backBoneBuildingBlock = backboneObject.generateBuildingBlock(numResidues, showBlockDirector=False)
    backBoneBuildingBlock.transformBBToLabFrame(director, startPos, rotation)
    backBoneBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    filename = "betastrand.txt"

    # create the backbone generator object using static file parameters
    backboneObject = peptideBackboneGenerator(filename)

    # generate backbone realtime parameters
    numResidues = 15
    startPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    
    backBoneBuildingBlock = backboneObject.generateBuildingBlock(numResidues)
    backBoneBuildingBlock.transformBBToLabFrame(director, startPos, rotation)
    backBoneBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))


    
    print("backbone done")