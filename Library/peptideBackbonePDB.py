import sys
from Utilities.keyProc import keyProc
import Utilities.cartesian as cart
import Utilities.fileIO as fIO
import numpy as np
from PDBProc.pdbLib import PDB as pdb
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG


class pdbPeptideBackboneGenerator(BBG):
    # This is a class for generating a building block with a peptide backbone
    # that is drawn from a PDB file.
    #
    # The input filename is the pdb filename. This is parsed
    # for the backbone automatically at loading.
    
    def __init__(self, pdbFilename):
        # This is a keyproc object because a BuildingBlockGenerator is 
        # a keyproc object. However, object has no parameters so no need 
        # to invoke the keyproc part of the object
        # keyProc.__init__(self, paramFilename)
        
        # create a pdb object which automatically loads and parses
        # the pdb file.
        self.pdb = pdb(pdbFilename)
        self.numPoints = len(self.pdb.atoms)
    
    def initialiseParameters(self):
        
        keyProc.initialiseParameters(self)
    
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for pdbPeptideBackboneGenerator"
            sys.exit()        

    def generateBuildingBlock(self, director=None, showBlockDirector=False, nameCA=False):
        self.backboneIndices = self.extractBackBoneIndices()
        self.numPoints = len(self.backboneIndices)
        self.nameCA = nameCA
        self.director = director
        self.dumpInterimFiles = 0
        minDistDummy = 0.0
        return BBG.generateBuildingBlock(self, self.numPoints, minDistDummy, showBlockDirector=showBlockDirector)

    def generateBuildingBlockXYZ(self):
        return  [ np.array([self.pdb.atoms[i][7], self.pdb.atoms[i][8], self.pdb.atoms[i][9]]) for i in self.backboneIndices ] 
        
    def generateBuildingBlockNames(self):
        # extract the atom names to match the backbone indices
        names = [ self.pdb.atoms[i][1] for i in self.backboneIndices ]
        
        # replace CA with C in the names array for blender - if requested (which it is by default)
        if not self.nameCA:
            names = [ name if not name =='CA' else 'C' for name in names]
        return names
    
    def generateBuildingBlockDirector(self):
        # check to see if there was an externally imposed director
        director = None
        if not self.director==None:
            director = self.director
        else:
            # calculate the eigenvectors of the inertial tensor of the protein backbone and 
            # return the unit vector along the principal axis ( assumes each point is of unit mass so not quite correct).
            pAxis = cart.getPrincipalAxis(self.buildingBlockXYZ)
            director = pAxis/np.linalg.norm(pAxis)
        return director

    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)
    
    def extractBackBoneIndices(self):
        # Extracts the index numbers of the backbone atoms
        return  [ atomIndex for atomIndex, atom in enumerate(self.pdb.atoms) if atom[1] in ['N', 'CA', 'C']]

        
if __name__ == "__main__":
    
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the pdb object generator
    pdbObject = pdbPeptideBackboneGenerator(filename)

    # generate backbone realtime parameters
    refPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    
    pdbBuildingBlock = pdbObject.generateBuildingBlock(showBlockDirector=True) 
    pdbBuildingBlock.transformBBToLabFrame(director, refPos, rotation)
    pdbBuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'pdb')) # adds the xyz extension automatically
    
    print "pdb building block done"