import sys
from Utilities.keyProc import keyProc
import Utilities.cartesian as cart
import Utilities.fileIO as fIO
import numpy as np
from PDBProc.pdbLib import PDB as pdb
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG


class pdbPeptideBackboneGenerator(BBG):
    # This is a class for generating a building block with a peptide backbone
    # that is drawn from a PDB file. Also allows to draw the entire PDB, 
    # and save a set of modified coords to file in the same PDB format 
    # that it was opened in.
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
            print("Critical Parameters are undefined for pdbPeptideBackboneGenerator")
            sys.exit()        

    def generateBuildingBlock(self, backboneOnly=True, director=None, showBlockDirector=False, nameCA=False):
        self.backboneOnly = backboneOnly
        self.backboneIndices = self.extractBackBoneIndices()
        self.numPoints = len(self.backboneIndices)
        self.nameCA = nameCA
        self.director = director
        self.dumpInterimFiles = 0
        minDistDummy = 0.0
        return BBG.generateBuildingBlock(self, self.numPoints, minDistDummy, showBlockDirector=showBlockDirector)

    def generateBuildingBlockXYZ(self):
        if self.backboneOnly:
            retXYZ = [ np.array([self.pdb.atoms[i][7], self.pdb.atoms[i][8], self.pdb.atoms[i][9]]) for i in self.backboneIndices ]
        else:
            retXYZ = [ np.array([atom[7], atom[8], atom[9]]) for atom in self.pdb.atoms]
        return retXYZ   
        
    def generateBuildingBlockNames(self):
        if self.backboneOnly:
            # extract the atom names to match the backbone indices
            names = [ self.pdb.atoms[i][1] for i in self.backboneIndices ]
            
            # replace CA with C in the names array for blender - if requested (which it is by default)
            if not self.nameCA:
                names = [ name if not name =='CA' else 'C' for name in names]
        else:
            names = [ atom[1] for atom in self.pdb.atoms]
        return names
    
    def generateBuildingBlockDirector(self):
        # check to see if there was an externally imposed director
        director = None
        if None in self.director:
            # calculate the eigenvectors of the inertial tensor of the protein and 
            # return the unit vector along the principal axis ( assumes each point is of unit mass so not quite correct).
            pAxis = cart.getPrincipalAxis(self.buildingBlockXYZ)
            director = pAxis/np.linalg.norm(pAxis)
        else:
            director = self.director

        return director

    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)
    
    def extractBackBoneIndices(self):
        # Extracts the index numbers of the backbone atoms
        return  [ atomIndex for atomIndex, atom in enumerate(self.pdb.atoms) if atom[1] in ['N', 'CA', 'C']]
    
    def exportPDBWithNewCoords(self, newXYZVals, filename):
        self.pdb.replacePdbAtoms(newXYZVals, filename)
        
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
    
    print("pdb building block done")