#!/usr/bin/env python
import os
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.peptideHairpin import peptideHairpinGenerator as PHPBBG
from PDBProc.pdbLib import PDB
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
import numpy as np
import random as rnd
import copy as cp

def CreatePDBFromSequenceString(sequenceFile, backboneFile, numCrankRandomizations=0, boxSize=None, forcefield='ff19SB', outputFilename=None, Term=True, minDist=20, postProcess=True, warp=False, **kwds): 
    print("Creating PDB from sequence file: ", sequenceFile)
    print("peptide parameters in file: ", backboneFile)
    if Term:
        print("Adding ACE and NME terminals")
    
    print("number of random cranks: ", numCrankRandomizations)
    print("warp: ", warp)

    # create a class for making PDBS
    PDBMaker = PDB()
    
    if not outputFilename:
        # generate the output filename
        outputFilename = sequenceFile[0:-4] + '.pdb'
    
    # generate a list of three letter sequence codes from the sequence file information
    seq = [ PDBMaker.ConvertAACodes(let) for let in fIO.readTextFile(sequenceFile)[0] ]
        

    # how many residues?
    N = len(seq)

    # increment number of residues by 2 if adding terminating residues ACE and NME.
    # also add the residue names to the list. 
    if Term:
        seq.insert(0, 'ACE')
        seq.append('NME')
        N = N + 2

    
    if boxSize:
        # box size should be 
        print("standard box size for system with ", N, " residues: ", np.power(N * 3 * 10, 1.0/3.0))
        print("limiting to specified box size: ", boxSize)    

    # Create the polymer building object with the parameters in the vesiform file 
    polymerBBG = PBG(backboneFile) # generate a polymer with the angles specified in the file.  

    if boxSize:        
        # create a hairpin generator object
        hairPinGen = PHPBBG(backboneFile)
        
        # generate 1 residue seed building block using the polymerBBG
        seedResidue = polymerBBG.generateBuildingBlock(1)
    
        # pick two points randomly well inside a box of boxSize centered at the origin. These are the start and end points of the peptide
        pointA = 0.9 * boxSize * np.array([ rnd.uniform(-.5, 0.5), rnd.uniform(-.5, 0.5), rnd.uniform(-.5, 0.5)])
        pnt = coords.pickRandomPointOnUnitSphere()
        pointADir = coords.polarToUnitSphereXYZ( pnt[0], pnt[1] )
        pointB = 0.9 * boxSize * np.array([ rnd.uniform(-.5, 0.5), rnd.uniform(-.5, 0.5), rnd.uniform(-.5, 0.5)])
        pnt = coords.pickRandomPointOnUnitSphere()
        pointBDir = coords.polarToUnitSphereXYZ( pnt[0], pnt[1] )
        
        # use the seedResidue to create the residue at a particular point and orientation.
        seedResidue.placeAtom(2, pointA)
        seedResidue.setBlockRefPoint(pointA)
        seedResidue.orientToDirector(pointADir)
        pointsA = cp.copy(seedResidue.blockXYZVals[:])
    
        # repeat with same seed residue but at different point and orientation.
        seedResidue.placeAtom(2, pointB)
        seedResidue.setBlockRefPoint(pointB)
        seedResidue.orientToDirector(pointBDir)
        pointsB = cp.copy(seedResidue.blockXYZVals[:])
    
        fIO.saveXYZList([pointsA[0], pointsA[1], pointsA[2], pointsB[0], pointsB[1], pointsB[2]], ['Ca', 'Ca', 'Ca', 'O', 'O', 'O'], 'labEndPoints.xyz')
        strBSize = str(boxSize/2.0)
        envelopeList=["cuboid -" + strBSize + " " + strBSize + " -" + strBSize + " " + strBSize + " -" + strBSize + " " + strBSize]
    
        # build building block 
        polymerBB = hairPinGen.generateBuildingBlock(N, pointsA, pointsB, minDist, numCrankRandomizations, envelopeList=envelopeList, warp=warp, **kwds)
        
        # extract the xyz vals as list
        xyzVals = polymerBB.blockXYZVals
    else:
        polymerBB = polymerBBG.generateBuildingBlock(N, showBlockDirector=False, nameCA=True, warp=warp, **kwds)
        
        # extract XYZ vals as a list and move the end point to the origin so all z vals are +ve 
        # avoids a formatting issue with large -ve numbers in the PDB file of a very long straight chain 
        xyzVals = [ r - polymerBB.blockXYZVals[-1] for r in polymerBB.blockXYZVals]
        
    # reformat the xyzVals list into a useful size 
    xyzVals = np.reshape(xyzVals, (N, 3, 3) )
    
    # set up output array.    
    atoms = []
    atomIndex = 1
    resIndex = 1

    # generate PDB atom information strings
    for resName, posns in zip(seq, xyzVals):

        if resName=='ACE':
            # only append the C of the coords 
            atoms.append([atomIndex, 'C', '', 'ACE', 'A', resIndex, '', posns[2][0], posns[2][1], posns[2][2], 0.0, 0.0, '', 'C', ''] )
            atomIndex += 1
            
        elif resName=='NME':
            # only append the N of the coords 
            atoms.append([atomIndex, 'N', '', 'NME', 'A', resIndex, '', posns[0][0], posns[0][1], posns[0][2], 0.0, 0.0, '', 'N', ''] )
            atomIndex += 1
            
        else:
            # create atoms in blocks of three for each residue.
            atoms.append([atomIndex, 'N', '', resName, 'A', resIndex, '', posns[0][0], posns[0][1], posns[0][2], 0.0, 0.0, '', 'N', ''] )
            atomIndex += 1
            
            atoms.append([atomIndex, 'CA', '', resName, 'A', resIndex, '', posns[1][0], posns[1][1], posns[1][2], 0.0, 0.0, '', 'C', ''] )
            atomIndex += 1
            
            atoms.append([atomIndex, 'C', '', resName, 'A', resIndex, '', posns[2][0], posns[2][1], posns[2][2], 0.0, 0.0, '', 'C', ''] )
            atomIndex += 1
            
        # increment residue number
        resIndex += 1
    
    print("Writing outputfile: ", outputFilename)
        
    # output the atoms list as a PDB file
    PDBMaker.writePDBAtomsFromList(atoms, outputFilename)
    
    if postProcess:
        # post process PDB
        l1 = "reduce -trim " + outputFilename + " > " + outputFilename[0:-4] + "_noh.pdb"
        os.system(l1)


        L2 = ["source leaprc.protein." + forcefield + "\n", 
              "mol = loadpdb " + outputFilename + '\n',
              "savepdb mol " + outputFilename[0:-4] + 't.pdb \n',
              "saveamberparm mol " + outputFilename[0:-4] + '.prmtop ' + outputFilename[0:-4] + ".inpcrd \n",
              "quit \n"]
        
        fIO.writeTextFile(L2, "tleap_source")
        os.system("tleap -f tleap_source") 

if __name__=="__main__":
    import sys
    import json

    try:      
        with open( sys.argv[2] ) as f:
            params = json.load(f)
    except IndexError as e:
        print("Usage: \n", 
              "CreatePDBFromSequenceString.py seqFile config.json\n\n",
              "Params in json are:\n",
              "numCrankRandomizations=0\n",
              "boxSize=None\n", 
              "outputFilename\n", 
              "Term= (True or False, 1, or 0)\n", 
              "minDist=20\n", 
              "postProcess=True\n", 
              "warp=False\n")
       
    CSQ = CreatePDBFromSequenceString(sys.argv[1], params['oldSchoolFilename'], **params)    
    
    print("done")