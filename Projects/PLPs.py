'''
Created on Aug 17, 2021

@author: Chris

Project classes for the Tau/PLP RO1 proposal.

'''
import os, sys
import json
from Projects.GeneralPolymerBrush import GeneralBlockPolymer as GBP
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
from PDBProc.pdbLib import PDB as PDBL
from PDBProc.Aligner import Aligner as Al
import numpy as np
import copy as cp

class AllthePLPS():
    def __init__(self, filename):
        # load the filename
        try:
            with open( filename) as f:
                plpStructs = json.load(f)        
        except Exception as e:
            print("Unable to load filename:", filename, e)
            sys.exit()
        # check the PLP Structs for errors before we do them for real
        self.cycleThroughPLPs(plpStructs, dummyRun=True)
        self.cycleThroughPLPs(plpStructs, dummyRun=False)
    
    def cycleThroughPLPs(self, plpStructs, dummyRun=False):
        # capture current working dir
        pwd = os.getcwd()        

        # loop through each PLP
        for plpStruct in plpStructs:
            if plpStruct not in ["choose", "oneToDo"]:
                if (plpStructs["choose"]=="doAll") or (plpStructs["choose"]=="doOne" and plpStructs["oneToDo"]==plpStruct):
                    print("*****" +  plpStruct + "* starting ****")
                    # make a directory and change to that directory
                    os.system("mkdir " + str(plpStruct))
                    os.chdir(str(plpStruct))
                    os.system("cp ../data/*.txt .")
                    os.system("cp ../data/*.pdb .")
                    
                    # save the dictionary for the current PLP to file as a json
                    self.saveObject(plpStruct + '.json', plpStructs[plpStruct])
            
                    # run the plp creator from the working dir with the structure
                    # if dummyRun is true then all it does it check the PLP structure for errors. 
                    GenericPLP_Sequence(plpStruct + '.json', dummyRun=dummyRun)
            
                    # change back up to the root directory
                    os.chdir(pwd)
                    
                    if dummyRun:
                        print("*****" +  plpStruct + " preRun Check complete. No Errors Detected *****")
                    else:
                        print("*****" +  plpStruct + " complete *****")
                else:
                    print("Not building PLP: ", plpStruct)
    

            

            
    
    def saveObject(self, configfile, objectId):
        """
        Saves the configuration data for the given object

        Args:
            configfile (File): The configuration filename

        """
        with open( configfile, "w" ) as f:
            json.dump( objectId, f, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__) 
            
    
    

class GenericPLP_Sequence(PDBL):
    def __init__(self, filename, dummyRun=False):
    
        # initialize the PDBObj parent class without a PDB input file 
        super().__init__( )

        # set up the plp structure dict    
        self.plpStruct = self.setUpPlpStruct(filename)

        try:
            chainLet =  self.plpStruct['chainLet']
        except KeyError as e:
            print('Chain Letter for PLP not specified. Defaulting to A', e)
            chainLet = 'A'

        if dummyRun==False:

            # create an aligner object for aligning vectors. 
            self.aligner = Al()
    
            # set up the stapler dictionary from the stapler file
            stapler = self.SetupStapler()
            
            # create a list of brush dictionaries
            brushList = self.GenerateBrushes()
            
            # Use a stapler to set up the sparse vesform polymers into roughly the right set up
            self.PolymerisePLPStapler(brushList, stapler)
            
            # write out the PDB to file with each brush as a separate chain.
            self.writePDB(brushList, self.plpStruct['name'] + '.pdb', terminateBrushes=True)
            
            # generate a leap script which populates each sparse vesiform shape with the right extra bits
            # and then runs sander to minimise each structure by itself in roughly the right place but not connected to backbone.
            tleapPDB, minPDB, bashscriptname = self.leapItUp(self.plpStruct['name'])
            
            os.system("./" + bashscriptname)
            # execute the process via wsl 
            # print(subprocess.run(['wsl',  "./" + bashscriptname], env={'PATH':'d:/box/Projects/amber/amber20/bin', 'AMBERHOME':'d:/box/Projects/amber/amber20'}))
    
            # load in the minimised brushes into a new brush list
            self.generateBrushListFromPDB(brushList, minPDB)
            
            # Repeat the "polymerise" step with energy relaxed brushes - no need for a stapler this time.
            # this time perform the final phase and azimuth rotations  
            self.PolymerisePLPRotate(brushList)
            
            # write out the PDB to file
            self.writePDB(brushList, self.plpStruct['name'] + '_final' + '.pdb', chainLet=chainLet)
   

    def setUpPlpStruct(self, filename):
        errVal = False
        # load the plp structure into the object
        try:
            with open( filename) as f:
                plpStruct = json.load(f)        
        except Exception as e:
            print("Unable to load filename:", filename, e)
            errVal = True
            sys.exit()

        # check the PLP structure for validity
        try:
            print("Structure Name: ", plpStruct['name'])
        except KeyError as e:
            print("Crucial Key Missing from PLP structure: ", e)
            errVal = True

        try:
            print("backboneDihedralMin: ", plpStruct['backBoneDihedralMin'])
        except KeyError as e:
            print("Crucial Key Missing from PLP structure: ", e)
            errVal = True                
        
        try:
            print("backboneDihedralMax: ", plpStruct['backBoneDihedralMax'])
        except KeyError as e:
            print("Crucial Key Missing from PLP structure: ", e)
            errVal = True
            
        try:
            print("backBoneBondAngleMin: ", plpStruct['backBoneBondAngleMin'])
        except KeyError as e:
            print("Crucial Key Missing from PLP structure: ", e)
            errVal = True
        
        try:
            print("backBoneBondAngleMax: ", plpStruct['backBoneBondAngleMax'])
        except KeyError as e:
            print("Crucial Key Missing from PLP structure: ", e)
            errVal = True                
        
        try:
            if plpStruct['stapler']=='LEU.pdb':
                stapleResName = 'L'
            elif plpStruct['stapler']=='VAL.pdb':
                stapleResName = 'V'
            else:
                print("Invalid StapleRes File Selected. Generate LEU.pdb or VAL.pdb using tleap. OR add your new file to the error check list in PLPs.py")
                errVal = True
        except KeyError as e:
            print("Crucial Key Missing from PLP structure: ", e)
            errVal = True
            
        try:
            for block in plpStruct['blocks']:
                try:
                    for i, seq in enumerate(plpStruct['blocks'][block]['sequences']):
                        if stapleResName not in seq:
                            print("Staple Res Missing from sequence: " + str(i + 1) + " in block " + block)
                            errVal = True
                        if seq[int(plpStruct['blocks'][block]['stapleRes'][i]) - 1]!=stapleResName:
                            print("Amino acid indicated by StapleRes entry is not a staple residue: sequence entry " + str(i + 1) + " in " + str(block))
                            errVal = True
                        if plpStruct['blocks'][block]['secondary'][i] not in ['coil', 'alphahelix', 'betastrand', 'file']:
                            print("Must specify 'coil', 'alphahelix' or 'betastrand' for secondary structure in sequence " + str(i) + " in " + str(block))
                            errVal = True
                        if plpStruct['blocks'][block]['secondary'][i]=='file':
                            try:
                                plpStruct['blocks'][block]['secondaryFile']
                            except KeyError as e:
                                print("Secondary Structure defined as filename, but filename not specified. Add attribute ", e)
                                errVal=True                                    
                except IndexError as e:
                    print("Input array length inconsistency detected in " + str(block) + ". ", e)
                    errVal=True
                if len(plpStruct['blocks'][block]['sequences'])!=len(plpStruct['blocks'][block]['secondary']):
                    print("Inconsistent number of brush sequences and secondary structure definitions " + str(block))
                    errVal = True
                if len(plpStruct['blocks'][block]['sequences'])!=len(plpStruct['blocks'][block]['azimuth']):
                    print("Inconsistent number of brush sequences and azimuth definitions " + str(block))
                    errVal = True
                if len(plpStruct['blocks'][block]['sequences'])!=len(plpStruct['blocks'][block]['phase']):
                    print("Inconsistent number of brush sequences and phase definitions " + str(block))
                    errVal = True
                print(str(block) + " will be repeated " + str(plpStruct['blocks'][block]['repeat']) + " times.")
        except KeyError as e:
            print("Crucial Key Missing from PLP structure:", e)
            errVal = True

        if errVal:
            print("Errors in PLP Structure definition detected.")
            sys.exit()
        
        # if we survived that lot then PLP struct is good. 
        return plpStruct

    def SetupStapler(self):
        # Set up the staple
        stapleAtoms = self.readPDBAtomsRet(self.plpStruct['stapler']) # read in the staple atoms
        staplerXYZ = [ np.array([atom[7], atom[8], atom[9]]) for atom in stapleAtoms ] # get the xyz points
        staplerNames = [ atom[1] for atom in stapleAtoms ] # get the names of the atoms

        # Set up the sets of stapler names that are useful and depend on which stapler we are using.
        if self.plpStruct['stapler']=='LEU.pdb':
            staplerNamesToInsert = ['CD2','CD1','CG','CB'] # stapler atoms to copy into brush array in reverse order of appearance
            staplerBBJoin = ['CD1', 'CG', 'CD2'] # the order of these names is crucial as it sets up the pairwise connectivity with the backbone
        elif self.plpStruct['stapler']=='VAL.pdb':
            staplerNamesToInsert = ['CG2','CG1','CB'] # must be inserted last first so write in reverse order they should appear in PDB
            staplerBBJoin = ['CG1', 'CB', 'CG2'] # order crucial
            
        # set up the indices - the order is important here or it tries to align the wrong atoms with the backbone carbons
        staplerBackBoneJoinIndices = [ staplerNames.index(staplerBBJoin[0]), 
                                       staplerNames.index(staplerBBJoin[1]),
                                       staplerNames.index(staplerBBJoin[2])]

        # N, CA and C are in the right order already in the PDB ordering so this works.
        staplerBrushJoinIndices = [ i for i, atom in enumerate(stapleAtoms) if atom[1] in ['N', 'CA', 'C'] ]
        
        # generate the symbols and tweak some of them to a different symbol to aid debugging output with blender visualisation
        staplerSymbols = [ self.lookupSym(atom[1]) for atom in stapleAtoms ]

        return {'atoms': stapleAtoms, 
                'staplerXYZ':staplerXYZ, 
                'staplerNames': staplerNames, 
                'staplerBackBoneJoinIndices': staplerBackBoneJoinIndices,
                'staplerBrushJoinIndices': staplerBrushJoinIndices, 
                'staplerSymbols': staplerSymbols,
                'staplerNamesToInsert': staplerNamesToInsert}


    def GenerateBrushes(self):
        # create a list to store the brush dictionaries
        brushStructList = []

        # set up a counter to count all the residues in previous brushes
        prevResidues = 0
        
        # loop through each block in the PLP        
        for block in self.plpStruct['blocks']:
            
            # repeat as specified
            for _ in range(self.plpStruct['blocks'][block]['repeat']):
                
                # for each repeating sequence pattern defined in this block generate a polymer N-CA-C backbone with the right attributes 
                for seq, secondary, Az, Phase, stapleRes in zip(self.plpStruct['blocks'][block]['sequences'], 
                                                                self.plpStruct['blocks'][block]['secondary'],
                                                                self.plpStruct['blocks'][block]['azimuth'],
                                                                self.plpStruct['blocks'][block]['phase'],
                                                                self.plpStruct['blocks'][block]['stapleRes']):
                    
                    # get the number of N-CA-C triplets in the backbone (num of residues)
                    numTriples = len(seq)
                    
                    # get the names of the backbone
                    names = numTriples * ['N', 'CA', 'C']
                    
                    # get the symbol of each atom
                    symbols = [ self.lookupSym(name) for name in names ]
                    
                    # set the residue number for each atom in the backbone (resid is a 1 based index)
                    resIds = [ [i, i, i] for i in range(1, numTriples + 1) ]
                    resIds = [ resNum + prevResidues for resId in resIds for resNum in resId ] # flatten the residue number array

                    # compute the joining indices for this brush backbone relative to the first back bone atom in the brush (N CA, C in the stapleRes'th residue
                    joiningIndices = [ (stapleRes - 1) * 3, (stapleRes - 1) * 3 + 1, (stapleRes - 1) * 3 + 2 ]  
                    
                    # set the residuenames for each atoms in the brush
                    resNames = [ [resName, resName, resName] for resName in seq]
                    resNames = [ resName for resNameT in resNames for resName in resNameT ] # flatten res names array
                    resNames = [ self.ConvertAACodes(res) for res in resNames ] # convert to the three letter names
                    
                    # generate the starting coords for each brush
                    if secondary=='coil':
                        xyz = self.makeRandomCoil( numTriples * 3 ) 
                    
                    elif secondary=='alphahelix':
                        xyz = self.makePeptideBackBone( numTriples, filename='alphahelix.txt' )
                    
                    elif secondary=='betastrand':
                        xyz = self.makePeptideBackBone( numTriples, filename='betastrand.txt' )
                        
                    elif secondary=='file':
                        secPDB = PDBL( self.plpStruct['blocks'][block]['secondaryFile'] ) 
                        xyz = secPDB.extractCoords()
            
                    # create a dictionary to hold each brush in
                    brush = {'names': names, 
                             'xyz': xyz, 
                             'symbols': symbols,
                             'resIds': resIds, 
                             'resNames': resNames, 
                             'sequence': seq, 
                             'type': secondary, 
                             'Az':  Az,
                             'Phase': Phase,
                             'joinIndices': joiningIndices,
                             'stapleRes': stapleRes}
            
                    # add the brush to the list of brushes 
                    brushStructList.append(brush) 
            
                    # increment the previous residue count by the number of residues we added to the dictionary 
                    prevResidues += numTriples
            
        return brushStructList


    def PolymerisePLPStapler(self, brushList, stapler):
        ''' This function creates a backbone polymer, and positions each Brush on it.
        
        The backbone polymer we create is just a dummy template.  That atoms that make the actual
        backbone are supplied by the brushes themselves - the side chain on the polymerising amino acid.   
        
        We are using valine or leucine for the purposes of the cartoon to side step having to create
        a new amino acid unit in tleap which is a mission. The valine and leucine's have a nice forked
        carbon on the end of the side chain and they can offer up these two carbons to form the back bone.  
         
        The polymer back bone for a valine or leucine linked system is a carbon chain that has 3x the number of 
        carbons as there are brushes. Each successive carbon triad in the polymer chain is taken to be 
        the CG1 , CB and CG2 of the valine linker unit or the CD1, CG and CD2 of the leucine. 
        
        We therefore create a pdb containing a single valine or leucine which acts as the stapler. 
        The function  'SetupStapler' has created a dictionary called stapler which has named all the attributes
        we need of the stapler in the dictionary, whether it is leucine or valine. That structure is passed 
        in as a variable along with a list of all the brushes we are working with. 
        
        First we position the staple so it's forked carbons aligns with the current triad of carbon atoms
        selected in the polymer backbone template.
        
        Then we position the brush so that the N, CA, C backbone in the linking valine or leucine, aligns with 
        the N, CA, C of the backbone aligned valine.
        
        Then we add the carbons of the stapler to the carbons of the brush.  This means Tleap 
        does not add the missing atoms for stapleRes residue in the brush in the wrong place. 
        
        We do this for each brush - including the little ones which are spacers between the big ones.
        
        The updated brush structures are modified in place in the brushList for the output.'''
        numBrushes = len(brushList)

        # generate a chain of points which is 3 times the number of brushes a carbon bond length apart.
        backBoneXYZ  = self.makeRandomCoil(3 * numBrushes, 
                                           'RandomPolymer.txt',
                                           bondLength=1.4, 
                                           alpha1=self.plpStruct['backBoneDihedralMin'], 
                                           alpha2=self.plpStruct['backBoneDihedralMax'], 
                                           beta1=self.plpStruct['backBoneBondAngleMin'], 
                                           beta2=self.plpStruct['backBoneBondAngleMax'])
        
        # work out way along the backbone in sets of three. Each triad forms the template against 
        # which we will align the forked carbon triad of the stapler.  
        # start with the first three backbones
        backBoneTriadIndices=[0,1,2] # the first three entries in the backbone 
        for brush in brushList:

            # Create some alternative name lists to help color the structures in blender. Useful debugging info. 
            bbNames = [ 'C' if not i in backBoneTriadIndices else 'B' for i  in range(len(backBoneXYZ)) ]
            stapBBNames = [ sym if not i in stapler['staplerBackBoneJoinIndices'] else 'Be' for i, sym in enumerate(stapler['staplerSymbols']) ]

            # create a copy of the stapler XYZ oriented in the current backbone position
            # creates a copy of the stapler positions    
            bbAlignedStaplerXYZ = self.aligner.align(backBoneXYZ, 
                                                     stapler['staplerXYZ'], 
                                                     backBoneTriadIndices, 
                                                     stapler['staplerBackBoneJoinIndices'],
                                                     alignPoint=1)
            # useful debug files
            fIO.saveXYZList(backBoneXYZ, bbNames, "bbBoneXYZ.xyz")
            fIO.saveXYZList(bbAlignedStaplerXYZ, stapBBNames, "alignedStapler.xyz")
            
            # align the join points of the brush backbone with join points of the aligned stapler, 
            # which is now the static party. 
            brush['xyz'] = self.aligner.align(bbAlignedStaplerXYZ,
                                              brush['xyz'], 
                                              stapler['staplerBrushJoinIndices'],
                                              brush['joinIndices'],
                                              alignPoint=1)

            # useful debug files
            fIO.saveXYZList(brush['xyz'], brush['symbols'], "alignedBrush.xyz")
            
            # The aligned staple points can now be added to the brush structure
            # in canonical form. The brush dictionary will be altered in place.
            self.InsertStaplerPointsIntoBrush(brush, stapler, bbAlignedStaplerXYZ)

            # dump the aligned full brush to debug file. 
            fIO.saveXYZList(brush['xyz'], brush['symbols'], "alignedBrushFullStruc.xyz")
            
            # move on to the next backbone triad
            backBoneTriadIndices = [ i + 3 for i in backBoneTriadIndices ]

     
    def leapItUp(self, name):
        rawPDBName = name + '.pdb'
        leapInFilename = "tleap.in"
        mdinFilename = "mdin"
        mdoutFilename = "mdout"
        prmtopFilename =  name + ".prmtop"
        coordsFilename = name + ".coords"
        rstFilename = name + ".rst"
        pdbleapFilename = name + "_tleap.pdb"
        pdbminFilename = name + "_min.pdb"
        bashscriptname = name + '.sh'
        
        vst=open(leapInFilename,'w', newline='\n')
        vst.write("source leaprc.protein.ff14SB\n")
        vst.write("mol=loadpdb " + rawPDBName + "\n")
        vst.write("saveamberparm mol " + prmtopFilename + " " + coordsFilename + "\n")
        vst.write("savepdb mol " + pdbleapFilename + "\n")
        vst.write("quit" + '\n')
        vst.close()
        
        vst=open(mdinFilename,'w', newline='\n')
        vst.write(" minimise " + rawPDBName + "\n")
        vst.write(" &cntrl\n")
        vst.write("  imin = 1,\n")
        vst.write("  maxcyc =500,\n")
        vst.write("  ncyc = 250,\n")
        vst.write("  ntb = 0,\n")
        vst.write("  igb = 0,\n")
        vst.write("  cut = 12\n")
        vst.write(" /\n")
        vst.close()

        vst=open(bashscriptname, "w", newline='\n')
        vst.write("#!/bin/bash\n")
        vst.write("tleap -f " + leapInFilename + "\n")
        vst.write("sander -O -i " + mdinFilename + " -o " + mdoutFilename + " -p " + prmtopFilename + " -c " + coordsFilename + " -r " + rstFilename + "\n")
        vst.write("ambpdb -p " + prmtopFilename + " -c " + rstFilename + " > " + pdbminFilename + "\n")
        vst.close()

        return pdbleapFilename, pdbminFilename, bashscriptname


    def generateBrushListFromPDB(self, brushList, filename):
        
        # read in the atoms from the given pdb file
        atoms = self.readPDBAtomsRet(filename)

        # trim the H's
        atoms = [ atom for atom in atoms if 'H' not in atom[1] ] # not in ['HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23']]
        
        # find the indices in the atom list where each brush ends. 
        atomTerIndices = [ i for i, atom in enumerate(atoms) if atom[1]=='OXT']
        
        # find the indices in the atom list where each brush starts.
        atomStartIndices = [ i + 1 for i in atomTerIndices ]
        atomStartIndices.insert(0, 0)
        
        # break up the atoms into subsets that correspond to the brushes.
        brushAtoms = [ atoms[i:j+1] for i, j in zip(atomStartIndices, atomTerIndices)]

        # now go through each brush and repopulate each brush with the new atom information.
        # can use the residue IDs from the newAtoms. Replaces the brushList in Place.          
        for brush, newAtoms in zip(brushList, brushAtoms):

            # copy across the information from newAtoms to brush one attribute at a time. 
            brush['names'] = [atom[1] for atom in newAtoms ]
            brush['xyz'] = [ np.array([float(atom[7]), float(atom[8]), float(atom[9])]) for atom in newAtoms ]
            brush['symbols'] = [ atom[13] for atom in newAtoms ]
            brush['resIds'] = [ atom[5] for atom in newAtoms ]                          
            brush['resNames'] = [ atom[3] for atom in newAtoms ]
            
            # figure out the ResId of the residue that we are fixing to the backbone
            attachmentResId = brush['stapleRes'] + brush['resIds'][0] - 1
            
            # get the index of the first atom with that residue id
            attachmentResFirstAtomIndex = brush['resIds'].index(attachmentResId) 

            # Set up the sets of names that are useful and depend on which residue is the linker
            if brush['resNames'][attachmentResFirstAtomIndex]=='LEU':
                joinNames = ['CD1', 'CG', 'CD2'] # the order of these names is crucial as it sets up the pairwise connectivity with the backbone
                azNames =['CD1', 'CD2'] # atoms in this residue for doing azimuthal rotation
                phaseNames =['CG', 'CB'] # atoms in this residue for doing phase rotation
                phaseFreezeNames = ['CD1', 'CG', 'CB', 'CD2'] # atoms in this residue that do not rotate when phase rotation happens
            elif brush['resNames'][attachmentResFirstAtomIndex]=='VAL':
                joinNames = ['CG1', 'CB', 'CG2'] # order crucial
                azNames =['CG1', 'CG2']
                phaseNames =['CB', 'CA']
                phaseFreezeNames = ['CG1', 'CB', 'CA', 'CG2']
            else:
                print("Identified Staple Res is not Leu or Val. Something's wrong.")
                sys.exit()
            
            # convert the names to indices
            # order crucial for joinIndices
            brush['joinIndices'] = [ i for jName in joinNames for i, name in enumerate(brush['names']) if name==jName and brush['resIds'][i]==attachmentResId ] 
            brush['azIndices'] = [ i for i, name in enumerate(brush['names']) for azName in azNames if name==azName and brush['resIds'][i]==attachmentResId]
            brush['phaseIndices'] = [ i for i, name in enumerate(brush['names']) for phaseName in phaseNames if name==phaseName and brush['resIds'][i]==attachmentResId ] 
            brush['phaseFreezeIndices'] = [ i for i, name in enumerate(brush['names']) for pFName in phaseFreezeNames if name==pFName and brush['resIds'][i]==attachmentResId ]
            brush['CAAttachResIndex'] = [ i for i, name in enumerate(brush['names']) if name=='CA' and brush['resIds'][i]==attachmentResId ][0] 
            #print(brush['CAAttachResIndex'], brush['names'][brush['CAAttachResIndex']], brush['resIds'][brush['CAAttachResIndex']])
            #print('OK')
    
            
    def PolymerisePLPRotate(self, brushList):
        ''' This function takes each of the brushes, which have now been energy minimised. 
        Each possesses the valine /leucine fork, which can be aligned on the back bone directly 
        in a single step.  
        
        We can then rotate the entire brush backbone about the CB-CA (val) or CG-CB (Leucin) axis of the 
        aligned valine or leucine. This sets the local phase rotation of the brush about an 
        sticking out from the backbone.
        
        Then we rotate the brush backbone about the CG1 and CG2 axis (valine), or 
        CD1 and CD2 axis (leucine) by the azimuth angle specied in the input numbers.  
        This allows us to control the azimuthal position of each brush on the polymer.  
        
        We do this for each brush - including the little ones which are spacers between the big ones.
        
        The updated brush structures are modified in place in the brushList for the output.'''
        numBrushes = len(brushList)

        # generate a chain of points which is 3 times the number of brushes a carbon bond length apart.
        backBoneXYZ  = self.makeRandomCoil(3 * numBrushes, 
                                           'RandomPolymer.txt',
                                           bondLength=1.4, 
                                           alpha1=self.plpStruct['backBoneDihedralMin'], 
                                           alpha2=self.plpStruct['backBoneDihedralMax'], 
                                           beta1=self.plpStruct['backBoneBondAngleMin'], 
                                           beta2=self.plpStruct['backBoneBondAngleMax'])
        
        # work out way along the backbone in sets of three. Each triad forms the template against 
        # which we will align the forked carbon triad of the stapler.  
        # start with the first three backbones
        backBoneTriadIndices=[0,1,2] # the first three entries in the backbone 
        for brush in brushList:

            # Create some alternative name lists to help color the structures in blender. Useful debugging info. 
            bbNames = [ 'C' if not i in backBoneTriadIndices else 'B' for i  in range(len(backBoneXYZ)) ]
            
            # align the join points of the brush with backbone triad 
            brush['xyz'] = self.aligner.align(backBoneXYZ,
                                              brush['xyz'], 
                                              backBoneTriadIndices,
                                              brush['joinIndices'],
                                              alignPoint=1)

            # useful debug files
            fIO.saveXYZList(backBoneXYZ, bbNames, "backBoneFinal.xyz")
            fIO.saveXYZList(brush['xyz'], brush['symbols'], "alignedBrushV2.xyz")
            
            # if the number of residues in the sequence is less than 2
            if len(set(brush['resIds']))<3:
                print('Warning: Sequence too short to orient effectively.')
            else:
                # Get the N-CA-C axis for the peptide
                backboneAxisBrush = [ xyz for i, xyz in enumerate(brush['xyz']) if brush['names'][i] in ['N', 'CA', 'C']]
                fIO.saveXYZ(backboneAxisBrush, "C", "backBoneBrush.xyz")
    
                # get some kind of unit vector along the brush in xyz space
                hAxis = self.getBrushVector(backboneAxisBrush)
                hAxisHat = hAxis/np.linalg.norm(hAxis)
    
                # get TNB at backbone - a local frame defined in terms of the coords of the backbone
                # this function was used to design TNB at end of polymer. Putting backbone
                # triad in, in unusual configuration to generate T aligned with the Azimuth azis.
                P0 = backBoneXYZ[backBoneTriadIndices[1]]
                P1 = backBoneXYZ[backBoneTriadIndices[0]]
                P2 = backBoneXYZ[backBoneTriadIndices[2]]
                TNB = coords.constructTNBFrame(P0, 
                                               P1, 
                                               P2)
    
                # user provides a vector in terms of thea and phi in the TNB axis locally. 
                # phi is rotation about T axis. 
                # theta is elevation above BN plane
                desiredBrushVecLabFrame = coords.generateTNBVecXYZ(TNB, 
                                                                   brush['Phase']* np.pi/180.0, 
                                                                   brush['Az']* np.pi/180.0) 
                dBVLabFrameHat = desiredBrushVecLabFrame/np.linalg.norm(desiredBrushVecLabFrame)
    
                # Rotation Axis is cross product of hAxis and desiredBrushVecLabFrame
                nHat = np.cross(hAxisHat, dBVLabFrameHat) 
    
                #Angle to Rotate is given by dot product of hAxis and desiredVector
                angle = np.arccos(np.dot(hAxisHat, dBVLabFrameHat))
    
                CAVec = brush['xyz'][brush['CAAttachResIndex']]
    
                vecs = [CAVec]
                names = ['H'] 
                
                hAxisList = [ (1 * 0.3 * i * hAxisHat) + CAVec for i in range(1,1000) ]
                hAxisNameList = 1000*['He']
                for v, n in zip(hAxisList, hAxisNameList):
                    vecs.append(v)
                    names.append(n)
    
                
                dbvList_CA = [ (0.3 * i * dBVLabFrameHat) + CAVec for i in range(1,10) ]
                dbvNameList1 = 10*['Li']
                for v, n in zip(dbvList_CA, dbvNameList1):
                    vecs.append(v)
                    names.append(n)
    
                nAxisList = [ (0.3 * i * nHat) + CAVec for i in range(1,10) ]
                nAxisNameList = 100*['Be']
                for v, n in zip(nAxisList, nAxisNameList):
                    vecs.append(v)
                    names.append(n)
                
                TNB0List = [ (0.3 * i * TNB[0]) + P1 for i in range(1,10) ]
                TNB1List = [ (0.3 * i * TNB[1]) + P1 for i in range(1,10) ]
                TNB2List = [ (0.3 * i * TNB[2]) + P1 for i in range(1,10) ]
                TNB0NameList = 10*['B']
                TNB1NameList = 10*['C']
                TNB2NameList = 10*['N']
                for v, n in zip(TNB0List, TNB0NameList):
                    vecs.append(v)
                    names.append(n)
                for v, n in zip(TNB1List, TNB1NameList):
                    vecs.append(v)
                    names.append(n)            
                for v, n in zip(TNB2List, TNB2NameList):
                    vecs.append(v)
                    names.append(n)            
    
                dbvList_TNB = [ (0.3*i * dBVLabFrameHat) + P1 for i in range(1,10) ]
                dbvNameList2 = 10*['O']
                for v, n in zip(dbvList_TNB, dbvNameList2):
                    vecs.append(v)
                    names.append(n)            
                
                fIO.saveXYZList(vecs, names, "debugPoints.xyz")
    
                # rotate all points in the brush except the frozen points in the backbone connector by 
                # an angle angle around the nHat vector.             
                for i,p in enumerate(brush['xyz']): 
                    if not i in brush['phaseFreezeIndices']:
                        brush['xyz'][i] = self.aligner.rotPAboutAxisAtPoint(p, brush['xyz'][brush['CAAttachResIndex']], nHat, angle) 
                        
                
                
                # Control the "phase rotation" of the brush - it's angle about a normal to the backbone. 
                # All points in the brush are rotated by the angle Phase 
                # about the axis formed by two atoms in the new part of the brush, that point away from the backbone. 
                # CG and CB of Leucine and CA and CB if valine is stapler. It is the rotation about the bond
                # linking the brush to the backbone.  
                #P0 = brush['xyz'][brush['phaseIndices'][0]]
                #P1 = brush['xyz'][brush['phaseIndices'][1]]
    
                # rotate the points in brush['xyz'] around the bond that joins to the BB.
                # don't rotate the anchor points involved in the connection to the backbone. 
                #for i,p in enumerate(brush['xyz']): 
                #    if not i in brush['phaseFreezeIndices']:
                #        brush['xyz'][i] = self.aligner.rotPAboutAxisBetweenPoints(p, P0, P1, brush['Phase']* np.pi/180.0) 
    
                # dump useful file for debugging            
                fIO.saveXYZList(brush['xyz'], brush['symbols'], "alignedBrushXYZ_AzEl.xyz")
    
                # Control the azimuthal rotation of the brush - its azimuthal position around the backbone. 
                # All points in the brush are rotate by the angle Az 
                # about the axis formed by the two valine fork carbons offered by the brush to the backbone. 
                # It is the rotation of the entire brush about the axis of the backbone. 
                #P0 = brush['xyz'][brush['azIndices'][0]]
                #P1 = brush['xyz'][brush['azIndices'][1]]
                #brush['xyz'] = [ self.aligner.rotPAboutAxisBetweenPoints(p, P0, P1, brush['Az']* np.pi/180.0) for p in brush['xyz'] ] 
    
                # dump useful file for debugging
                #fIO.saveXYZList(brush['xyz'], brush['symbols'], "alignedBrushXYZ_phaseAzRot.xyz")
                
            # move on to the next backbone triad
            backBoneTriadIndices = [ i + 3 for i in backBoneTriadIndices ]

    
    # The aligned staple points can now be added to the brush structure
    # in canonical PDBorder. The brush dictionary will be altered in place.
    def InsertStaplerPointsIntoBrush(self, brush, stapler, bbAlignedStaplerXYZ):
        # need to update the names arrays, XYZ arrays  symbols, resIds, resnames arrays, joiningIndices)
        # Must insert the CB, CG1, CG2, CD1, CD2 between CA and CB in the stapleRes'th residue
        brushIndexToInsert = 3 * (brush['stapleRes'] - 1) + 2
        brushIndexResIdsToCopy = 3 * (brush['stapleRes'] - 1) 
        
        for name in stapler['staplerNamesToInsert']: # names are in correct order already
            staplerIndex = cp.copy(stapler['staplerNames'].index(name))
            brush['xyz'].insert(brushIndexToInsert, cp.copy(bbAlignedStaplerXYZ[staplerIndex]))
            brush['names'].insert(brushIndexToInsert, cp.copy(stapler['staplerNames'][staplerIndex]))
            brush['symbols'].insert(brushIndexToInsert, cp.copy(stapler['staplerSymbols'][staplerIndex]))
            brush['resIds'].insert(brushIndexToInsert, cp.copy(brush['resIds'][brushIndexResIdsToCopy]))
            brush['resNames'].insert(brushIndexToInsert, cp.copy(brush['resNames'][brushIndexResIdsToCopy]))
        brush['joinIndices'][0] = brush['names'].index('N') 
        brush['joinIndices'][1] = brush['names'].index('CA')
        brush['joinIndices'][2] = brush['names'].index('C')
    
    def getBrushVector(self, vecList):
        
        if len(vecList)<=5:
            # outvec is the mean of the vectors from the first point to all the others. 
            outVec = np.mean([ a - vecList[0] for a in vecList[5:] ], axis=0)
        else:
            # take average of first 5 points in the list and compute 
            # mean to the rest of the list. 
            P0 = np.mean([ a for a in vecList[0:5] ], axis=0)
            outVec = np.mean([ a - P0 for a in vecList[5:] ], axis=0)
        # take average vector between P0 points to get some sort of axial vector
        return outVec
        
    def makeRandomCoil(self, N, filename='RandomPolymer.txt', minDist=1.5, bondLength=1.5, alpha1=30, alpha2=60, beta1=140, beta2=180):
        brushBBG = RPPBBG(filename)
        Z1 = minDist
        R1 = 3 * minDist
        R2 = N * 3 * minDist
        Z2 = N * bondLength
        
        polyStart = np.array([ 0.0, 0.0, Z1])
        envelopeList = ['frustum ' + str(Z1 - minDist) + ' ' + str(R1) + ' ' + str(Z2) + ' ' + str(R2)]
        brushBB = brushBBG.generateBuildingBlock( N,
                                                  polyStart,
                                                  alpha1,
                                                  alpha2,
                                                  beta1,
                                                  beta2,
                                                  minDist,
                                                  bondLength,
                                                  envelopeList=envelopeList,
                                                  visualiseEnvelope=(0, 100, 'envelope.xyz'))
        return brushBB.blockXYZVals
            
    def makePeptideBackBone(self, N, filename='alphahelix.txt'):
        brushBBG = PBG(filename) # generate a betastrand or an alphahelix (if you want hair pin, just generate a random polymer)
        brushBB = brushBBG.generateBuildingBlock(N, showBlockDirector=False, nameCA=True)
        return brushBB.blockXYZVals

    def lookupSym(self, name):
        retVal = ''
        if 'H' in name:
            retVal='H'
        elif 'C' in name:
            retVal='C'
        elif 'O' in name:
            retVal='O'
        elif 'N' in name:
            retVal='N'
        else:
            retVal=''
        return retVal    

    # list of individual brushes and filename to write as a single PDB model            
    def writePDB(self, brushList, filename, conectList=[], terminateBrushes = False, chainLet=' '):
        brushIndex = 0
        previousAtoms = 0
        atomList = []
        # loop through each chain in the list (each brush)
        for brush in brushList:
            # count the atoms in the chain
            numAtomsInBrush = len(brush['names'])
            
            # create an atom array
            for atomIndex in range(numAtomsInBrush):
                atom = [ previousAtoms + atomIndex + 1,
                         brush['names'][atomIndex],
                         ' ',
                         brush['resNames'][atomIndex],
                         chainLet,# brushLets[brushIndex],
                         brush['resIds'][atomIndex],
                         ' ',
                         brush['xyz'][atomIndex][0],
                         brush['xyz'][atomIndex][1],
                         brush['xyz'][atomIndex][2],
                         1.0,
                         0.0,
                         ' ',
                         brush['symbols'][atomIndex],
                         0]
                atomList.append(atom)
            previousAtoms += numAtomsInBrush 
            if terminateBrushes:
                # output a record that will trigger a TER statement at the end of each brush
                # remembers the atomIndex and brush information from immediately prior loop.
                # helps tleap out and gets round the fact that we have more chains than 26 letters
                atomList.append(['TER',
                                 previousAtoms + 1,  
                                 brush['resNames'][atomIndex],
                                 ' ',
                                 brush['resIds'][atomIndex]])
                previousAtoms += 1
            # increment the brush index
            brushIndex += 1
            print("Building: ", brushIndex)
        # write PDB file
        self.writePDBAtomsFromList(atomList, filename, conectList=conectList)
    
        return atomList
    
class GenericPLP():
    def __init__(self, backboneLengths, backbones, sequences, secondaries, PLPName, doProlineKink=False, BBAlphas=[(40,50)], BBBetas=[(165,185)], BrAlphas=[(40,50)], BrBetas=[(145,155)], outMode='blender_xyz' ):

        polymerBrushDict={}
        polymerBrushDict['backbones'] = []
        polymerBrushDict['brushes'] = []
        polymerBrushDict['connectors'] = []
        
        if len(BBAlphas)==1:
            BBAlphas = len(backbones) * BBAlphas
        if len(BBBetas)==1:
            BBBetas = len(backbones) * BBBetas
        if len(BrAlphas)==1:
            BrAlphas = len(backbones) * BrAlphas
        if len(BrBetas)==1:
            BrBetas = len(backbones) * BrBetas
        
        for sequence, backbone, length, alphaBBRange, betaBBRange, alphaBrRange, betaBrRange, secondary in zip(sequences, backbones, backboneLengths, BBAlphas, BBBetas, BrAlphas, BrBetas, secondaries):
            
            # set up the dictionaries
            backboneDict = {}
            brushDict = {}

            # pack the dictionaries
            backboneDict['filename'] = 'RandomPolymer.txt'  
            backboneDict['mode'] = 'Polymer'
            backboneDict['name'] = 'BlockA'

            if outMode=='VMD_pdb':
                # here can set up the 'polymer unit'
                backboneDict['residueLength'] = 2
                backboneDict['numMonomers'] = backboneDict['residueLength'] * length
                
                backboneDict['monomerNames'] = length * ['C' + str(a + 1) for a in range(backboneDict['residueLength']) ]
                backboneDict['residueNames'] = length * [backbone]
                brushDict['indexOfFirstBrush'] = 1 # sets which backbone atom the first brush comes out of  
                numUnitsToSkip = backboneDict['residueLength']
                brushDict['spacing'] = 1 + 2 * numUnitsToSkip   # spacing is number of particles between particles with a brush    
                
            elif outMode=='blender_xyz':
                backboneDict['residueLength'] = 1
                backboneDict['numMonomers'] = length
                backboneDict['monomerNames'] = length * [backbone]
            else:
                backboneDict['residueLength'] = 1
                backboneDict['numMonomers'] = length 
                backboneDict['monomerNames'] = length * [backbone]

            backboneDict['alpha1'] = alphaBBRange[0]
            backboneDict['alpha2'] = alphaBBRange[1]
            backboneDict['beta1'] = betaBBRange[0]
            backboneDict['beta2'] = betaBBRange[1]
            backboneDict['minDist'] = 1.3
            backboneDict['bondLength'] = 1.4
            backboneDict['Z1'] = 2
            backboneDict['R1'] = 20
            backboneDict['Z2'] = 1.5 *backboneDict['numMonomers'] * backboneDict['bondLength'] + backboneDict['Z1']
            backboneDict['R2'] = 70
        
            brushDict['name'] ='brush'
 
            # choose whether to output brush as a pdb for vmd or an xyz for blender
            if outMode=='VMD_pdb':
                if secondary=='Beta':
                    brushDict['filename'] = 'betastrand.txt'
                    brushDict['mode'] = 'PeptideRandom'
                elif secondary=='Alpha':
                    brushDict['filename'] = 'alphahelix.txt'
                    brushDict['mode'] = 'PeptideRandom'
                else:
                    brushDict['filename'] = 'alphahelix.txt'
                    brushDict['mode'] = 'PeptideAlternate'

                # in VMD mode, numMonomers = numResidues and generates an N, CA, C triplet in configuration specified in file
                brushDict['residueLength'] = 3
                brushDict['numMonomers'] = len(sequence)         
                brushDict['monomerNames'] =  len(sequence) * ['N', 'CA', 'C'] # set up to be names of final sequence
                brushDict['residueNames'] = [ residue for residue in sequence ] # add this to generate the PDB
                # set up brushes on the CAs. Can mess with these to force it to skip every other residue etc. 
                brushDict['indexOfFirstBrush'] = 1
                # spacing is number of particles between particles with a brush 
                brushDict['spacing'] = 2
            elif outMode=='blender_xyz':
                brushDict['residueLength'] = 1
                brushDict['filename'] = 'RandomPolymer.txt'  
                brushDict['mode'] = 'Polymer'
                # in blender mode, eachpoint represents a single residue
                brushDict['numMonomers'] = len(sequence)         
                brushDict['monomerNames'] = [ self.convertPeptideToElementsForBlender(peptide) for peptide in sequence ] # if specified in peptide mode then overrides the C-CA-N peptide names            
                brushDict['indexOfFirstBrush'] = 0
                brushDict['spacing'] = 1
            else:
                brushDict['residueLength'] = 1
                brushDict['filename'] = 'RandomPolymer.txt'  
                brushDict['mode'] = 'Polymer'
                brushDict['numMonomers'] = len(sequence)         
                brushDict['monomerNames'] = [ self.convertPeptideToElementsForBlender(peptide) for peptide in sequence ] # if specified in peptide mode then overrides the C-CA-N peptide names
                brushDict['indexOfFirstBrush'] = 0
                brushDict['spacing'] = 1
            
            brushDict['phaseRange'] = 180.0 
            brushDict['alpha1'] = alphaBrRange[0]
            brushDict['alpha2'] = alphaBrRange[1]
            brushDict['beta1'] = betaBrRange[0]
            brushDict['beta2'] = betaBrRange[1]
            brushDict['minDist'] = 1.0
            brushDict['bondLength'] = 2.0
            brushDict['Z1'] = 2
            brushDict['R1'] = 50
            brushDict['Z2'] = 1.5 * brushDict['numMonomers'] * brushDict['bondLength'] + brushDict['Z1']
            brushDict['R2'] = 50
            brushDict['doProlineKink']=doProlineKink
            brushDict['prolineKinkAngle']=90 * np.pi/180

            polymerBrushDict['backbones'].append(cp.copy(backboneDict))
            polymerBrushDict['brushes'].append(cp.copy(brushDict))

        # set up any connectors that we need (loop through the first N-1 backbones
        for backbone in (polymerBrushDict['backbones'])[0:-1]:
            connectorDict = {}
            connectorDict['displacement'] = backbone['bondLength']
            connectorDict['alpha'] = 0.0
            connectorDict['beta'] = 170.0

            polymerBrushDict['connectors'].append(cp.copy(connectorDict))
        
        self.xyz, self.names, self.numBrushesPerBlock = GBP(polymerBrushDict)

        # save an xyz file in any case for blender even if in PDB mode. 
        fIO.saveXYZList(self.xyz, self.names, PLPName + ".xyz")

        residueLengths = [backboneDict['residueLength'] for backboneDict in polymerBrushDict['backbones']]
        for brushDict in polymerBrushDict['brushes']:
            residueLengths.append(brushDict['residueLength'])


        # if in vmd mode, then write a pdb to file
        if outMode=='VMD_pdb':
            self.writePDB(self.xyz, self.names, backboneLengths, backbones, sequences, residueLengths, self.numBrushesPerBlock, secondaries, PLPName + ".pdb")
            
    def writePDB(self, xyz, names, backBoneLengths, backboneNames, brushResidueNames, residueLengths, numBrushesPerBlock, secondaries, filename):
        # create empty PDB obj
        PDBObj = PDBL()
        numBlocks = len(backBoneLengths)
        
        # for each atom figure out the name of the residue it lives in 
        residues = [ bbLength * 2 * [ bbName ] for bbLength, bbName in zip(backBoneLengths, backboneNames) ]  
        for brushNames, numBrushes in zip(brushResidueNames, numBrushesPerBlock):
            rNames = []
            for resName in brushNames:
                rNames.append(resName)
                rNames.append(resName)
                rNames.append(resName)
            residues.append( numBrushes * rNames)
        residues = [ item for sublist in residues for item in sublist ]
        
        # for each atom figure out the number of the residue it lives in
        resNums = []
        resNum = 1
        chainLets='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        chains = []
        chainIndex = 0
        # first deal with the backbones
        for blockNum, blockLength in enumerate(backBoneLengths):
            for _ in range(blockLength):
                for _ in range(residueLengths[blockNum]):
                    resNums.append(resNum)
                    chains.append(chainLets[chainIndex])
                resNum+=1
            chainIndex+=1
        
        blockNum = numBlocks
        # loop through each block, getting the brush length for that block
        for brushNum, nameList in zip(numBrushesPerBlock, brushResidueNames):
            # loop through each brush
            for _ in range(brushNum):
                for _ in range(len(nameList)):
                    # loop through each residue outputting the resnum
                    for _ in range(residueLengths[blockNum]):
                        resNums.append(resNum)
                        chains.append(chainLets[chainIndex])
                    # increment res num at end of res
                    resNum+=1
                chainIndex+=1
            blockNum+=1
        print(len(resNums), len(residues))
        
        # get the symbol name for each atom
        sym = [ self.lookupSym(name) for name in names ]

        # figure out the chain for each atom
        

        
        # create each atom line
        atoms = [ [i+1, names[i], ' ', PDBObj.ConvertAACodes(residues[i]), chains[i], resNums[i], ' ', xyz[i][0], xyz[i][1], xyz[i][2], 1.0, 0.0, ' ', sym[i], 0] for i in range(len(names))]


        # create each helix and sheet in the system

        # for each block we know the number of helices = numBrushes in that block
        helices = []
        sheets = []
        helNum = 1
        sheetNum = 1
        blockAlpha = ['A', 'B', 'C', 'D', 'E']
        blockIndex = 0 
        # for the PLPs the convention is 
        # block A: Backbone residues 1 to backboneLengths[0]
        # block B: Backbone residues backboneLengths[0]+1 to backboneLengths[0] + backboneLengths[1]
        # then the brushes start
        # 
        # get the residue index of the first brush
        brushesResNumsStart = np.sum(backBoneLengths) + 1
        # step through the number of blocks 
        for secondary, numBrushes in zip(secondaries, numBrushesPerBlock):
            # check to see what kind of secondary structure we are adding
            if secondary=='Alpha':
                for brushNum in range(numBrushes):
                    helData = [helNum, 
                               blockAlpha[blockIndex] + str(brushNum),
                               PDBObj.ConvertAACodes(brushResidueNames[blockIndex][0]),
                               'A',
                               brushesResNumsStart + brushNum * len(brushResidueNames[blockIndex]),
                               ' ',
                               PDBObj.ConvertAACodes(brushResidueNames[blockIndex][-1]),
                               'A',
                               brushesResNumsStart + (brushNum + 1) * len(brushResidueNames[blockIndex]) - 2,
                               ' ',
                               1,
                               '0123456789012345678901234',
                               len(brushResidueNames[blockIndex])]
                    helices.append(helData)
                    helNum += 1
            if secondary=='Beta':
                pass
            
            blockIndex += 1
            
        # create empty PDB object (no filename)
        
        
        # write PDB file
        PDBObj.writePDBAll(atoms, helices, sheets, "bb_" + filename)

        # write without backbone
        PDBObj.writePDBAll([ atom for atom in atoms if atom[3] not in backboneNames ], helices, sheets, filename )

        # generate a tleaped file        
        leap = leapify()
        print(leap.leapItUp(filename[0:-4]))
    

    def lookupSym(self, name):
        retVal = ''
        if 'C' in name:
            retVal='C'
        elif 'O' in name:
            retVal='O'
        elif 'N' in name:
            retVal='N'
        else:
            retVal=''
        return retVal
    
    def getVals(self):
        return self.xyz, self.names

    def convertPeptideToElementsForBlender(self, peptide):
        peptides = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B']
        elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'Fl', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ga', ]
        return elements[ peptides.index(peptide) ]

if __name__=='__main__':
    AllthePLPS(sys.argv[1])
    # gPLP = GenericPLP_Sequence(sys.argv[1])
    # GenericPLP([15,30], ['PNa','PNb'], ['YYYYYPYYYYYYYYRRR', 'GSGSGSGSGS' ], ['Alpha', 'Beta'], 'delta', doProlineKink=True, outMode='VMD_pdb')
    # GenericPLP([10,10], ['PNa', 'PNb'], ['GGGGGGGGGG', 'VVVVVVVVVV'], ['Alpha', 'Beta'], 'om', doProlineKink=False, outMode='VMD_pdb')
    print("example done")