#!/usr/bin/env python
import sys
import numpy as np
import Utilities.fileIO as fIO
import copy as cp

class PDB():
    ''' A class of utility functions for handling PDB files.'''
    
    def __init__(self, filename):
        #read in the atoms from the pdb
        print "Processing pdb filename: ", filename
        self.filename = filename
        self.readPDBAtoms(filename)
        
    def readPDBAtoms(self, filename):
        # parse a pdb file and saves all kinds of ways of looking at it as member variables that
        # are accessible for the operational functions
        self.rawPDB = fIO.readTextFile(filename)
        self.atoms = self.extractAtomsFromPDB()
        self.numAtoms = len(self.atoms)
        self.backBoneIndices = self.extractBackBoneIndices()
        
    def writePDBAtoms(self, filename):
        #open data file
        try:
            vst = open(filename, 'w')
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            raise Exception, "Unable to open output file: "+filename
    
        #parse data 
        for atom in self.atoms:
            l = self.pdbLineFromAtom(atom)
            vst.write(l)
        vst.close()
        return

    def extractHetAtomsFromPDB(self):
        #parse data 
        atoms = []
        for line in self.rawPDB:
            if line[0:6]=="HETATM":
                try:
                    a = self.parsePdbLine(line)
                    if not a in atoms:
                        atoms.append(a)
                except:
                    print "line: " + line + " not understood"
                    exit(0)
        print len(atoms), " HETATOM lines read from pdb file."                    
        return atoms


    def extractAtomsFromPDB(self):
        #parse data 
        atoms = []
        for line in self.rawPDB:
            if line[0:4]=="ATOM":
                try:
                    a = self.parsePdbLine(line)
                    if not a in atoms:
                        atoms.append(a)
                except:
                    print "line: " + line + " not understood"
                    exit(0)
                    
        print len(atoms), " ATOM lines read from pdb file."
        return atoms

    def parsePdbLine(self, line):
        l=line
        atom_seri = int(l[6:11])
        atom_name = l[12:16].split()[0]
        alte_loca = l[16]
        resi_name = l[17:20].split()[0]
        chai_iden = l[21]
        resi_numb = int(l[22:26])
        code_inse = l[26]
        atom_xcoo = float(l[30:38])
        atom_ycoo = float(l[38:46])
        atom_zcoo = float(l[46:54])
        try:
            atom_occu = float(l[54:60])
        except:
            atom_occu=0.0
    
        try:
            atom_bfac = float(l[60:66])
        except:
            atom_bfac=0.0    
        
        try:
            seg_id = l[72:76]
        except:
            seg_id=' '
    
        try:
            atom_symb = l[76:78].split()[0]
        except:
            try:
                atom_symb = l[68]
            except:
                atom_symb= ' '
    
        try:
            charge=l[78:80]
        except:
            charge=' '
    
        return [atom_seri, atom_name, alte_loca, resi_name, chai_iden, resi_numb, code_inse, atom_xcoo, atom_ycoo, atom_zcoo, atom_occu, atom_bfac,seg_id,atom_symb,charge]

    # write a pdb line from an atom string
    def pdbLineFromAtom(self, atom):
        try:
            l = 'ATOM {: >06d} {: <4}{:1}{:3} {:1}{: >4d}{:1}   {: >8.3f}{: >8.3f}{: >8.3f}{: >6.2f}{: >6.2f}      {: <4}{: >2}{: >2}\n'.format(int(atom[0]), atom[1], atom[2], atom[3], atom[4], int(atom[5]), atom[6], float(atom[7]), float(atom[8]), float(atom[9]), float(atom[10]), float(atom[11]), atom[12],atom[13],atom[14])
        except:
            print "unable to write atom to string: "
            print atom
            exit(0)
        return l

    def extractBackBoneIndices(self):
        # Extracts the index numbers of the backbone atoms
        return  [ atomIndex for atomIndex, atom in enumerate(self.atoms) if atom[1] in ['N', 'CA', 'C']]


    def makeXYZForBlenderFromPDB(self, outfile, backBoneOnly=True, residueTypeCheck=False):
        
        if backBoneOnly:
            # filter out the backbone atoms
            atom2XYZ = [ atom for atom in self.atoms if atom[1] in ['CA', 'C', 'N']]
            for atom in atom2XYZ:
                if atom[1]=='CA':
                    atom[1]='C'
        else:
            atom2XYZ = self.atoms
    
        # check to see if markGroups flag is set.
        if residueTypeCheck:
            # convert the atom type to an atom colour that reflects the functional group of the residue.
            for atom in atom2XYZ:
                if atom[3] in ['ARG', 'HIS', 'LYS']:
                    atom[1] ='P' # make the positive charge residues gold
                elif atom[3] in ['ASP', 'GLU']:
                    atom[1] ='O' # make the negative charge residues brown
                elif atom[3] in ['SER', 'THR', 'ASN', 'GLN']:
                    atom[1] ='C' # make the polar residues purple
                elif atom[3] in ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']:
                    atom[1] ='N' # make the hydrophobic residues blue
                elif atom[3] in ['CYS', 'SEC','GLY', 'PRO', 'HYP']:
                    atom[1] ='S' # make the special residues green
                else:
                    print 'Unrecognised residue name:', atom[3], ' line:', atom[0]
                    atom[1] ='Pb' # dark grey
    
        self.saveXYZFromPDB(atom2XYZ, outfile)
        return          

    def saveXYZFromPDB(self, atoms, outfile):
        lines = [str(len(atoms))+'\n']
        lines.append('XYZ data written by pdbProc\n')
        for atom in atoms:
            lines.append( str(atom[1]) + ' ' + str(atom[7]) + ' ' + str(atom[8]) + ' ' + str(atom[9]) +'\n')
    
        fIO.writeTextFile(lines, outfile)

    def insertXYZIntoPdb(self, xyzFilename, outfile):
        '''function reads in an xyz file and replaces the coords in the pdb with the coords in the xyz file.
           Generates am output PDB file for each frame in the xyz file.'''
        frameNum=0
        for frame in fIO.loadXYZFrames(xyzFilename):
            self.replacePdbAtoms(frame[1], outfile + '.' + str(frameNum) + '.pdb')
            frameNum+=1
        return

    def replacePdbAtoms(self, newCoords, outfile):
        ''' Creates a verbatim copy of the pdb data but with the coords replaced by the list of coords in newCoords.
            Assumes the coords are in the same order as the pdb file atom lines and are in a list of xyz triplets.'''
    
        #check the lists are compatible
        if len(self.atoms)!=len(newCoords):
            print "atom lists incompatible sizes" 
            sys.exit(0)
    
        #seed the output array
        outputlines=[]
        curCoord=0
    
        #loop through the raw input
        for line in self.rawPDB:
            #copy the line for output
            newLine=line[:]
    
            #check to see if the current line in the pdb is an atom or a hetatom
            if line[0:4]=="ATOM" or line[0:6]=="HETATM":
                #if it is then parse the input atom line
                try:
                    atom = self.parsePdbLine(newLine)
                except:
                    print "line: " + line + " not understood"
                    exit(0)
                #replace the coords
                atom[7] = newCoords[curCoord][0]
                atom[8] = newCoords[curCoord][1]
                atom[9] = newCoords[curCoord][2]
     
                curCoord += 1
                
                #generate the new line
                newLine = self.pdbLineFromAtom(atom)
           
            outputlines.append(cp.copy(newLine))
    
        fIO.writeTextFile(outputlines, outfile)
    
        return

    def replaceAtoms(self, newCoords):
        ''' Replaces the atomic coords in the atoms array with a new set of coords.'''
    
        if len(self.atoms)!=len(newCoords):
            print "atom lists incompatible sizes" 
            sys.exit(0)
    
        newAtoms = []
        
        #loop through the raw input
        for atom, newCoord in zip(self.atoms, newCoords):
            #copy the line for output
            newAtom = cp.copy(atom)
    
            #replace the coords
            newAtom[7] = newCoord[0]
            newAtom[8] = newCoord[1]
            newAtom[9] = newCoord[2]
             
            #append the data to the output array
            newAtoms.append(newAtom)
    
        return newAtoms

    def extractCoords(self):
        # returns numpy array of atoom coords from PDB
        return [np.array([atom[7],atom[8],atom[9]]) for atom in self.atoms]

    def extractAllAtomNames(self):
        return  [ atom[1] for atom in self.atoms]

    def breakChainIntoResidues(self, chain):
    
        # Create the output array
        residues=[]
    
        # Set the current Residue to the value of the first residue in the chain.
        curRes=chain[5]
        
        # Create the dummy variable to store a list of atoms for each residue.
        atomsInRes=[]
        
        # Go through the chain appending each atom to the dummy list.
        # If we find the first atom in a new residue then increment the residue number and 
        # add the previous residue to the outgoing list.
        # Reset the list for the new current residue
        for atom in chain:
            if atom[5]!=curRes:
                curRes+=1
                residues.append(atomsInRes)
                atomsInRes=[]
            atomsInRes.append(atom)
        
        #append the last residue
        residues.append(atomsInRes)
        return residues

    def readSequence(self, mode, outfile):
        residues = self.findResidues(self.atoms)
        
        #mode 1 numbered with three letter codes on separate lines
        if mode==1:
            print "mode: 1 selected"
            fIO.writeTextFile([str(res[0])+' '+str(res[1])+'\n' for res in residues],outfile)
        #mode 1 unumbered with three letter codes on separate lines - suitable for a modify sequence command
        if mode==2:
            print "mode: 2 selected"
            fIO.writeTextFile([str(res[1])+'\n' for res in residues],outfile)
        #mode 3 string of first letters only
        if mode==3:
            print "mode: 3 selected"
            l=''
            for res in residues:
                #if not a special case then copy first letter otherwise deal with each residue type on its own.
                if not res[1] in ['ARG','LYS','ASP','GLU','ASN','GLN','SEC','PHE','TYR','TRP']:
                    l=l+str(res[1][0])
                else:
                    if res[1]=='ARG':
                        l+='R'
                    if res[1]=='LYS':
                        l+='K'
                    if res[1]=='ASP':
                        l+='D'
                    if res[1]=='GLU':
                        l+='E'
                    if res[1]=='ASN':
                        l+='N'        
                    if res[1]=='GLN':
                        l+='Q'        
                    if res[1]=='SEC':
                        l+='U'
                    if res[1]=='PHE':
                        l+='F'                
                    if res[1]=='TYR':
                        l+='Y'
                    if res[1]=='TRP':
                        l+='W'
            l+='\n'
            fIO.writeTextFile(l,outfile)
        return
    
    def modifySequence(self, newSequence, startResidue, outFile):
        
        newSeq = fIO.readTextFile(newSequence)
        startResidue = int(startResidue)
    
        #trim off carriage returns - might be operating system dependent. works today!
        newSeq=[ newS[0:-1] for newS in newSeq]
        
        #open the outputfile
        fO=open(outFile,'w')
    
        #There are two list of atoms: an input (atoms) and an output (fO).
        curAtomInputIndex = 0
        curAtomOutputIndex = 1
     
        #Create an index for the replacement sequence (newSeq) 
        newResSeqIndex = 0  #always zero based for first entry
    
        #count the number of existing residues we have already processed. Initialise to zero.
        resCount = 0
    
        #keep note of number of the last residue that we processed
        lastResNumber = self.atoms[0][5] #initialise to the residue number of the first atom in the input list
    
        #create a variable to store the name of the last residue processed
        lastResName = self.atoms[0][3] #initialise to the residue name of the first atom in the input list
        
        #initialy we are not renaming until a set number of residues has been processed 
        renaming = False
        
        #copy across the raw data verbatim and dump to file.
        #if we hit an ATOM or HETATM decide whether or not to output the atom.
            #if we output the atom then replace the residue name and atom number 
        #if we hit a TER we must output with the name of the last redisue and increment the atomnumber    
        for line in self.rawPDB:
            #split the current line into tokens
            tokens=line.split()
    
            #if we encounter an atom line then decide if we need to modify it or ignore it.
            if tokens[0] in ['ATOM','HETATM']:
      
                #extract the data for the current atom
                atom = self.atoms[curAtomInputIndex]
                
                #increment the input index to the next atom for the next cycle
                curAtomInputIndex += 1
    
                #if the residue number has changed from last cycle then we have finished processing the previous residue.
                #Increment the number of residues processed.
                if atom[5]!=lastResNumber:
                    #make a note of the new res number and name
                    lastResNumber = atom[5]
                    lastResName = atom[3]
                    
                    #increment the number of residues processed
                    resCount+=1
                    
                    #if we are renaming then also increment the new sequence residue number
                    if renaming:
                        newResSeqIndex+=1
                        if newResSeqIndex>=len(newSeq):
                            renaming = False
                        else:
                            lastResName = newSeq[newResSeqIndex]                    
    
                #check the conditions for jumping into renaming mode.        
                if ((resCount>=startResidue-1) and (newResSeqIndex<len(newSeq))):
                    renaming=True
    
                #debug statement
                #print atom[5], newResSeqIndex, resCount, lastResNumber, lastResProc
    
                #Do we need to reprocess the current residue?
                if (renaming):
                    #Is the new sequence name different from the existing sequence?
                    if atom[3]!=newSeq[newResSeqIndex]:
                        #Only output atoms if they are backbone atoms, or belongs to an NME or ACE residue
                        if (str(atom[1]) in ['C','N','O','H','CA']) or (atom[3] in ['NME','ACE']):
                            #give the sequence a new name
                            atom[3] = newSeq[newResSeqIndex]
                            
                            #give the atom a new number
                            atom[0] = int(curAtomOutputIndex)
                            
                            #output the atom to file
                            fO.write(self.pdbLineFromAtom(atom))
                            
                            #increment the atom output index
                            curAtomOutputIndex+=1
                    else:
                        #new sequence is the same as the old one so just output the atom as is but with a new atom number
                        #give the atom a new number
                        atom[0]=int(curAtomOutputIndex)
                        #output the atom to file
                        fO.write(self.pdbLineFromAtom(atom))
                        #increment the atom output index
                        curAtomOutputIndex+=1
                else:
                    #We are not renaming the current residue so just output the atom with a new number
                    atom[0]=int(curAtomOutputIndex)
                    #output the line
                    fO.write(self.pdbLineFromAtom(atom))
                    #increment the atom output index
                    curAtomOutputIndex+=1
            else:
                if tokens[0] in ['TER']:
                    l='TER  {: >06d}      {:3} {:1}{: >4d}\n'.format(curAtomOutputIndex, lastResName, tokens[3], int(tokens[4]))
                    fO.write(l)
                    curAtomOutputIndex+=1       
                else:
                    #write the input line to file without mods
                    fO.write(line)
    
        #close file when there is no more input data
        fO.close()
    
        return

    def atomsToCOM(self):
        '''sets an atoms array to its centre of mass. Returns the COM and the new array'''
        newAtoms = []
        COM = np.array([0.0,0.0,0.0])
        for atom in self.atoms:
            COM += np.array([atom[7], atom[8], atom[9]])
        COM /= len(self.atoms)
        for atom in self.atoms:
            newAtom=[item for item in atom]
            newAtom[7]-=COM[0]
            newAtom[8]-=COM[1]
            newAtom[9]-=COM[2]
            newAtoms.append(newAtom)
            
        self.Com = COM
        return newAtoms
    
    
if __name__=="__main__":
    
    pdbfilename = sys.argv[1]
    
    # create instance of PDB object
    pdbObject = PDB(pdbfilename)
    infile = fIO.fileRootFromInfile(pdbfilename, "pdb")
    pdbObject.makeXYZForBlenderFromPDB(infile + ".xyz", backBoneOnly=True)
    print"Done"
    
    
    
        
    
    
