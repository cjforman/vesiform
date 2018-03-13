'''
Created on 14 Dec 2017

@author: chris forman
'''
import numpy as np

# ***** IO Functions ******

def ct2xyz(filename):
    atomLines = readTextFile(filename)
    
    # get the num of atoms in the ct file
    numAtoms = int(atomLines[1].split()[0])
        
    # parse each line in the ct file into separate values 
    xyzListAll = [ line.split() for line in atomLines[2:numAtoms+2]]
        
    # extract the names and xyz positions as floats into separate arrays
    atomNameList  = [ xyzEntry[3] for xyzEntry in xyzListAll]
    xyzList  = [ np.array([float(xyzEntry[0]), float(xyzEntry[1]), float(xyzEntry[2])])  for xyzEntry in xyzListAll]

    # check that the number of lines read in match the stated length in the file                 
    if (len(xyzList) != numAtoms):
            print("Warning: Num atoms read doesn't match num atoms specified in file")
            
    saveXYZList(xyzList, atomNameList, filename[0:-3]+'.xyz')

def saveXYZ(vectors, atom, filename):
    zLines = [str(len(vectors)) + '\n']
    zLines.append("Test Block\n")
    for A in vectors:
        zLines.append(atom + ' ' + str(A[0]) + ' ' + str(A[1]) + ' ' + str(A[2]) +   '\n')
    writeTextFile(zLines, filename)

def loadXYZ(filename):
    atomLines = readTextFile(filename)
    
    # get the length of the xyz file
    numAtoms = int(atomLines[0])
        
    # parse each line in the xyz file into separate values 
    xyzListAll = [ line.split() for line in atomLines[2:]]
        
    # extract the names and xyz positions as floats into separate arrays
    atomNameList  = [ xyzEntry[0] for xyzEntry in xyzListAll]
    xyzList  = [ np.array([float(xyzEntry[1]), float(xyzEntry[2]), float(xyzEntry[3])])  for xyzEntry in xyzListAll]

    # check that the number of lines read in match the stated length in the file                 
    if (len(xyzList) != numAtoms):
            print("Warning: Num atoms read doesn't match num atoms specified in file")

    return atomNameList, xyzList

def loadXYZFrames(xyz):
    ''' reads an xyz file and returns an array of frames. Each frame containes the list of atoms names and a list of corresponding coords as numpy coords.'''
    #read in the xyz file
    xyzData=readTextFile(xyz)

    #initialise frames output
    frames=[]

    curline=0
    while curline<len(xyzData):
        #read the frame length from the current line in xyzData
        framelen = int(xyzData[curline])
    
        #get the data for the current frame
        xyzList=xyzData[curline+2:curline+2+framelen]
    
        #increment the cur line index
        curline+=2+framelen
    
        #create the output lists for the atom names and atom coords
        atomNames=[]
        atomCoords=[]
        for atomData in xyzList:
            atom=atomData.split()
            atomNames.append(atom[0])
            atomCoords.append(np.array([float(atom[1]),float(atom[2]),float(atom[3])]))
        
        frames.append([atomNames,atomCoords])

    return frames

def saveEllipsoidXYZList(atomNames, ePositions, eSizes, eRs, eRotVecs, filename):
    zLines = [' ' + str(len(ePositions)) + '\n']
    zLines.append("Energy of minimum      1=      -12.8344616458 first found at step       95\n")
    for atomName, ePosition, eSize, eR, eRotVec in zip(atomNames, ePositions, eSizes, eRs, eRotVecs):
        line = ' ' + atomName 
        line += ' ' + str(ePosition[0]) + ' ' + str(ePosition[1]) + ' ' + str(ePosition[2])
        line += ' ellipse ' + str(eSize[0]) + ' ' + str(eSize[1]) + ' ' + str(eSize[2])
        line += ' ' + str(eR[0][0]) + ' ' + str(eR[1][0]) + ' ' + str(eR[2][0])
        line += ' ' + str(eR[0][1]) + ' ' + str(eR[1][1]) + ' ' + str(eR[2][1])
        line += ' ' + str(eR[0][2]) + ' ' + str(eR[1][2]) + ' ' + str(eR[2][2])
        line += ' atom_vector ' + str(eRotVec[0]) + ' ' + str(eRotVec[1]) + ' ' + str(eRotVec[2]) + '\n'
        zLines.append(line)
    writeTextFile(zLines, filename)


def saveXYZList(vectors, atomList, filename):
    zLines = [str(len(vectors)) + '\n']
    zLines.append("Test Block\n")
    for v, atom in zip(vectors, atomList):
        zLines.append(atom + ' ' + str(v[0]) + ' ' + str(v[1]) + ' ' + str(v[2]) +   '\n')
    writeTextFile(zLines, filename)

def readTextFile(filename):
    #read line data in from file
    try:
        vst = open(filename, 'r')
        lines = vst.readlines()
        vst.close()
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        return False

    return lines

def writeTextFile(lines, filename):
    # print("writing:", filename, " with ", len(lines)-2, " monomers.")
    # write line data to file
    try:
        with open(filename, 'w') as f:
            for line in lines:
                a=line
                f.write(a)
    except:
        pass
    return

def fileRootFromInfile(infile, root):
    if root in infile:
        return infile[0:-len(root)-1]
    else:
        return None
        
def dumpVecList(inpList, filename):
    writeTextFile([str(pos[0]) + ' ' + str(pos[1]) + ' ' + str(pos[2]) + '\n' for pos in inpList], filename) 

