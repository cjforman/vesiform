import sys
from Utilities.keyProc import keyProc
import Utilities.cartesian as cart
import Utilities.coordSystems as coords 
import Utilities.fileIO as fIO
import numpy as np
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG

class spidroinBackboneGenerator(BBG):
    # This is a class for generating a spidroin coarse grain backbone. 
    
    def __init__(self, paramFilename):
        # initialise the parameter dictionary
        keyProc.__init__(self, paramFilename)    
    
    def initialiseParameters(self):

        keyProc.initialiseParameters(self)
    
        self.PPbondLength = self.getParam('PPbondLength')
        self.QQBondLength = self.getParam('QQbondLength')
        self.GG1bondLength = self.getParam('GG1bondLength')
        self.GG2bondLength = self.getParam('GG2bondLength')
        self.CNBondLength = self.getParam('CNbondLength')
        self.CAngle = self.getParam('CAngle')
        self.NAngle = self.getParam('NAngle')
        self.GPPsi = self.getParam('GPPsi')
        self.GPPhi = self.getParam('GPPhi')
        self.PGPsi = self.getParam('PGPsi')
        self.PGPhi = self.getParam('PGPhi')
        self.GQPsi = self.getParam('GQPsi')
        self.GQPhi = self.getParam('GQPhi')
        self.QGPsi = self.getParam('QGPsi')
        self.QGPhi = self.getParam('QGPhi')
        self.SP1NumGUnits = self.getParam('SP1NumGUnits')
        self.SP2NumGUnits = self.getParam('SP2NumGUnits')
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')
        
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for spidroin BackboneGenerator"
            sys.exit()        

    def generateBuildingBlock(self, minDist, species='SP1', showBlockDirector=False):
    
        self.species = species
        
        if species=='SP1':
            self.numGUnits = self.SP1NumGUnits 
            self.numPQUnits = self.SP1NumGUnits
        else:
            self.numGUnits = self.SP2NumGUnits 
            self.numPQUnits = self.SP2NumGUnits + 1
    
        self.numPoints = self.numGUnits * 3 + self.numPQUnits * 3  
        self.directorHat = np.array([0.0, 0.0, 1.0])
        self.minDist = minDist
        
        return BBG.generateBuildingBlock(self, self.numPoints, minDist, showBlockDirector=showBlockDirector)

    def generateBuildingBlockXYZ(self):
        # fundamental unit is a G-PQ unit.

        # create a start unit. In SP1 we start with a G-Q, in SP2 we start with a P-G-P 
        strand = self.startGPQUnit()

        # tally the number of GUnits we have added to strand        
        numGUnits = 1

        # loop until we've added the right number of GUnits
        while numGUnits < self.numGUnits:
            
            # create new G-PQ from the last G-PQ Unit (6 atoms)
            G_PQUnit = self.generateGPQUnit(strand[-6:])
            
            # add the new Residues to the strandVec list
            [strand.append(pos) for pos in G_PQUnit]
            
            numGUnits += 1
        
        # perform final orientation depending on whether or not a seed residue was provided.
        # If seed was provided then don't fiddle with orientation at all.
        # if seed wasn't provided then return the strand with it's centre of mass
        # at the origin and it's helical axis pointing along the z-axis
        if self.numPoints > 9:
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
    
    def startGPQUnit(self):
        # generates the first G-PQ unit of the spidroin Backbone
        
        # A G or PQ unit has three nodes: the NPos, the FPos and the CPos.
        # The NPos is the N terminal of the unit, The FPos is the force centre and the CPos is the C terminal.
        # Together these three things define an ellipsoid and a central repulsive/attractive field.
        
        if self.species == 'SP2':
            # First Unit in SP2 is a P-Unit
            # set up a dummy pos to sort out first dihedral
            DPos = np.array([1.0, 0.0, 0.0])
            
            # set first NPos to be at origin.
            NPosP1 = np.array([0.0, 0.0, 0.0])
          
            # Assume FPos is along director axis from NPos - bad assumption - we'll correct it later
            FPosP1 = NPosP1 + self.PPbondLength/2.0 * self.directorHat
            
            # Make the CPos at the end of the GUnit
            CPosP1 = NPosP1 + self.PPbondLength/2.0 * self.directorHat
    
            # set up the TNB frame to define the alpha and betas for the first part of the P to G Unit
            TNB1 =  coords.constructTNBFrame(DPos, NPosP1, CPosP1)
            NPosG = CPosP1 + self.CNbondLength *  coords.generateTNBVecXYZ(TNB1, self.angleC, self.PGPsi)

            # set up the TNB frame to define the alpha and betas for the second part of the P to G Unit
            TNB2 =  coords.constructTNBFrame(NPosP1, CPosP1, NPosG)

            # set up FPosG and CPosG
            dirHat = coords.generateTNBVecXYZ(TNB2, self.angleN, self.PGPhi)
            FPosG = NPosG + self.GG2bondLength/2.0 * dirHat  
            CPosG = FPosG + self.GG2bondLength/2.0 * dirHat

            # now generate NPosP2, FPosP2, CPosP2
            # set up the TNB frame to define the alpha and betas for the first part of the G to P Unit
            TNB3 =  coords.constructTNBFrame(CPosP1, NPosG, CPosG)
            NPosP2 = CPosG + self.CNbondLength *  coords.generateTNBVecXYZ(TNB3, self.angleC, self.GPPsi)

            # set up the TNB frame to define the alpha and betas for the second part of the G to P Unit
            TNB2 =  coords.constructTNBFrame(NPosG, CPosG, NPosP2)
            dirHat = coords.generateTNBVecXYZ(TNB3, self.angleN, self.GPPhi)
            
            # Assume FPos is along director axis from NPos - bad assumption - we'll correct it later
            FPosP2 = NPosP2 + self.PPbondLength/2.0 * dirHat 
            
            # Make the CPos at the end of the GUnit
            CPosP2 = NPosP1 + self.PPbondLength/2.0 * dirHat
        
            GPQUnit = [NPosP1, FPosP1, CPosP1, NPosG, FPosG, CPosG, NPosP2, FPosP2, CPosP2]
        
        # In SP1 we start with a G-Q, in SP2 we start with a P-G-P 
        if self.species == 'SP1':
            # First unit in SP1 is a G-Q unit.
            
            # set up a dummy pos to sort out first dihedral
            DPos = np.array([1.0, 0.0, 0.0])
            
            # set first NPos to be at origin.
            NPosG = np.array([0.0, 0.0, 0.0])
          
            # Assume FPos is along director axis from NPos - bad assumption - we'll correct it later
            FPosG = NPosG + self.GG1bondLength/2.0 * self.directorHat
            
            # Make the CPos at the end of the GUnit
            CPosG = FPosG + self.GG1bondLength/2.0 * self.directorHat

            # set up the TNB frame to define the alpha and betas for the first part of the G to P Unit
            TNB1 =  coords.constructTNBFrame(DPos, NPosG, CPosG)
            NPosQ = CPosG + self.CNBondLength *  coords.generateTNBVecXYZ(TNB1, self.CAngle, self.GQPsi)

            # set up the TNB frame to define the alpha and betas for the second part of the G to Q Unit
            TNB2 =  coords.constructTNBFrame(NPosG, CPosG, NPosQ)

            # set up FPosPQ and CPosPQ
            dirHat = coords.generateTNBVecXYZ(TNB2, self.NAngle, self.GQPhi)
            FPosQ = NPosQ + self.QQBondLength/2.0 * dirHat  
            CPosQ = FPosQ + self.QQBondLength/2.0 * dirHat
    
            GPQUnit = [NPosG, FPosG, CPosG, NPosQ, FPosQ, CPosQ]
    
        return GPQUnit
                
    def generateGPQUnit(self, prevGPQ):
        
        # given a list of 6 xyz positions for the previous G-PQ Unit
        # construct the next G-PQ unit with appropriate bond angles
        # and dihedral angles.
        
        # construct the first G Unit from the previous G-PQ unit
        TNB1 = coords.constructTNBFrame(prevGPQ[2], prevGPQ[3], prevGPQ[5])
        NPosG = prevGPQ[5] + self.CNBondLength * coords.generateTNBVecXYZ(TNB1, self.CAngle, self.PGPsi)
        
        TNB2 = coords.constructTNBFrame(prevGPQ[3], prevGPQ[5], NPosG)
        dirHat = coords.generateTNBVecXYZ(TNB2, self.NAngle, self.PGPhi)
        if self.species=='SP1':
            FPosG = NPosG + self.GG1bondLength/2.0 * dirHat
            CPosG = FPosG + self.GG1bondLength/2.0 *  dirHat
        else:
            FPosG = NPosG + self.GG2bondLength/2.0 * dirHat
            CPosG = FPosG + self.GG2bondLength/2.0 * dirHat
        
        
        TNB3 = coords.constructTNBFrame(prevGPQ[5], NPosG, CPosG)
        
        # now construct the P or Q unit 
        if self.species=='SP1':
            dirHat = coords.generateTNBVecXYZ(TNB3, self.CAngle, self.GPPsi)
            NPosP =  CPosG + self.CNBondLength * dirHat 
            TNB4 = coords.constructTNBFrame(NPosG, CPosG, NPosP)
            dirHat = coords.generateTNBVecXYZ(TNB4, self.NAngle, self.GPPhi)
            FPosP = NPosP + self.PPbondLength/2.0 * dirHat
            CPosP = FPosP + self.PPbondLength/2.0 * dirHat
            G_PQUnit = [NPosG, FPosG, CPosG, NPosP, FPosP, CPosP]
        else:
            dirHat = coords.generateTNBVecXYZ(TNB3, self.CAngle, self.GQPsi)
            NPosQ =  CPosG + self.CNBondLength * dirHat 
            TNB4 = coords.constructTNBFrame(NPosG, CPosG, NPosQ)
            dirHat = coords.generateTNBVecXYZ(TNB4, self.NAngle, self.GQPhi)
            FPosQ = NPosQ + self.QQbondLength/2.0 * dirHat
            CPosQ = FPosQ + self.QQbondLength/2.0 * dirHat
            G_PQUnit = [NPosG, FPosG, CPosG, NPosQ, FPosQ, CPosQ]
            
        return G_PQUnit

    def generateBuildingBlockNames(self):
        if self.species=='SP1':
            names = ['N', 'B', 'C', 'N', 'P', 'C'] * self.numGUnits
        else:
            names = ['N', 'P', 'C'] + ['N', 'B', 'C', 'N', 'P', 'C'] * self.numGUnits
        return names
    
    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)
    
    def generateBuildingBlockDirector(self):
        if self.numPoints > 9:
            director = coords.axisFromHelix(self.buildingBlockXYZ)
        else:
            director = self.buildingBlockXYZ[-1] - self.buildingBlockXYZ[0]
        directorHat = director/np.linalg.norm(director)   
        return directorHat 

    def generateBuildingBlockConnectors(self):
        return [np.array([2, 1, 0]), np.array([self.numPoints - 3, self.numPoints - 2, self.numPoints - 1])]
                
if __name__ == "__main__":
    
    filename = "spidroinBackbone.txt"

    # create the backbone generator object using static file parameters
    backboneObject = spidroinBackboneGenerator(filename)

    # generate backbone realtime parameters
    startPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 0 * np.pi/180
    
    SP1BuildingBlock = backboneObject.generateBuildingBlock('SP1')
    SP1BuildingBlock.transformBBToLabFrame(director, startPos, rotation)
    SP1BuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))

    SP2BuildingBlock = backboneObject.generateBuildingBlock('SP2')
    SP2BuildingBlock.transformBBToLabFrame(director, startPos, rotation)
    SP2BuildingBlock.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))


    print "spidroin backbone done"