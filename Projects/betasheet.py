import sys
import numpy as np
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.peptideHairpin3 import peptideHairpinGenerator as PHG
from Utilities import cartesian as cart

class betasheetGen(BBG):
    # Class generates a beta sheet building block by constructing a polymer as a 
    # base line with a separation of betaStrandDistance and then distributing betaStrand 
    # objects using the seedBlock functionality
    
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        
        self.betaStrandSeparation = self.getParam('betaStrandSeparation')
        self.CCbondLength = self.getParam('CCbondLength')
        self.CNbondLength = self.getParam('CNbondLength')        
        self.phi = self.getParam('phi') * np.pi / 180.0
        self.psi = self.getParam('psi') * np.pi / 180.0
        self.omega = self.getParam('omega') * np.pi / 180.0
        self.angleN = self.getParam('angleN') * np.pi / 180.0
        self.angleCA = self.getParam('angleCA') * np.pi / 180.0
        self.angleC = self.getParam('angleC') * np.pi / 180.0

    def generateBuildingBlock(self, numStrands, strandLength, inStrandDirector, crossStrandDirector, offset, polarity='NC', parallel=True, loopedEnds=True):
    
        self.startPos = np.array([0.0, 0.0, 0.0])
        self.numStrands = numStrands
        self.numPoints = strandLength * 3 # number of points per strand.
        self.strandLength = strandLength
        self.inStrandDirectorHat = inStrandDirector/np.linalg.norm(inStrandDirector)
        self.crossStrandDirectorHat = crossStrandDirector/np.linalg.norm(crossStrandDirector)
        self.oStrandDirectorHat = np.cross(self.inStrandDirectorHat, self.crossStrandDirectorHat)
        self.offsetXYZ = offset[2] * self.inStrandDirectorHat + offset[0] * self.crossStrandDirectorHat + offset[1] * self.oStrandDirectorHat
        self.polarity = polarity
        self.parallel = parallel
        self.loopEnds = loopedEnds 
        
        return BBG.generateBuildingBlock(self, self.numPoints)
    
    def generateBuildingBlockXYZ(self):

        # create a strand and loop generator
        strandGen =PBG(self.paramFilename)
        loopGen = PHG(self.paramFilename)

        loopPoints =[]
        sphereCentrePoints = []
        dummyConnector = [[2, 1, 0]]

        # construct the base line for the strand start points
        baseLine = [ self.startPos + n * self.betaStrandSeparation * self.crossStrandDirectorHat for n in range(self.numStrands)]

        # construct a prototype beta strand building Block in the defined polarity
        strandA_BB = strandGen.generateBuildingBlock(self.numPoints, self.polarity, dummyConnector)
        
        # assume loop length needs to be long enough to stretch a whole beta strand.
        loopLength = int(np.ceil(1.5 * self.strandLength))
        
        # In a not parallel situation construct the antiparallel strand - StrandB
        # also add the offset to every other strand starting on the first or second strand depending on the polarity
        # also compute the alternative strand building block
        if not self.parallel:
            loopLength = 2 # in a parallel situation the loop is always 2 (residues) long. I just said so.
            if self.polarity=='NC':
                baseLine = [ basePoint + self.offsetXYZ if ( n % 2)==0 else basePoint for n, basePoint in enumerate(baseLine)]
                strandB_BB = strandGen.generateBuildingBlock(self.numPoints, 'CN', dummyConnector)
            else: 
                baseLine = [ basePoint + self.offsetXYZ if ( n % 2)==1 else basePoint for n, basePoint in enumerate(baseLine)]
                strandB_BB = strandGen.generateBuildingBlock(self.numPoints, 'NC', dummyConnector)
     
        # compute the relative vector from the end of the first chain to the next atom in the hairpin. 
        TNB = coords.constructTNBFrame(strandA_BB.xyzVals[-3], 
                                       strandA_BB.xyzVals[-2],
                                       strandA_BB.xyzVals[-1])

        if self.polarity == 'NC':
            # if polarity is NC then the last three values are an NCC triad.  
            beta =  self.angleC
            alpha = self.phi
            bondLength = self.CNbondLength
        else:
            # if polarity is CN then the last three values are CCN triad.  
            beta =  self.angleN
            alpha = self.psi
            bondLength = self.CNbondLength

        # the vector from the end of the first strand to its hairpin is defined as the positive connectorVector.
        ConnectorVectorPos = bondLength * coords.generateTNBVecXYZ(TNB, beta, alpha)
        
        # -ve connector vector is the reciprocal of that 
        ConnectorVectorNeg = -1 * ConnectorVectorPos

        # for parallel all CVs are always the Positive from end of strand to hairpin
        # and from hairpin to end of strand
        connectorVectorStrandToHairpin = ConnectorVectorPos
        connectorVectorHairpinToStrand = ConnectorVectorPos
        
        # for parallel strands are always the same set of xyz vals
        curStrandXYZ = strandA_BB.xyzVals

        # now the preparatory stuff is complete begin constructing the beta sheet
        # copy across the first strand to get the structure going.
        buildingBlockXYZ = [ baseLine[0] + xyz for xyz in strandA_BB.xyzVals]

        innerRadiusParallel = 0.75 * self.strandLength * 3 * bondLength/2
        outerRadiusParallel = 40 * self.strandLength * 3 * bondLength/2
        innerRadiusAntiParallel = 0.9 * self.betaStrandSeparation/2
        outerRadiusAntiParallel = 40 * self.betaStrandSeparation
        print "Radii: ", innerRadiusParallel, outerRadiusParallel, innerRadiusAntiParallel, outerRadiusAntiParallel
        print "Diameter: ", 2*innerRadiusParallel, 2*outerRadiusParallel, 2*innerRadiusAntiParallel, 2*outerRadiusAntiParallel
            
        # contruct the remaining strands following this procedure
        curStrand = 1
        while curStrand < self.numStrands:

            # if parallel, strands and connectorvectors are always the same so do nothing
            # to them.
            
            # if not parallel work out the connectorVectors for attaching hairpin between 
            # the previous strand and the current strand.
            if not self.parallel:
                if ( curStrand % 2 )==0:
                    connectorVectorStrandToHairpin = ConnectorVectorNeg
                    connectorVectorHairpinToStrand = ConnectorVectorPos
                    curStrandXYZ = strandA_BB.xyzVals
                if ( curStrand % 2)==1:
                    connectorVectorStrandToHairpin = ConnectorVectorPos
                    connectorVectorHairpinToStrand = ConnectorVectorNeg
                    curStrandXYZ = list(reversed(strandB_BB.xyzVals))

            # loop start is end of building block + connectorvector
            loopStartPoint = buildingBlockXYZ[-1] + connectorVectorStrandToHairpin
            
            loopPoints.append(loopStartPoint)
            
            # compute values for new strand
            newStrand = [ baseLine[curStrand] + xyz for xyz in curStrandXYZ]

            # loop end point is start of next strand - connector vector.
            loopEndPoint = newStrand[0] - connectorVectorHairpinToStrand
            
            loopPoints.append(loopEndPoint)
            sphereCentrePoint = (loopEndPoint+loopStartPoint)/2
            sphereCentrePoints.append(sphereCentrePoint)
            
            print "loopPointdist: ", np.linalg.norm(loopEndPoint - loopStartPoint)
            
            # now we have the length of the loop and the start and end points
            # makes and attempt to stop the hairpin from entering
            # a sphere at the mid point of the start and end points
            if self.loopEnds:
                if self.parallel:
                    
                    loop = loopGen.generateBuildingBlock( loopLength * 3, 
                                                          loopStartPoint,
                                                          loopEndPoint, 
                                                          0, 
                                                          -180, 
                                                          180, 
                                                          beta*180/np.pi - 40, 
                                                          beta*180/np.pi + 40, 
                                                          0.9*bondLength, 
                                                          bondLength, 
                                                          innerRadiusParallel, 
                                                          outerRadiusParallel, 
                                                          sphereCentrePoint,
                                                          self.polarity)
                else:
                    loop = loopGen.generateBuildingBlock( loopLength * 3, 
                                                          loopStartPoint,
                                                          loopEndPoint, 
                                                          0, 
                                                          -180, 
                                                          180, 
                                                          beta*180/np.pi - 40, 
                                                          beta*180/np.pi + 40, 
                                                          0.9*bondLength, 
                                                          bondLength,
                                                          innerRadiusAntiParallel, 
                                                          outerRadiusAntiParallel,
                                                          sphereCentrePoint, 
                                                          self.polarity)
                # append the loop to the array
                buildingBlockXYZ = buildingBlockXYZ + loop.xyzVals 
                
                 
            # append each vector to the current output array
            buildingBlockXYZ = buildingBlockXYZ + newStrand
        
            
        
            #next strand
            curStrand += 1
        
        #buildingBlockXYZ = buildingBlockXYZ + loopPoints
        #buildingBlockXYZ = buildingBlockXYZ + sphereCentrePoints
        
        self.numPoints = len(buildingBlockXYZ) 
        
        return buildingBlockXYZ 

    def generateBuildingBlockNames(self):
        #totalNumResidues = len(self.buildingBlockXYZ)/3
        totalNumResidues = self.strandLength * self.numStrands 
        if (self.polarity == "NC"):
            names = ['N', 'C', 'C'] * totalNumResidues
        if (self.polarity == "CN"):
            names = ['C', 'C', 'N'] * totalNumResidues

        #names = names + ['O'] * (numStrands - 1) * 2
        
        #names = names + ['H'] * (numStrands - 1)
        
        return names

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])

    def generateBuildingBlockDirector(self):
        return self.inStrandDirectorHat  
    
if __name__=="__main__":
    # get the file name from the command line
    filename = sys.argv[1]

    # create the backbone generator object using static file parameters
    betasheetGenerator = betasheetGen(filename)

    # generate backbone realtime parameters
    numStrands = 4
    lengthStrand = 7
    startPos = np.array([0.0, 0.0, 0.0])
    globalDirector = np.array([1.0, 1.0, 1.0])
    inStrandDirector = np.array([0.0, 0.0, 1.0])
    crossStrandDirector = np.array([1.0, 0.0, 0.0])
    rotation = 0 * np.pi/180
    offset =  np.array([1.5, 0.0, 1.5])
    polarity = 'NC'
    
    betasheetBB = betasheetGenerator.generateBuildingBlock(  numStrands, 
                                                             lengthStrand, 
                                                             inStrandDirector, 
                                                             crossStrandDirector, 
                                                             offset, 
                                                             polarity, 
                                                             parallel = False, 
                                                             loopedEnds = False) 

    betasheetBB.transformFromBlockFrameToLabFrame(globalDirector, startPos, rotation)
    betasheetBB.exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print "backbone done"
    
    
        
    