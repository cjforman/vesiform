import numpy as np
import copy as cp
import random as rnd
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
import Utilities.cartesian as carts 
import Utilities.coordSystems as coords  
import Utilities.fileIO as fIO

def GeneralBlockPolymer(polymerDict):

    # unpack the brush dictionary
    connectorDictList = polymerDict['connectors']
    backBoneDictList = polymerDict['backbones']
    brushDictList = polymerDict['brushes']

    # initialise output
    polymerBrushXYZ = []
    polymerBrushNames = []
    
    if ( len(backBoneDictList)==len(brushDictList)==len(connectorDictList)+1)==False:
        print("polymer Component Dictionaries inconsistent") 
        return polymerBrushXYZ, polymerBrushNames
    
    # loop throught the backbone dictionary and generate the backbone of the generalised Block Co-Polymer
    backboneInfo = makeBlockCopolymerBackbone(backBoneDictList, connectorDictList)
    
    # decorate the polymer backbone with brushes using the brush information
    polymerBrushXYZ, polymerBrushNames, brushesPerBlock = decoratePolymerBackbone(brushDictList, backboneInfo)

    return polymerBrushXYZ, polymerBrushNames, brushesPerBlock
        
def makeBlockCopolymerBackbone(backboneDictList, connectorDictList):
    ''' function generates a block Copolymer with the parameters of each block defined in
        each entry of a list of parameter dictionaries.  The segments are then joined
        together in the order of the list according to the parameters in the list of connector
        dictionaries. '''    
    backboneSet = makeBackbones(backboneDictList)
    return connectBackboneSegments(backboneSet, connectorDictList)
    
def makeBackbones(backboneDictList):
    ''' Takes the backbone dictionary and makes a list of backbone building blocks defined in their own block frames ''' 
    
    backboneSet=[]
    
    for backbone in backboneDictList:
    
        # create a backbone generator object for the polymer
        backboneBBG = RPPBBG(backbone['filename'])
        
        # create the starting point of the polymer
        try:
            startPointBackBone = backbone['pointToStart']
        except KeyError:
            startPointBackBone = np.array([ 0.0, 0.0, backbone['Z1']])
    
        # generate the frustrum containing the polymer
        envelopeList = ['frustum ' + str(backbone['Z1']) + ' ' + str(backbone['R1']) + ' ' + str(backbone['Z2']) + ' ' + str(backbone['R2'])]
    
        try:
            pointsToAvoid = backbone['PointsToAvoid']
        except KeyError:
            pointsToAvoid = []
    
        # assume failure
        rejectPolymer = True

        # loop until we get a viable polymer    
        while rejectPolymer:
            # generate polymerBB of appropriate length. 
            backboneBB = backboneBBG.generateBuildingBlock( backbone['numMonomers'],  
                                                            startPointBackBone,
                                                            backbone['alpha1'],
                                                            backbone['alpha2'],
                                                            backbone['beta1'],
                                                            backbone['beta2'],
                                                            backbone['minDist'],
                                                            backbone['bondLength'],
                                                            envelopeList=envelopeList,
                                                            visualiseEnvelope=(0, 2*backbone['R2'] , backbone['name'] + '_envelope.xyz'),
                                                            pointsToAvoid=pointsToAvoid)
            
            backboneBB.blockAtomNames = backbone['monomerNames']
        
            # check that Z coord of point 0 is below z last point of polymer  otherwise reject polymer and choose another one.  
            if backboneBB.blockXYZVals[0][2] < backboneBB.blockXYZVals[-1][2]:
                rejectPolymer = False
            else:
                print("Polymer Rejected. Doing another.")
    
            # add building block to list of block copolymer segments
            if not backboneSet:
                backboneSet = [cp.copy(backboneBB)]
            else:
                backboneSet.append(cp.copy(backboneBB))

    return backboneSet

def connectBackboneSegments(backboneSet, connectorDictList):    
    ''' connects together a set of polymer building block objects according to the
    connection information in the connectorDictList. The order of the segments is defined by the order 
    they appear in the backboneSet list. 
    If there are N segments there must be N-1 connector objects with the link information.
    The connector information contains the angle between polymer directors, and the azimuthal rotation between
    the directors of each segment, as well as the distance from the last monomer of the first polymer
    and the first monomer of the second polymer, in the direction of the first polymer's 
    labdirector. The first two segments are always contained in the zx plane, and the first polymer is
    always aligned with the z axis. Setting all the azimuthal rotations to zero thus keeps all subsequent 
    segments in the zx plane which is useful for 2D figures, and for defining a zero plane
    to put side brushes in later on.'''
    
    # intialise the first labdirector, TNB frame and start point
    labDirector = np.array([0.0, 0.0, 1.0])
    labBidirector  = np.array([1.0, 0.0, 0.0])  # initial bi-director along the x axis.
    labNormDirector = np.cross(labDirector, labBidirector)
    zeroPoint = np.array([0.0, 0.0, 0.0])
    
    # initialise the output arrays with data from the first segment. 
    backboneXYZ = backboneSet[0].blockXYZVals
    backboneNames = backboneSet[0].blockAtomNames
    backboneIndexList = [len(backboneXYZ)]
    backboneTNBList = [ (cp.copy(labDirector), cp.copy(labNormDirector), cp.copy(labBidirector) ) ]  
    
    # find the principal axis of the first segment
    blockDirector = coords.axisFromHelix(backboneXYZ)
    
    # rotate the xyz vals to the z-axis lab director and translate the first point to the zeroPoint as defined
    backboneXYZ = coords.transformFromBlockFrameToLabFrame(labDirector, zeroPoint, 0.0, blockDirector, backboneXYZ[0], backboneXYZ)     
     
    # loop through the connectors and remaining segments orienting them according to the instructions in the connector list
    for connector, newBackboneSegmentBB in zip(connectorDictList, backboneSet[1:]):
        
        # generate unit TNB frame for the last segment 
        TNB = [labDirector, labNormDirector, labBidirector] 
        
        # compute the lab director for the new block and connector information, which is defined relative to end of last polymer.
        newLabDirector = coords.generateTNBVecXYZ(TNB, connector['beta']*np.pi/180.0, connector['alpha']*np.pi/180.0)

        # compute the zero point for the new block to attach to the end of the polymer backbone as it is just now
        zeroPoint = backboneXYZ[-1] + connector['displacement'] * newLabDirector 

        # find the current principal axis of the new segment
        blockDirector = coords.axisFromHelix(newBackboneSegmentBB.blockXYZVals)
               
        # rotate the xyz vals to the z-axis lab director and translate the first point to the zeroPoint as defined
        newBackboneXYZ = coords.transformFromBlockFrameToLabFrame(newLabDirector, 
                                                                  zeroPoint, 
                                                                  0.0, 
                                                                  blockDirector, 
                                                                  newBackboneSegmentBB.blockXYZVals[0], 
                                                                  newBackboneSegmentBB.blockXYZVals)     
        # generate new TNB frame for the next segment
        labNormDirector = np.cross(labDirector, newLabDirector)
        labBidirector = np.cross(labNormDirector, newLabDirector)
        labDirector = newLabDirector

        # add the output information to the output arrays
        backboneXYZ = np.concatenate((backboneXYZ, newBackboneXYZ), 0)
        backboneNames = np.concatenate((backboneNames, newBackboneSegmentBB.blockAtomNames), 0)
        backboneIndexList.append(len(newBackboneXYZ))
        backboneTNBList.append((cp.copy(newLabDirector), cp.copy(labNormDirector), cp.copy(labBidirector)))
        
    # return all the information that describes the backbone necessary for decorating it with side chains        
    return (backboneXYZ, backboneNames, backboneIndexList, backboneTNBList)

def breakListToSubList(inpList, indexPoints):
        return [ inpList[index:nextIndex] for index, nextIndex in zip(indexPoints[0:-1],indexPoints[1:])] 

def rotateLastNEntries(XYZList, index, angle):
    newXYZList = cp.copy(XYZList)
    p1 = newXYZList[index] - newXYZList[index - 1]
    p2 = newXYZList[index + 1] - newXYZList[index]
    n = np.cross(p2, p1)
    n = n/np.linalg.norm(n)
    for i,p in enumerate(newXYZList):
        if i>=index:
            newXYZList[i] = carts.rotPAboutAxis(p - newXYZList[index], n, angle) + newXYZList[index] 
    return newXYZList
            
def decoratePolymerBackbone(brushDictList, backboneInfo):
    # each dictionary defines the nature of the side chains on one block of the main polymer.
    # the overall behaviour is controlled by use of the mode keyword. 
    # THe mode Keyword can contain any combination of the following control phrases in any order
    #
    # Peptide or Polymer  (One of these is mandatory and you cannot have both).
    # Alternate   (the side chains are on alternate sides of the peptide (adds 180 to the phase)
    # Mirror (The phase of half of the side chains are randomly flipped by 180 degrees)
    # Random (side chains are dispersed statistically over the backbone)
    # Residues (filters a peptide for CA positions only to represent each residue).
    #
    # other keywords may be necessary and are defined in the brush dictionary separately, and not in the 'mode' keyword.     

    # unpack input
    backboneXYZ = backboneInfo[0]
    backboneNames = backboneInfo[1]
    backboneSegmentLengths = backboneInfo[2]
    backboneTNBList = backboneInfo[3]

    # break the list of XYZ points into sublists for the backbone points that belong to each block
    backBoneXYZList = breakListToSubList(backboneXYZ, np.cumsum(np.concatenate((np.array([0]), backboneSegmentLengths), 0)))

    # set up output - it is important to ensure the correlation between name and xyz points is always maintained. 
    # Should probably set them up as tuples but meh.
    # begin the final output which is just a list of points and point names for the xyz file 
    polymerXYZ = cp.copy(backboneXYZ)
    polymerNames = cp.copy(backboneNames)
    polymerBrushesPerBlock = []

    # loop through each block and add the side chain info to the output lists.     
    for brush, blockLength, TNB, points in zip(brushDictList, backboneSegmentLengths, backboneTNBList, backBoneXYZList):
        
        # make a 0-based list of backbone indices from indexOfFirstBrush to N - 1 where N is the blockLength.
        backboneIndexList = np.arange(brush['indexOfFirstBrush'], blockLength)
        
        # make a sublist starting at 0 to the end of the backboneIndex list with an increment of spacing
        brushPointIndexList = backboneIndexList[ 0 :: int( brush['spacing'] + 1 ) ]
        
        # count the number of brushes
        numBrushes = len(brushPointIndexList)

        # record number of brushes        
        polymerBrushesPerBlock.append(numBrushes)
        
        # now we now the number of brushes from the spacing information for this brush,
        # if we need them randomly distributed we can just choose the right number of indices randomly from the original list. 
        # Can't choose an index before indexOfFirstBrush because they aren't in the list!
        if 'Random' in brush['mode']:
            brushPointIndexList = np.sort(np.random.choice(backboneIndexList, size=numBrushes, replace=False))

        # only add a brush if the numMonomers is >0
        if brush['numMonomers']>0:
            # if the brush is a peptide sequence then it will be the same for all the brushes,
            # i.e. an alpha helix or a beta strand. For a hairpin, generate a random polymer. 
            if 'Peptide' in brush['mode']:
                # numMonomers in brush dict refers to number of residues!
                brushBBG = PBG(brush['filename']) # generate a betastrand or an alphahelix (if you want hair pin, just generate a random polymer)
                brushBB = brushBBG.generateBuildingBlock(brush['numMonomers'], showBlockDirector=False, nameCA=True)
                
                # if requested filter the peptide for only CA positions to represent a residue rather than backbone atoms. 
                if 'Residues' in brush['mode']:
                    brushNames = [ name for name in brushBB.blockAtomNames if name=='CA' ]
                    brushXYZ = [ pos for pos, name in zip(brushBB.blockXYZVals, brushBB.blockAtomNames) if name=='CA' ]
                else:
                    brushNames = brushBB.blockAtomNames
                    brushXYZ = brushBB.blockXYZVals
                try:
                    if brush['doProlineKink']:
                        prolinePoints= [i for i, x in enumerate(brush['monomerNames']) if x == "P"]
                        for prolinePoint in prolinePoints:
                            brushXYZ = rotateLastNEntries(brushXYZ, prolinePoint, brush['prolineKinkAngle'])
                except KeyError:
                    pass
            else:
                # otherwise we are using a polymer. Create a generator. Each will be different. 
                brushBBG = RPPBBG(brush['filename'])

            # Set the angle of each brush around the block director defined relative to the Binormal in the TNB frame for that polymer.
            # The range of angles is defined in the dictionary.
            brushPhaseAngles = [ rnd.uniform(0, brush['phaseRange'] * np.pi/180.0) for _ in range(numBrushes) ]

            # if the word Alternate is in the mode than add 180 to every other phase        
            if "Alternate" in brush['mode']:
                for brushNumber in range(numBrushes):
                    if brushNumber % 2==0:
                        brushPhaseAngles[brushNumber] += np.pi
            
            # if the word "mirror" is in the mode then add 180 randomly to half the brush phases 
            if "Mirror" in brush['mode']:
                for index in np.random.choice(np.arange(numBrushes), int(np.floor(numBrushes/2)), replace=False):
                    brushPhaseAngles[index] += np.pi
                    
            # generate directors in direction of phase angles
            brushDirectors = [ np.cos(angle) * TNB[2] + np.sin(angle) * TNB[1] for angle in brushPhaseAngles] 
            brushDirectorsNorm = [ d/np.linalg.norm(d) for d in brushDirectors]
    
            # for each of the points in the backbone for which we are generating a brush, create a brush 
            # and translate and rotate it to the correct place. For peptides we already generated a base brush and names, which we will clone.     
            for brushIndex, backBoneIndex in enumerate(brushPointIndexList):

                if "Polymer" in brush['mode'] and not "Peptide" in brush['mode']:
                    
                    # if we're doing a polymer then generate a new polymer with the given polymer parameters for each pass.
                    polyStart = np.array([ 0.0, 0.0, brush['Z1']])
                    envelopeList = ['frustum ' + str(brush['Z1'] - brush['minDist']) + ' ' + str(brush['R1']) + ' ' + str(brush['Z2']) + ' ' + str(brush['R2'])]
                
                    brushBB = brushBBG.generateBuildingBlock( brush['numMonomers'],
                                                              polyStart,
                                                              brush['alpha1'],
                                                              brush['alpha2'],
                                                              brush['beta1'],
                                                              brush['beta2'],
                                                              brush['minDist'],
                                                              brush['bondLength'],
                                                              envelopeList=envelopeList,
                                                              visualiseEnvelope=(0, 100, brush['name'] + '_envelope.xyz'))
                    
                    brushXYZ = brushBB.blockXYZVals
                    brushNames = brushBB.blockAtomNames
                        
                # overrides the names created in the building block generator with monomerNames in BrushDictionary
                try:
                    if len(brushXYZ)==len(brush['monomerNames']):
                        brushNames = brush['monomerNames']
                    else:
                        print("Warning: Number of XYZ vals and number of specified monomer names inconsistent. Not overriding generator.") 
                except KeyError:
                    pass

                # work out the brush director from the brushXYZ vals
                if brush['numMonomers']>2:
                    brushDirector = coords.axisFromHelix(brushXYZ)
                elif brush['numMonomers']==2:
                    brushDirector = brushXYZ[1] - brushXYZ[0]
                    brushDirector = brushDirector/np.linalg.norm(brushDirector)
                else:
                    brushDirector = np.array([0.0, 0.0, 1.0])

                # flip the director in Peptide mode. Work around some bug to do with the way the peptide chain is constructed  
                if "Peptide" in brush['mode'] and "Residues" not in brush['mode']:
                    brushDirector *= -1 
    
                try:
                    if brush['doProlineKink']:
                        prolinePoints= [i for i, x in enumerate(brush['monomerNames']) if x == "P"]
                        for prolinePoint in prolinePoints:
                            brushXYZ = rotateLastNEntries(brushXYZ, prolinePoint, brush['prolineKinkAngle'])
                except KeyError:
                    print("Skipping kink")
    
                # get the target point and director
                targetDirector = brushDirectorsNorm[brushIndex]
                targetPoint = points[backBoneIndex]

                # rotate, translate the coords and concatenate the brush positions and brush names to the output arrays
                polymerXYZ = np.concatenate( (polymerXYZ, 
                                              coords.transformFromBlockFrameToLabFrame( targetDirector,
                                                                                        targetPoint + brush['bondLength'] * targetDirector,
                                                                                        0.0,
                                                                                        brushDirector,
                                                                                        brushXYZ[0],
                                                                                        brushXYZ)), 0)
                polymerNames = np.concatenate( (polymerNames, brushNames), 0)
        
    return polymerXYZ, polymerNames, polymerBrushesPerBlock

if __name__=="__main__":

    # set up the dictionaries
    backboneDict1 = {}
    backboneDict2 = {}
    connectorDict12 = {}
    brushDict1 = {}
    brushDict2 = {}
    
    # pack the dictionaries
    backboneDict1['filename'] = 'RandomPolymer.txt'  
    backboneDict1['mode'] = 'Polymer'
    backboneDict1['name'] = 'BlockA'
    backboneDict1['numMonomers'] = 18 
    backboneDict1['monomerNames'] = backboneDict1['numMonomers'] * ['P']
    backboneDict1['alpha1'] = 40
    backboneDict1['alpha2'] = 50
    backboneDict1['beta1'] = 165
    backboneDict1['beta2'] = 185
    backboneDict1['minDist'] = 2.0
    backboneDict1['bondLength'] = 2.0
    backboneDict1['Z1'] = 2
    backboneDict1['R1'] = 20
    backboneDict1['Z2'] = 1.5 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
    backboneDict1['R2'] = 70

    backboneDict2['filename'] = 'RandomPolymer.txt'  
    backboneDict2['mode'] = 'Polymer'
    backboneDict2['name'] = 'BlockB'
    backboneDict2['numMonomers'] = 6
    backboneDict2['monomerNames'] = backboneDict2['numMonomers'] * ['C']
    backboneDict2['alpha1'] = 40
    backboneDict2['alpha2'] = 50
    backboneDict2['beta1'] = 165
    backboneDict2['beta2'] = 185
    backboneDict2['minDist'] = 2.0
    backboneDict2['bondLength'] = 2.0
    backboneDict2['Z1'] = 2
    backboneDict2['R1'] = 20
    backboneDict2['Z2'] = 1.5 * backboneDict2['numMonomers'] * backboneDict2['bondLength'] + backboneDict2['Z1']
    backboneDict2['R2'] = 70

    connectorDict12['displacement'] =  1.0 * backboneDict1['bondLength']
    connectorDict12['alpha'] = 0.0
    connectorDict12['beta'] = 175.0

    brushDict1['filename'] = 'RandomPolymer.txt'  
    brushDict1['mode'] = 'PolymerAlternate'
    brushDict1['name'] ='brush1'
    brushDict1['numMonomers'] = 1 
    brushDict1['monomerNames'] = brushDict1['numMonomers'] * ['Fe']
    brushDict1['phaseRange'] = 5.0 
    brushDict1['alpha1'] = 40
    brushDict1['alpha2'] = 50
    brushDict1['beta1'] = 145
    brushDict1['beta2'] = 155
    brushDict1['minDist'] = 1.0
    brushDict1['bondLength'] = 2.5
    brushDict1['Z1'] = 2.0
    brushDict1['R1'] = 10.0
    brushDict1['Z2'] = brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
    brushDict1['R2'] = 30.0

    brushDict2['filename'] = 'RandomPolymer.txt'  
    brushDict2['mode'] = 'PolymerAlternate'
    brushDict2['name'] ='brush2'
    brushDict2['numMonomers'] = 8 
    brushDict2['monomerNames'] = ['Na', 'Zn', 'O', 'Pb', 'N', 'B', 'H', 'S']
    brushDict2['phaseRange'] = 10.0 
    brushDict2['alpha1'] = 40
    brushDict2['alpha2'] = 50
    brushDict2['beta1'] = 155
    brushDict2['beta2'] = 175
    brushDict2['minDist'] = 1.0
    brushDict2['bondLength'] = 2.0
    brushDict2['Z1'] = 2
    brushDict2['R1'] = 50
    brushDict2['Z2'] = 1.5 * brushDict2['numMonomers'] * brushDict2['bondLength'] + brushDict2['Z1']
    brushDict2['R2'] = 50
              
    polymerBrushDict={}
    polymerBrushDict['backbones'] = [backboneDict1, backboneDict2]#, backboneDict1, backboneDict2]
    polymerBrushDict['brushes'] = [brushDict1, brushDict2]#, brushDict1, brushDict2]
    polymerBrushDict['connectors'] = [connectorDict12]#, connectorDict12, connectorDict12]
       
    strand = GeneralBlockPolymer(polymerBrushDict)
    fIO.saveXYZList(strand[0], strand[1], "Polymer.xyz")
    
    print("example done")