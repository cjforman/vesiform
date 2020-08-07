import numpy as np
import copy as cp
import random as rnd
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
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

    addBrushList = []    
    for backBone in backBoneDictList:
        addBrush = []
        for i in range(len(backboneInfo[0])):
            if i % backBone['spacing'] == backBone['spacingPhase']: 
                addBrush.append(1)
            else:
                addBrush.append(0)
    
        addBrushList.append(addBrush)
    
    backboneInfo = backboneInfo + (addBrushList,)
    
    # decorate the polymer backbone with brushes using the brush information
    polymerBrushXYZ, polymerBrushNames = decoratePolymerBackbone(brushDictList, backboneInfo)

    return polymerBrushXYZ, polymerBrushNames
        
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
        startPointBackBone = np.array([ 0.0, 0.0, backbone['Z1']])
    
        # generate the frustrum containing the polymer
        envelopeList = ['frustum ' + str(backbone['Z1']) + ' ' + str(backbone['R1']) + ' ' + str(backbone['Z2']) + ' ' + str(backbone['R2'])]
    
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
                                                            visualiseEnvelope=(0, 100, backbone['name'] + '_envelope.xyz'))
            
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
            
def decoratePolymerBackbone(brushDictList, backboneInfo): 

    # unpack input
    backboneXYZ = backboneInfo[0]
    backboneNames = backboneInfo[1]
    backboneSegmentLengths = backboneInfo[2]
    backboneTNBList = backboneInfo[3]
    addBrushList = backboneInfo[4]
    

    # break the list of XYZ points into list of points for each segment
    backBoneXYZList = breakListToSubList(backboneXYZ, np.cumsum(np.concatenate((np.array([0]), backboneSegmentLengths), 0)))

    # set up output
    polymerXYZ = cp.copy(backboneXYZ)
    polymerNames = cp.copy(backboneNames)
    
    for brush, numBrushes, TNB, points, addBrush in zip(brushDictList, backboneSegmentLengths, backboneTNBList, backBoneXYZList, addBrushList):
            if brush['numMonomers']>0:
                # create the brush generator. Default is polymer. Check to see if we want a peptide brush 
                if 'Peptide' in brush['mode']:
                    brushBBG = PBG(brush['filename'])
                    brushBB = brushBBG.generateBuildingBlock(brush['numMonomers']) # only need to create this once for peptide backbones as they are either alpha or beta.
                    brushBB.blockAtomNames = brush['monomerNames']
                else:
                    brushBBG = RPPBBG(brush['filename'])
            
                # Set the angle of each brush around the segment director defined relative to the Binormal in the TNB frame for that polymer.
                # The range of angles is defined in the dictionary.
                brushPhaseAngles = [ rnd.uniform(0, brush['phaseRange'] * np.pi/180.0) for _ in range(0, numBrushes) ]
    
                # if the word Alternate is in the mode than add 180 to every other phase        
                if "Alternate" in brush['mode']:
                    for index in range(0, numBrushes):
                        if index % 2 ==0:
                            brushPhaseAngles[index] += np.pi
                
                # if the word "mirror" is in the mode then add 180 randomly to half the brush phases 
                if "Mirror" in brush['mode']:
                    numFlipped = 0
                    indexList = list(range(0, numBrushes))
                    while numFlipped < numBrushes/2:
                        index = rnd.choice(indexList)
                        brushPhaseAngles[index] += np.pi
                        numFlipped += 1
                        indexList.remove(index)    
            
                # generate directors in direction of phase angles
                brushDirectors = [ np.cos(angle) * TNB[2] + np.sin(angle) * TNB[1] for angle in brushPhaseAngles] 
                brushDirectorsNorm = [ d/np.linalg.norm(d) for d in brushDirectors]
        
                # for each of the points in the backbone generate a brush as defined    
                for point, director, doBrush in zip(points, brushDirectorsNorm, addBrush):
                    if doBrush:
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
                            # work out the block director
                            if brush['numMonomers']>2:
                                blockDirector = coords.axisFromHelix(brushBB.blockXYZVals)
                            else:
                                if brush['numMonomers']==2:
                                    blockDirector = brushBB.blockXYZVals[1] - brushBB.blockXYZVals[0]
                                    blockDirector = blockDirector/np.linalg.norm(blockDirector)
                                if brush['numMonomers']==1:
                                    blockDirector = np.array([0.0, 0.0, 1.0])
            
                            # rotate the brush and set it at the right place
                            brushXYZVals = coords.transformFromBlockFrameToLabFrame(director, point + brush['bondLength'] * director, 0.0, blockDirector, brushBB.blockXYZVals[0], brushBB.blockXYZVals)
                    
                            # concatenate the brush positions and the brush names to the output arrays
                            polymerXYZ = np.concatenate( (polymerXYZ, brushXYZVals), 0)
                            polymerNames = np.concatenate( (polymerNames, brush['monomerNames']), 0)
        
    return polymerXYZ, polymerNames 

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
    backboneDict1['numMonomers'] = 21
    backboneDict1['monomerNames'] = backboneDict1['numMonomers'] * ['P']
    backboneDict1['spacing'] = 3
    backboneDict1['spacingPhase'] = 1
    backboneDict1['alpha1'] = 40
    backboneDict1['alpha2'] = 60
    backboneDict1['beta1'] = 140
    backboneDict1['beta2'] = 160
    backboneDict1['minDist'] = 1.0
    backboneDict1['bondLength'] = 2.0
    backboneDict1['Z1'] = 2
    backboneDict1['R1'] = 20
    backboneDict1['Z2'] = 15 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
    backboneDict1['R2'] = 240

    brushDict1['filename'] = 'RandomPolymer.txt'  
    brushDict1['mode'] = 'PolymerAlternate'
    brushDict1['name'] ='brush1'
    brushDict1['numMonomers'] = 10 
    brushDict1['monomerNames'] = ['Fe', 'C', 'O', 'S', 'Ca', 'K', 'Na', 'Zn', 'Mg', 'Cu']
    brushDict1['phaseRange'] = 40.0 
    brushDict1['alpha1'] = 30
    brushDict1['alpha2'] = 60
    brushDict1['beta1'] = 145
    brushDict1['beta2'] = 155
    brushDict1['minDist'] = 1.0
    brushDict1['bondLength'] = 2.0
    brushDict1['Z1'] = 2.0
    brushDict1['R1'] = 10.0
    brushDict1['Z2'] = brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
    brushDict1['R2'] = 20.0

                 
    polymerBrushDict={}
    polymerBrushDict['backbones'] = [backboneDict1]
    polymerBrushDict['brushes'] = [brushDict1]
    polymerBrushDict['connectors'] = []
       
    brush1 = GeneralBlockPolymer(polymerBrushDict)
    brush2 = GeneralBlockPolymer(polymerBrushDict)
    brush3 = GeneralBlockPolymer(polymerBrushDict)
    brush4 = GeneralBlockPolymer(polymerBrushDict)
    brush5 = GeneralBlockPolymer(polymerBrushDict)
    
    #fIO.saveXYZList(brush1[0], brush1[1], 'brush1.xyz')
    #fIO.saveXYZList(brush2[0], brush2[1], 'brush2.xyz')
    #fIO.saveXYZList(brush3[0], brush3[1], 'brush3.xyz')
    #fIO.saveXYZList(brush4[0], brush4[1], 'brush4.xyz')
    #fIO.saveXYZList(brush5[0], brush5[1], 'brush5.xyz')
    
    from Library.VolumePackEllipsoid import VolumePackEllipsoidBBG as VPEBBG
    SphereBBG = VPEBBG('VolumePackEllipsoid.txt')
    unimerBBG = RPPBBG('RandomPolymer.txt')
    
    numBrushes =  50 #185 #430
    numUnimers = 200
    brushRadius = 50
    unimerRadius = 20
    SphereRadius = 600 
    phiMin = -180
    phiMax = 180    
    thetaMin = -90
    thetaMax = 0
    
    centerPos = np.array([0.0, 0.0, 0.0])

    # generate the XYZVals packed in the outer cylinder
    BrushSphereBB = SphereBBG.generateBuildingBlock(numBrushes, SphereRadius, SphereRadius, SphereRadius, thetaMin, thetaMax, phiMin, phiMax, brushRadius)
    BrushSphereBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
    BrushSphereBB.exportBBK("sphereBasePoints")
    BrushSpherePoints = BrushSphereBB.blockXYZVals 

    # find a random director at each point
    brushAngles = [ coords.pickRandomPointOnUnitSphere() for director in BrushSpherePoints]
    brushDirectors = [ coords.sphericalPolar2XYZ( np.array([1.0, angle[0], angle[1]]) ) for angle in brushAngles ] 

    # generate unique polymer strands
    for i in range(numBrushes):
        print(i, "unimers out of: ", numBrushes)
        if i==0:
            Brushes=[GeneralBlockPolymer(polymerBrushDict)]
        else:
            Brushes.append(GeneralBlockPolymer(polymerBrushDict))

    # transform the brushes to the sphere points
    curStrand = 0
    for directorHat, pos, strand in zip(brushDirectors, BrushSpherePoints, Brushes):
        labRotation = rnd.uniform(0, 2 * np.pi)
        xyzVals = coords.transformFromBlockFrameToLabFrame(directorHat, pos, labRotation, np.array([0.0, 0.0, 1.0]), strand[0][0], strand[0])
        if curStrand==0:
            xyzValsList = xyzVals
            allNames = strand[1]
        else:
            xyzValsList = np.concatenate((xyzValsList, xyzVals), 0)
            allNames = np.concatenate( (allNames, strand[1]), 0)
        curStrand += 1

    # generate the XYZVals for unimers packed in the outer sphere
    UnimerSphereBB = SphereBBG.generateBuildingBlock(numUnimers, SphereRadius, SphereRadius, SphereRadius, thetaMin, thetaMax, phiMin, phiMax, unimerRadius)
    UnimerSphereBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
    UnimerSphereBB.exportBBK("sphereBasePoints")
    UnimerSpherePoints = UnimerSphereBB.blockXYZVals 

    # find a random director at each point
    unimerAngles = [ coords.pickRandomPointOnUnitSphere() for director in UnimerSpherePoints]
    unimerDirectors = [ coords.sphericalPolar2XYZ( np.array([1.0, angle[0], angle[1]]) ) for angle in unimerAngles ] 

    # generate unique unimer strands
    Unimers=[]
    for i in range(numUnimers):
        print(i, "unimers out of: ", numUnimers)
        envelopeList = ['frustum ' + str(brushDict1['Z1'] - brushDict1['minDist']) + ' ' + str(brushDict1['R1']) + ' ' + str(brushDict1['Z2']) + ' ' + str(brushDict1['R2'])]
        polyStart = np.array([ 0.0, 0.0, brushDict1['Z1']])
                        
        UnimerBB = unimerBBG.generateBuildingBlock( brushDict1['numMonomers'],
                                                    polyStart,
                                                    brushDict1['alpha1'],
                                                    brushDict1['alpha2'],
                                                    brushDict1['beta1'],
                                                    brushDict1['beta2'],
                                                    brushDict1['minDist'],
                                                    brushDict1['bondLength'],
                                                    envelopeList=envelopeList,
                                                    visualiseEnvelope=(0, 100, brushDict1['name'] + '_envelope.xyz'))

        fIO.saveXYZList(UnimerBB.blockXYZVals, brushDict1['monomerNames'], 'Unimer_'+str(i)+'.xyz')
            
        Unimers.append((UnimerBB.blockXYZVals, brushDict1['monomerNames']))

    for directorHat, pos, strand in zip(unimerDirectors, UnimerSpherePoints, Unimers):
        labRotation = rnd.uniform(0, 2 * np.pi)
        xyzVals = coords.transformFromBlockFrameToLabFrame(directorHat, pos, labRotation, np.array([0.0, 0.0, 1.0]), strand[0][0], strand[0])
        xyzValsList = np.concatenate((xyzValsList, xyzVals), 0)
        allNames = np.concatenate( (allNames, strand[1]), 0)
    
    fIO.saveXYZList(xyzValsList, allNames, "allVolume.xyz")

    print("example done")
    
