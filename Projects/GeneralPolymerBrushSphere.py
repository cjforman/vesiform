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

    # break the list of XYZ points into list of points for each segment
    backBoneXYZList = breakListToSubList(backboneXYZ, np.cumsum(np.concatenate((np.array([0]), backboneSegmentLengths), 0)))

    # set up output
    polymerXYZ = cp.copy(backboneXYZ)
    polymerNames = cp.copy(backboneNames)
    
    for brush, numBrushes, TNB, points in zip(brushDictList, backboneSegmentLengths, backboneTNBList, backBoneXYZList):
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
            for point, director in zip(points, brushDirectorsNorm):
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
    brushDict2['mode'] = 'Polymer'
    brushDict2['name'] ='brush2'
    brushDict2['numMonomers'] = 7 
    brushDict2['monomerNames'] = ['Zn', 'O', 'Pb', 'N', 'B', 'H', 'S']
    brushDict2['phaseRange'] = 360.0 
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
       
    from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
    SphereBBG = SPSBBG('SurfacePackSphere.txt')

    numPolymersPerSphere =  185 #430
    UnimerBaseRadius = 2
    SphereRadius = 0.15 * (backboneDict1['Z2'] + backboneDict2['Z2'] - backboneDict1['Z1'] - backboneDict2['Z1'])
    phiMin = -90
    phiMax = 180    
    thetaMin = -90
    thetaMax = 90
    
    centerPos = np.array([0.0, 0.0, 0.0])

    # generate the XYZVals packed in the outer cylinder
    Polymer1SphereBB = SphereBBG.generateBuildingBlock(numPolymersPerSphere, SphereRadius, thetaMin, thetaMax, phiMin, phiMax, UnimerBaseRadius)
    Polymer1SphereBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
    Polymer1SphereBB.exportBBK("sphereBasePoints")
    Polymer1SpherePoints = Polymer1SphereBB.blockXYZVals 

    # find the radial director vector at each point
    directorsHat = [ director/np.linalg.norm(director) for director in Polymer1SpherePoints]

    # generate unique polymer strands
    for i in range(numPolymersPerSphere):
        print(i, "unimers out of: ", numPolymersPerSphere)
        if i==0:
            UnimerStrands=[GeneralBlockPolymer(polymerBrushDict)]
        else:
            UnimerStrands.append(GeneralBlockPolymer(polymerBrushDict))

    # transform the strands to the sphere points
    curStrand = 0
    for directorHat, pos, strand in zip(directorsHat, Polymer1SpherePoints, UnimerStrands):
        labRotation = rnd.uniform(0, 2 * np.pi)
        xyzVals = coords.transformFromBlockFrameToLabFrame(directorHat, pos, labRotation, np.array([0.0, 0.0, 1.0]), strand[0][0], strand[0])
        if curStrand==0:
            xyzValsList = xyzVals
            allNames = strand[1]
        else:
            xyzValsList = np.concatenate((xyzValsList, xyzVals), 0)
            allNames = np.concatenate( (allNames, strand[1]), 0)
        curStrand += 1

    fIO.saveXYZList(xyzValsList, allNames, "polymerSphere.xyz")

    print("example done")