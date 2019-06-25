import numpy as np
#import networkx as nx
from Library.branchedPolymerNetwork2 import BranchedPolymerPackBBG as BPPBBG
from Library.PackNetwork import PackNetworkBBG as PNBBG
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
import random as rnd  

def makeWormNetwork(numGen1Worms, 
                    xSize, ySize, zSize, 
                    alpha1, alpha2, beta1, beta2, 
                    wormRadius, leafletThickness, backboneSegmentLength, unimerMinDist, packingFraction,
                    minSegmentsPerWorm, coneAngleWorm, wormDirector, wormAlignmentAngularRange, 
                    numWormGenerations, minBranchAngle, wormBranchDensity, maxNumAttempts, epsilon,
                    selfAvoid=False, pointsToAvoid=[], visualiseEnvelope=(0, 20, 'envelope.xyz'), envelopeList=["None"]):

    # create the work network generator
    wormNetworkBackboneBBG = BPPBBG('branchedPolymerNetwork.txt')
    packNetworkBBG = PNBBG('networkPack.txt')
        
    # generate a set of points and create a random polymer pack network in a region of space
    BranchedPolymerPackBB, networkGraph = wormNetworkBackboneBBG.generateBuildingBlock(  numGen1Worms, 
                                                                               -xSize/2, xSize/2,
                                                                               -ySize/2, ySize/2,
                                                                               -zSize/2, zSize/2, 
                                                                               2 * wormRadius,
                                                                               alpha1, alpha2,
                                                                               beta1, beta2,
                                                                               0.9 * backboneSegmentLength, 
                                                                               backboneSegmentLength,
                                                                               minSegmentsPerWorm,
                                                                               maxNumAttempts,
                                                                               selfAvoid,
                                                                               coneAngleWorm, 
                                                                               wormAlignmentAngularRange,
                                                                               wormDirector, 
                                                                               numWormGenerations, 
                                                                               wormBranchDensity,
                                                                               minBranchAngle, 
                                                                               visualiseEnvelope=visualiseEnvelope,
                                                                               pointsToAvoid=pointsToAvoid,
                                                                               envelopeList=envelopeList,
                                                                               SpaceCurveTransform=False)

    fIO.saveXYZ(BranchedPolymerPackBB.blockXYZVals, 'C', "wormBackBoneNetwork.xyz")
    
 
    # dump the network to file to help judge the output
    #fIO.visualiseNetwork(networkGraph, BranchedPolymerPackBB.blockXYZVals, 'visualisedNetwork.xyz', particleName='P', dZ=0.1)

    # guestimate number of points
    nEdges = len(networkGraph.edges)
    
    #numpoints = total surface area of network surface/ area of unimer. Just an approximationto total surface area.
    # can adjust packingRatio to get desired results. factor of 2pi cancels out. 
    numPoints = int(packingFraction * nEdges * wormRadius * backboneSegmentLength / unimerMinDist)
       
    # generate the networkpacked building block - get a set of points and director for each point
    (NetworkPackBB, directors) = packNetworkBBG.generateBuildingBlock(numPoints, unimerMinDist, 
                                                         networkGraph, BranchedPolymerPackBB.blockXYZVals, 
                                                         [wormRadius, leafletThickness], 'bilayer', epsilon)
    
    fIO.saveXYZ(NetworkPackBB.blockXYZVals, 'C', "networkPack.xyz")
    
    dirPosns = [ pos + director for (pos, director) in zip(NetworkPackBB.blockXYZVals, directors) ]
    
    fIO.saveXYZ(dirPosns, 'O', "networkPackDirectorPoints.xyz")
    
    return NetworkPackBB.blockXYZVals, directors

def computeSurfaceNormal(points):
    return [ np.array([0.0, 0.0, 1.0] for _ in points) ]

if __name__=="__main__":

    # set up the dictionaries describing the polymer brush
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
    backboneDict1['minDist'] = .4
    backboneDict1['bondLength'] = .4
    backboneDict1['Z1'] = .4
    backboneDict1['R1'] = 4
    backboneDict1['Z2'] = 1.5 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
    backboneDict1['R2'] = 14.0

    backboneDict2['filename'] = 'RandomPolymer.txt'  
    backboneDict2['mode'] = 'Polymer'
    backboneDict2['name'] = 'BlockB'
    backboneDict2['numMonomers'] = 6
    backboneDict2['monomerNames'] = backboneDict2['numMonomers'] * ['C']
    backboneDict2['alpha1'] = 40
    backboneDict2['alpha2'] = 50
    backboneDict2['beta1'] = 165
    backboneDict2['beta2'] = 185
    backboneDict2['minDist'] = 0.4
    backboneDict2['bondLength'] = 0.4
    backboneDict2['Z1'] = .4
    backboneDict2['R1'] = 4.0
    backboneDict2['Z2'] = 1.5 * backboneDict2['numMonomers'] * backboneDict2['bondLength'] + backboneDict2['Z1']
    backboneDict2['R2'] = 14.0

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
    brushDict1['minDist'] = .2
    brushDict1['bondLength'] = .5
    brushDict1['Z1'] = 0.4
    brushDict1['R1'] = 2.0
    brushDict1['Z2'] = brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
    brushDict1['R2'] = 6.0

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
    brushDict2['minDist'] = 0.2
    brushDict2['bondLength'] = 0.4
    brushDict2['Z1'] = 0.4
    brushDict2['R1'] = 10.0
    brushDict2['Z2'] = 1.5 * brushDict2['numMonomers'] * brushDict2['bondLength'] + brushDict2['Z1']
    brushDict2['R2'] = 10.0
              
    polymerBrushDict={}
    polymerBrushDict['backbones'] = [backboneDict1, backboneDict2]#, backboneDict1, backboneDict2]
    polymerBrushDict['brushes'] = [brushDict1, brushDict2]#, brushDict1, brushDict2]
    polymerBrushDict['connectors'] = [connectorDict12]#, connectorDict12, connectorDict12]
       
    numGen1Worms = 6
    xSize = 850
    ySize = 850
    zSize = 200
    alpha1 = -5 
    alpha2 = 5
    beta1 = 160
    beta2 = 180
    wormRadius = 10.0
    leafletThickness = 1.0
    segmentLength = 20
    unimerMinDist = .2
    packingFraction = 0.6
    minSegmentsPerWorm = 8
    coneAngleWorm = 45
    wormDirector = np.array([1.0, 0.0, 0.0])
    wormAlignmentAngularRange = 15
    numWormGenerations = 1
    minBranchAngle = 40
    wormBranchDensity = 0.08
    maxNumAttempts = 10
    epsilon = 1.0
    
    WormNetworkXYZPoints, WormNetworkXYZDirectors = makeWormNetwork( numGen1Worms, 
                                                                    xSize, ySize, zSize, 
                                                                    alpha1, alpha2, beta1, beta2, 
                                                                    wormRadius, leafletThickness, segmentLength, unimerMinDist, packingFraction,
                                                                    minSegmentsPerWorm, coneAngleWorm, wormDirector, wormAlignmentAngularRange, 
                                                                    numWormGenerations, minBranchAngle, wormBranchDensity, maxNumAttempts, epsilon,
                                                                    selfAvoid=False, pointsToAvoid=[], visualiseEnvelope=(0, 220, 'envelope.xyz'), envelopeList=["None"])
    # generate the XYZVals packed in the outer cylinder
#    fIO.saveXYZList(WormNetworkXYZPoints, ['C'] * len(WormNetworkXYZPoints), "wormNetwork.xyz")

    from Projects.GeneralPolymerBrush import GeneralBlockPolymer as GBCP
    # generate unique polymer strands
    UnimerStrands = []
    UnimerNum = 1
    for point, director in zip(WormNetworkXYZPoints, WormNetworkXYZDirectors):
        print(" generating unimer ", UnimerNum, " of ", len(WormNetworkXYZPoints))
        UnimerStrands += [ GBCP(polymerBrushDict)]
        UnimerNum += 1

    # transform the strands to the cylinder points
    xyzValsList = []
    allNames = []
    curStrand = 0
    for directorHat, pos, strand in zip(WormNetworkXYZDirectors, WormNetworkXYZPoints, UnimerStrands):
        labRotation = 0.0 # rnd.uniform(0, 2 * np.pi)
        xyzVals = coords.transformFromBlockFrameToLabFrame(directorHat, pos, labRotation, np.array([0.0, 0.0, 1.0]), strand[0][0], strand[0])
        xyzValsList += xyzVals
        allNames += strand[1].tolist()
        curStrand += 1

    # transform the strands to the inverse cylinder points
    curStrand = 0
    for directorHat, pos, strand in zip(WormNetworkXYZDirectors, WormNetworkXYZPoints, UnimerStrands):
        labRotation = 0.0 # rnd.uniform(0, 2 * np.pi)
        xyzVals = coords.transformFromBlockFrameToLabFrame(directorHat, pos, labRotation, np.array([0.0, 0.0, 1.0]), strand[0][0], strand[0])
        xyzValsList += xyzVals
        allNames += strand[1].tolist()
        curStrand += 1
        
        
 
    fIO.saveXYZList(xyzValsList, allNames, "WholeThing.xyz")

    print("example done")
