import numpy as np
import sys
import networkx as nx
from Library.SurfacePackCylinder import SurfacePackCylinderBBG as SPCBBG
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.branchedPolymerNetwork2 import BranchedPolymerPackBBG as BPPBBG
from Library.PackNetwork import PackNetworkBBG as PNBBG
import Utilities.fileIO as fIO
 

def makeWormNetwork(numGen1Worms, 
                    xSize, ySize, zSize, 
                    alpha1, alpha2, beta1, beta2, 
                    wormRadius, backboneSegmentLength, unimerMinDist, packingFraction,
                    minSegmentsPerWorm, coneAngleWorm, wormDirector, wormAlignmentAngularRange, 
                    numWormGenerations, minBranchAngle, wormBranchDensity, maxNumAttempts,
                    selfAvoid=False, pointsToAvoid=[], visualiseEnvelope=(0,20,'envelope.xyz'), envelopeList=["None"]):

    # create the work network generator
    wormNetworkBackboneBBG = BPPBBG('branchedPolymerNetwork.txt')
    packNetworkBBG = PNBBG('packnetwork.txt')
        
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
    
    # takes the graph of the network and the set of points to which the graph refers (the node names in the graph
    # are indices in the point list) and returns a set of xyzPoints and directors that populate the surface
    # of a cylindrical tube a distance worm radius from the backdone described by each individual cluster in the network.
    # the density of the xyzPoints are defined by unimerMinDist
    return populateClusters(networkGraph, BranchedPolymerPackBB.blockXYZVals, wormRadius, unimerMinDist, packingFraction)
    
def populateClusters(graph, points, radius, minDist, packingFraction):
    # analyse the graph and create a list of individual subgraphs which are connected.
    clusters = [graph.subgraph(c) for c in nx.connected_components(graph) ]
    
    # for each disjoint subgraph in the graph (a cluster) generate a series of unimer base 
    # points and directors across the surface of that cluster.
    clusterXYZ = []
    clusterDirectors = []
    for clusterNum, cluster in enumerate(clusters):
        print("Processing Cluster:", clusterNum + 1, " of ", len(clusters))
        # make cluster returns a tuple of points and directors that surround the cluster
        clusterInfo = makeCluster(cluster, points, radius, minDist, packingFraction)
        clusterXYZ += clusterInfo[0]
        clusterDirectors += clusterInfo[1]
        
    return clusterXYZ, clusterDirectors

# Generates an array of points, and surface normal directors for each point, which are on the surface 
# of a complex constructions of spheres and cylinders which form a connected worm like
# network. 
#
# The wormlike network is defined by a graph. The nodes in the graph refer to indices in the list of points 
# which are the centres of spheres of a given radius. The edges in the graph refer to cylinders of the same radius 
# connecting the nodes. 
#
# A list of random points and directors are generated for each sphere and cylinder. 
# Any points which are within the boundary defined by the union of the spheres and cylinders are deleted. 
# This is achieved by performing a pair wise analysis of each of the individual sphere and cylinders. 
# the points packed onto the surfaces are of uniform density defined by minDist.   
# All points on individual surfaces which are inside the composite surface are removed.
# 
# The Output is a flat list of points and their corresponding directors which define the wormlike network 
#
def makeCluster(cluster, points, radius, minDist, packingFraction):
    CylinderBBG = SPCBBG('SurfacePackCylinder.txt')
    SphereBBG = SPSBBG('SurfacePackSphere.txt')
    
    # determine the point density so it is the same on spheres and cylinders
    # number of points on Sphere = Area of sphere divided by area of circles of radius minDist.
    # np.pi's cancel out. Apply the packing fraction is the ratio relative to perfect packing
    nS = int(packingFraction * 4.0 * radius**2/ minDist**2)
    
    # compute number of points per unit length for each cylinder
    # area of cross section * 1 unit length / area of each circle
    # the np.pi's cancel out
    nCUnit_float = packingFraction * radius**2 / minDist**2

    # for each node in the graph generate a sphere of nS points of radius radius and density minDist 
    #spheres = [ SphereBBG.generateBuildingBlock(nS, radius, -90.0, 90.0, -180.0, 180.0, minDist) for _ in cluster ] 
    
    # compute an array of the sphere directors. Each sphere centred at origin so each pos is its own radial vector.
    # just normalise it and return in seperate array.
    #sphereDirectors = [ [ pos/np.linalg.norm(pos) for pos in sphere.blockXYZVals] for sphere in spheres ]

    # translate the sphere points so the centre of each sphere is at the given node.
    # Return an array the same shape as the director array for each point
    #sphereXYZPoints = [ [points[node] + pos for pos in sphere.blockXYZVals] for node, sphere in zip(cluster, spheres) ]
    
    # define the length of each cylinder in the cluster - one for each edge in the graph
    cylinderLengths = [ np.linalg.norm(points[edge[1]] - points[edge[0]]) for edge in cluster.edges ]

    # define the unit vector between each connected pairs of nodes
    cylinderAxes = [ (points[edge[1]] - points[edge[0]])/length for edge, length in zip(cluster.edges, cylinderLengths) ]

    # populate each cylinder with a specified number of points no less than minDist apart    
    cylinders = [ CylinderBBG.generateBuildingBlock(int(length * nCUnit_float), radius, radius, 0, length, -180.0, 180.0, minDist) for length in cylinderLengths ]
    
    # transform the cylinder points so the cylinder axis is aligned with the director and placed at the first node in the edge.  
    [ cylinderBB.transformBBToLabFrame(director, points[edge[0]], 0.0) for cylinderBB, director, edge in zip(cylinders, cylinderAxes, cluster.edges)]
    
    # extract the cylinder points from the building block object
    cylinderXYZPoints = [ cylinderBB.blockXYZVals for  cylinderBB in cylinders]
    
    # compute the director of each point (radial vector pointing away from the cylinder axis (rodrigues formula) (see page 183 in book 5).
    cylinderDirectors = [  [pos - points[edge[0]] - np.dot(director, pos - points[edge[0]]) * director for pos in cylinder] for cylinder, director, edge in zip(cylinderXYZPoints, cylinderAxes, cluster.edges) ]
    
    # normalise the directors
    cylinderDirectors = [  [pos/np.linalg.norm(pos) for pos in cylinder] for cylinder in cylinderDirectors ]
    
    # For each list of points that we have constructed (spheres and cylinders) 
    # go through and remove any of those points that are inside the *cylindrical* shells.
    # None of the spheres would ever have other sphere points or cylinder points inside them so no need to check them.
    # The cylinders dock around the spheres tangentially. 
    
    # first deal with the points belonging to the cylindrical surfaces
    
    # assume that all the points will survive - set up a flag array for each point that is same dims as director and points array
    cylinderPointFlagsList = [ len(cylinderPoints) * [True] for cylinderPoints in cylinderXYZPoints  ]
    
    # loop through the sets of points for each cylinder
    for indexInner, cylinderPoints in enumerate(cylinderXYZPoints):
        
        # loop through the geometric regions
        for director, length, indexOuter, edge in zip(cylinderAxes, cylinderLengths, range(0, len(cylinderLengths)), cluster.edges):
            # don't check a cylinder against itself
            if indexInner!=indexOuter:
                # test the current list of points against the cylinder
                # Set flags in the flag array for points that need to be removed   
                FlagPointsInsideCylinder(cylinderPointFlagsList[indexInner], cylinderPoints, points[edge[0]], director, length, radius)

    # apply the flags
    cylinderXYZPoints = [ [point for point, flag in zip(pointList, cylinderPointFlagsList[index]) if flag==True] for index, pointList in enumerate(cylinderXYZPoints)]
    cylinderDirectors = [ [director for director, flag in zip(directorList, cylinderPointFlagsList[index]) if flag==True] for index, directorList in enumerate(cylinderDirectors)]
    # reconstruct the flag list as all true
    cylinderPointFlagsList = [ len(cylinderPoints) * [True] for cylinderPoints in cylinderXYZPoints  ]
                
    # do the same for the sphere points
    
    # assume that all the points will survive - set up a flag array for each point that is same dims as director and points array
    #spherePointFlagsList = [ len(spherePoints) * [True] for spherePoints in sphereXYZPoints  ]
 
    # loop through the sets of points for each sphere
    #for indexInner, spherePoints in enumerate(sphereXYZPoints):
                
        # loop through the geometric regions
        #for director, length, edge in zip(cylinderAxes, cylinderLengths, cluster.edges):
            # test the current list of points against the cylinder 
            # set flags in the flag array for points that need to be removed   
           # FlagPointsInsideCylinder(spherePointFlagsList[indexInner], spherePoints, points[edge[0]], director, length, radius)
 
    # apply the flags
    #sphereXYZPoints = [ [point for point, flag in zip(pointList, spherePointFlagsList[index]) if flag==True] for index, pointList in enumerate(sphereXYZPoints)]
    #sphereDirectors = [ [director for director, flag in zip(directorList, spherePointFlagsList[index]) if flag==True] for index, directorList in enumerate(sphereDirectors)]
    # reconstruct the flag list as all true
   # spherePointFlagsList = [ len(spherePoints) * [True] for spherePoints in sphereXYZPoints  ]

    # now flag any points in the cylinders that are within minDist of the remaining points of the other cylinders. 
    # Don't check a cylinder against a cylinder that has already been checked.
    for indexTestSet, testPoints in enumerate(cylinderXYZPoints): 
        for refPoints in cylinderXYZPoints[indexTestSet + 1:]:
            FlagPointsTooClose(cylinderPointFlagsList[indexTestSet], testPoints, refPoints, minDist)

    # now flag any points in the *cylinders* that are within minDist of the remaining points of their end spheres 
    #for indexTestSet, testPoints in enumerate(cylinderXYZPoints): 
    #    for refPoints in np.concatenate((sphereXYZPoints[0], sphereXYZPoints[indexTestSet + 1]), 0):
     #       FlagPointsTooClose(cylinderPointFlagsList[indexTestSet], testPoints, refPoints, minDist)

    # apply the flags to the cylinders
    #cylinderXYZPoints = [ [point for point, flag in zip(pointList, cylinderPointFlagsList[index]) if flag==True] for index, pointList in enumerate(cylinderXYZPoints)]
    #cylinderDirectors = [ [director for director, flag in zip(directorList, cylinderPointFlagsList[index]) if flag==True] for index, directorList in enumerate(cylinderDirectors)]
    
    # add the cylinderPoints to the Sphere Points
    allPoints = cylinderXYZPoints #+ sphereXYZPoints
    allDirectors = cylinderDirectors #+ sphereDirectors
    
    # Return the flattened lists of points and directors, which by some minor miracle, should all correspond to each other
    return ([ point for minorList in allPoints for point in minorList], [ point for minorList in allDirectors for point in minorList]) 

    
# loops through a list of test points and flags for removal any that are within minDist of any of a set of reference points
def FlagPointsTooClose(FlagList, testPoints, refPoints, minDist):
    # Check the Flaglist and test points are equal lengths
    if len(FlagList)==len(testPoints):
        # loop through the test points 
        for index, testPoint in enumerate(testPoints):
            # if the test point is already flagged for removal then skip it
            if FlagList[index]==True:
                # check the test point against every refPoint.
                for refPoint in refPoints:
                    # if the distance between the test point and ref point is < minDist then remove it
                    if np.linalg.norm(testPoint - refPoint) < minDist:
                        FlagList[index]=False
    else:
        print("Flag list of unequal length to test list in FlagPointsTooClose")
        sys.exit()

    
# given a cloud of points determine which ones are inside a cylinder defined
# by a basepoint, a length and an axis.  (cylinder is circularly symmetric so rotation about axis is irrelevant).
# A list of flags is provided which is same length as the point cloud. Mark the corresponding flag as false if it is inside the cylinder
def FlagPointsInsideCylinder(FlagList, pointList, basePoint, axis, length, radius):
    
    if len(FlagList)==len(pointList):
        # loop through the points
        for index, point in enumerate(pointList):
            # find z component of point in cylindrical system
            zComponent = np.dot(axis, point - basePoint)
            
            # if z component is greater than zero or less than length it could be inside the cylinder
            # do the epsilon thing for rounding errors etc. assume anything within 1e-10 is actually same value 
            if zComponent > 1.0e-10 or zComponent <  length - 1.0e-10 :
                # only check radial component if z component doesn't rule it out
                if np.linalg.norm( point - basePoint - zComponent * axis ) < 0.9 * radius:
                    # if radial component is also smaller than radius then point is inside cylinder so flag it for removal 
                    FlagList[index]=False
    else:
        print("Flag list of unequal length to point list in FlagPointsInsideCylinder")
        sys.exit()

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
       
    numGen1Worms = 1
    xSize = 1000
    ySize = 1000
    zSize = 1000
    alpha1 = -5 
    alpha2 = 5
    beta1 = 160
    beta2 = 180
    wormRadius = 8
    segmentLength = 1
    unimerMinDist = 1.5
    packingFraction = 0.4
    minSegmentsPerWorm = 5
    coneAngleWorm = 45
    wormDirector = np.array([1.0, 0.0, 0.0])
    wormAlignmentAngularRange = 10
    numWormGenerations = 1
    minBranchAngle = 25
    wormBranchDensity = 0.1
    maxNumAttempts = 100
    
    WormNetworkXYZPoints, WormNetworkXYZDirectors= makeWormNetwork( numGen1Worms, 
                                                                    xSize, ySize, zSize, 
                                                                    alpha1, alpha2, beta1, beta2, 
                                                                    wormRadius, segmentLength, unimerMinDist, packingFraction,
                                                                    minSegmentsPerWorm, coneAngleWorm, wormDirector, wormAlignmentAngularRange, 
                                                                    numWormGenerations, minBranchAngle, wormBranchDensity, maxNumAttempts,
                                                                    selfAvoid=False, pointsToAvoid=[], visualiseEnvelope=(0,20), envelopeList=["None"])
    # generate the XYZVals packed in the outer cylinder
    fIO.saveXYZList(WormNetworkXYZPoints, ['C'] * len(WormNetworkXYZPoints), "wormNetwork.xyz")

#     from Projects.GeneralPolymerBrush import GeneralBlockPolymer as GBCP
# 
#     # generate unique polymer strands
#     for i in range(numPolymersPerCylinder):
#         print(i, "unimers out of: ", numPolymersPerCylinder)
#         if i==0:
#             UnimerStrands=[GeneralBlockPolymer(polymerBrushDict)]
#         else:
#             UnimerStrands.append(GeneralBlockPolymer(polymerBrushDict))
# 
#     # transform the strands to the cylinder points
#     curStrand = 0
#     for directorHat, pos, strand in zip(directorsHat, basePoints, UnimerStrands):
#         labRotation = rnd.uniform(0, 2 * np.pi)
#         xyzVals = coords.transformFromBlockFrameToLabFrame(directorHat, pos, labRotation, np.array([0.0, 0.0, 1.0]), strand[0][0], strand[0])
#         if curStrand==0:
#             xyzValsList = xyzVals
#             allNames = strand[1]
#         else:
#             xyzValsList = np.concatenate((xyzValsList, xyzVals), 0)
#             allNames = np.concatenate( (allNames, strand[1]), 0)
#         curStrand += 1
# 
#     fIO.saveXYZList(xyzValsList, allNames, "wormNetwork.xyz")

    
    print("example done")
