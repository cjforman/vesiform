import numpy as np
import random as rnd
import sys
import Utilities.cartesian as cart
from Library.VolumePackCuboid import VolumePackCuboidBBG as VPCBBG
from Library.SurfacePackCylinder import SurfacePackCylinderBBG as SPCBBG
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords

def makeWormNetwork(numNodesInit, nodeDist, xSize, ySize, zSize, WormRadius, MaxTubeLength, MaxNumAttempts, ValidNodesPerCluster, minDist, packingFraction):

    VolumeBBG = VPCBBG('VolumePackCuboid.txt')
  
    # generate the XYZVals of the network nodes packed in the volume
    NetworkNodesXYZ = VolumeBBG.generateBuildingBlock(numNodesInit, -xSize/2, xSize/2, -ySize/2, ySize/2, -zSize/2, zSize/2, nodeDist)

    fIO.saveXYZ(NetworkNodesXYZ.blockXYZVals, 'O', 'cubepackedPoints.xyz')
        
    # create a list of groups of indices which together form a good network
    ListOfNetworkGroups = groupNodesInNetwork(NetworkNodesXYZ.blockXYZVals, WormRadius, MaxTubeLength, ValidNodesPerCluster, MaxNumAttempts)

    # dump output as file
    outputWormNetworkAsAtoms(ListOfNetworkGroups, NetworkNodesXYZ.blockXYZVals, TubeRadius)


    # for each cluster in the network generate a series of points and directors 
    networkPoints = populateClusters(ListOfNetworkGroups, NetworkNodesXYZ.blockXYZVals, WormRadius, minDist, packingFraction)

    # repack each list
    xyzVals = [ points[0] for points in networkPoints]
    directors = [ points[1] for points in networkPoints]
    names = [ points[2] for points in networkPoints]

    # return flattened lists
    return ([ point for minorList in xyzVals for point in minorList ], [ director for minorList in directors for director in minorList ], [ name for minorList in names for name in minorList ] )


# debugging routine for visualising the network 
def outputWormNetworkAsAtoms(ListOfNetworkGroups, XYZList, radius):
    dZ = 0.1 * radius # ten particles per radius 
    xyzVals = np.array([])
    for cluster in ListOfNetworkGroups:
        rootNodePos = XYZList[cluster[0]]
        for endNode in cluster[1:]:
            endNodePos = XYZList[endNode]
            director = endNodePos - rootNodePos
            length = np.linalg.norm(director)
            director = director/length
            numPoints = int(length/dZ)
            if xyzVals.size==0:
                xyzVals = [  rootNodePos + float(zIndex) * dZ * director for zIndex in range(0, numPoints)   ]
                xyzVals = np.concatenate( (xyzVals, [endNodePos]), 0 )
            else:
                xyzVals = np.concatenate( (xyzVals, [rootNodePos + float(zIndex) * dZ * director for zIndex in range(0, numPoints)]), 0)

    fIO.saveXYZ(xyzVals, 'P', 'rawNetwork.xyz')
        

def populateClusters(clusterIndexList, spherePoints, radius, minDist, packingFraction):
    # take the list of clusters and generate a series of unimer base points and directors 
    # across the surface of all the clusters in the network
    clusterXYZ = []
    for cluster in clusterIndexList:
        clusterPoints = [ spherePoints[index] for index in cluster]
        if not clusterXYZ:
            clusterXYZ = [ makeCluster(clusterPoints, radius, minDist, packingFraction) ]
        else:
            clusterXYZ.append(makeCluster(clusterPoints, radius, minDist, packingFraction))
    
    return clusterXYZ


# Generates an array of points, and surface normal directors for each point, which are on the surface 
# of a complex constructions of spheres and cylinders which form a connected worm like
# network. 


# The wormlike network is defined by a list of nodes which are the centres of spheres of
# a given radius. The node in the list is the root node and all the other nodes
# in the list are connected by a cylinder of the same given radius.
#
# A list of random points is generated across the segmented surface - defined by the union of the spheres around each node
# and the connecting cylinders - with a uniform surface density (defined by minDist).  
# All points on individual surfaces which are inside the composite surface are removed.
# 
# The Output is a flat list of points and their corresponding directors which define the wormlike network 
#
def makeCluster(nodes, radius, minDist, packingFraction):
    CylinderBBG = SPCBBG('SurfacePackCylinder.txt')
    SphereBBG = SPSBBG('SurfacePackSphere.txt')
    
    # useful debug tool
    removeClosePointsFlag = True
    removeInnerPointsFlag = True
    includeSpheres = True
    includeCylinders = True
    
    # determine the point density so it is the same on spheres and cylinders
    # number of points on Sphere = Area of sphere divided by area of circles of radius minDist.
    # np.pi's cancel out. Apply the packing fraction is the ratio relative to perfect packing
    nS = int(packingFraction * 4.0 * (radius**2)/ minDist**2)
    
    # compute number of points per unit length for each cylinder
    # circumference * 1 unit length / base area of each point on surface
    # the np.pi's cancel out
    nCUnit_float = packingFraction * (2.0 * radius) / minDist**2

    # set the root point of the cluster
    rootPoint = nodes[0]
    
    if includeSpheres:
        # populate each sphere with specified number of points minDist apart
        spheres = [ SphereBBG.generateBuildingBlock(nS, radius, -90.0, 90.0, -180.0, 180.0, minDist) for _ in range(0,len(nodes)) ] 
        
        # compute an array of the sphere directors. Each sphere centred at origin so each pos is its own radial vector.
        # just normalise it and return in seperate array.
        sphereDirectors = [ [ pos/np.linalg.norm(pos) for pos in sphere.blockXYZVals] for sphere in spheres ]

        # translate the sphere points so the centre of each sphere is at the given node.
        # Return an array the same shape as the director array for each point
        sphereXYZPoints = [ [nodes[index] + pos for pos in sphere.blockXYZVals] for index, sphere in enumerate(spheres)]

        # generate sphere names
        sphereNames = [  ['S'] * len(sphere.blockXYZVals) for sphere in spheres] 

    if includeCylinders:    
        # define the length of each cylinder in the cluster    
        cylinderLengths = [ np.linalg.norm(rootPoint - node) for node in nodes[1:] ]
    
        # define the axial orientation of each cylinder as the unit vector from the rootPoint to the outer node
        cylinderAxes = [ (node - rootPoint)/length for node, length in zip(nodes[1:], cylinderLengths) ]
    
        # populate each cylinder with a specified number of points no less than minDist apart    
        cylinders = [ CylinderBBG.generateBuildingBlock(int(length * nCUnit_float), radius, radius, 0, length, -180.0, 180.0, minDist) for length in cylinderLengths ]
        
        # transform the cylinder points so the cylinder axis is aligned with the director away from the central node to the outer node 
        [ cylinderBB.transformBBToLabFrame(director, rootPoint, 0.0) for cylinderBB, director in zip(cylinders, cylinderAxes)]
        
        # extract the cylinder points from the building block object
        cylinderXYZPoints = [ cylinderBB.blockXYZVals for  cylinderBB in cylinders]
        
        # compute the director of each point (radial vector pointing away from the cylinder axis (rodrigues formula) (see page 183 in book 5).
        cylinderDirectors = [  [pos - rootPoint - np.dot(director, pos - rootPoint) * director for pos in cylinder] for cylinder, director in zip(cylinderXYZPoints, cylinderAxes) ]
        
        # normalise the directors
        cylinderDirectors = [  [pos/np.linalg.norm(pos) for pos in cylinder] for cylinder in cylinderDirectors ]

        # generate cylinder names
        cylinderNames = [  ['C'] * len(cylinder.blockXYZVals) for cylinder in cylinders]
    
    # For each list of points that we have constructed (spheres and cylinders) 
    # go through and remove any of those points that are *inside* the cylindrical shells.
    # None of the spheres would have other sphere points or cylinder point inside them so no need to check them.
    
    # first deal with the points belonging to the cylindrical surfaces

    if includeCylinders:    
        # assume that all the points will survive - set up a flag array for each point that is same dims as director and points array
        cylinderPointFlagsList = [ len(cylinderPoints) * [True] for cylinderPoints in cylinderXYZPoints  ]
    
    if removeInnerPointsFlag and includeCylinders:
        # loop through the sets of points for each cylinder
        for indexInner, cylinderPoints in enumerate(cylinderXYZPoints):
            
            # loop through the geometric regions
            for director, length, indexOuter in zip(cylinderAxes, cylinderLengths, range(0, len(cylinderLengths))):
                # don't check a cylinder against itself
                if indexInner!=indexOuter:
                    # test the current list of points against the cylinder
                    # Set flags in the flag array for points that need to be removed   
                    print("Removing points from cylinder:", indexInner + 1, " of ", len(cylinderXYZPoints), " that are inside cylinder: ", indexOuter + 1, "of ", len(cylinderLengths))
                    FlagPointsInsideCylinder(cylinderPointFlagsList[indexInner], cylinderPoints, rootPoint, director, length, radius)

    if includeCylinders:
        # apply the flags
        cylinderXYZPoints = [ [point for point, flag in zip(pointList, cylinderPointFlagsList[index]) if flag==True] for index, pointList in enumerate(cylinderXYZPoints)]
        cylinderDirectors = [ [director for director, flag in zip(directorList, cylinderPointFlagsList[index]) if flag==True] for index, directorList in enumerate(cylinderDirectors)]
        cylinderNames = [ [name for name, flag in zip(nameList, cylinderPointFlagsList[index]) if flag==True] for index, nameList in enumerate(cylinderNames)]
        # reconstruct the flag list as all true
        cylinderPointFlagsList = [ len(cylinderPoints) * [True] for cylinderPoints in cylinderXYZPoints  ]
                
    # do the same for the sphere points
    if includeSpheres:
        # assume that all the points will survive - set up a flag array for each point that is same dims as director and points array
        spherePointFlagsList = [ len(spherePoints) * [True] for spherePoints in sphereXYZPoints  ]

    if removeInnerPointsFlag and includeSpheres and includeCylinders:
        # loop through the sets of points for each sphere
        for indexInner, spherePoints in enumerate(sphereXYZPoints):
                    
            # loop through the geometric regions
            for director, length, indexOuter in zip(cylinderAxes, cylinderLengths, range(0,len(cylinderLengths))):
                print("Removing points from sphere:", indexInner + 1, " of ", len(sphereXYZPoints), " that are inside cylinder: ", indexOuter + 1, "of ", len(cylinderLengths))
                # test the current list of points against the cylinder 
                # set flags in the flag array for points that need to be removed   
                FlagPointsInsideCylinder(spherePointFlagsList[indexInner], spherePoints, rootPoint, director, length, radius)
 
    if includeSpheres:
        # apply the flags
        sphereXYZPoints = [ [point for point, flag in zip(pointList, spherePointFlagsList[index]) if flag==True] for index, pointList in enumerate(sphereXYZPoints)]
        sphereDirectors = [ [director for director, flag in zip(directorList, spherePointFlagsList[index]) if flag==True] for index, directorList in enumerate(sphereDirectors)]
        sphereNames = [ [name for name, flag in zip(nameList, spherePointFlagsList[index]) if flag==True] for index, nameList in enumerate(sphereNames)]
        # reconstruct the flag list as all true
        spherePointFlagsList = [ len(spherePoints) * [True] for spherePoints in sphereXYZPoints  ]

    if removeClosePointsFlag and includeCylinders:
        # now flag any points in the cylinders that are within minDist of the remaining points of the other cylinders. 
        # Don't check a cylinder against a cylinder that has already been checked.
        for indexTestSet, testPoints in enumerate(cylinderXYZPoints): 
            for indexRefSet, refPoints in enumerate(cylinderXYZPoints[indexTestSet + 1:]):
                print("Removing points in cylinder:", indexTestSet + 1, " of ", len(cylinderXYZPoints), " due to proximity with cylinder: ", indexRefSet + indexTestSet + 1)
                FlagPointsTooClose(cylinderPointFlagsList[indexTestSet], testPoints, refPoints, minDist)
    
    if removeClosePointsFlag and includeCylinders and includeSpheres:
        # now flag any points in the *cylinders* that are within minDist of the remaining points of their end spheres 
        for indexTestSet, testPoints in enumerate(cylinderXYZPoints): 
            for indexRefSet, refPoints in enumerate([sphereXYZPoints[0], sphereXYZPoints[indexTestSet + 1]]):
                print("Removing points in cylinder:", indexTestSet + 1, " of ", len(cylinderXYZPoints), " due to proximity with sphere: ", indexRefSet + 1, "of 2.")
                FlagPointsTooClose(cylinderPointFlagsList[indexTestSet], testPoints, refPoints, minDist)

    if includeCylinders:
        # apply the flags to the cylinders
        cylinderXYZPoints = [ [point for point, flag in zip(pointList, cylinderPointFlagsList[index]) if flag==True] for index, pointList in enumerate(cylinderXYZPoints)]
        cylinderDirectors = [ [director for director, flag in zip(directorList, cylinderPointFlagsList[index]) if flag==True] for index, directorList in enumerate(cylinderDirectors)]
        cylinderNames = [ [name for name, flag in zip(nameList, cylinderPointFlagsList[index]) if flag==True] for index, nameList in enumerate(cylinderNames)]    

    allPoints = []
    allDirectors = []
    allNames = []
    
    # compile output arrays
    if includeCylinders:
        # add the cylinderPoints to output arrays
        allPoints += cylinderXYZPoints
        allDirectors += cylinderDirectors
        allNames += cylinderNames
    
    if includeSpheres:
        allPoints += sphereXYZPoints
        allDirectors += sphereDirectors
        allNames += sphereNames
        
    # Return the flattened lists of points and directors, which by some minor miracle, should all correspond to each other
    return ([ point for minorList in allPoints for point in minorList], [ point for minorList in allDirectors for point in minorList], [ name for minorList in allNames for name in minorList]) 



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
            
            # if z component is greater than zero and less than length it could be inside the cylinder
            # do the epsilon thing for rounding errors etc. assume anything within 1e-6 is actually same value 
            if zComponent > 1.0e-6 and zComponent <  length - 1.0e-6 :
                # only check radial component if z component doesn't rule it out
                if np.linalg.norm( point - basePoint - zComponent * axis ) - radius < -1.0e-6:
                    # if radial component is also smaller than radius then point is inside cylinder so flag it for removal 
                    FlagList[index]=False
    else:
        print("Flag list of unequal length to point list in FlagPointsInsideCylinder")
        sys.exit()

#### checked #####

def groupNodesInNetwork(xyzList, tubeRadius, maxTubeLength, validClusterLengths, numAttempts):

    # generate outlist
    outList = []
    rootNodesRemaining = [ i for i in range(0, len(xyzList))]
    branchNodesRemaining = [ i for i in range(0, len(xyzList))]
    numLoopsWithoutANewCluster = 0

    # keep looping until we stop finding clusters. 
    # This will be a bit hit and miss what a good threshold is. but we don't care that
    # much as long as we get a decent number of valid clusters.
     
    while numLoopsWithoutANewCluster < numAttempts:
        
        # pick a point at random in the root nodes list
        rootNodeIndex = rnd.choice(rootNodesRemaining)
        
        # choose cluster length randomly from list of valid cluster length 
        numNodes = rnd.choice(validClusterLengths)
        
        # create a test node. Makes sure all the points are within maxTubeLength of index node
        testCluster = pickTestCluster(rootNodeIndex, branchNodesRemaining, xyzList, maxTubeLength, numNodes)

        # check to see if we found enough branch nodes to make a valid cluster
        if len(testCluster) in validClusterLengths: 
            # we produced a long enough cluster, does it intersect 
            # existing clusters. 
            if checkTestClusterAgainstList(testCluster, outList, xyzList, tubeRadius):
                # cluster does not intersect existing clusters so we found a good one.
                # add list of indices forming the cluster to the output list
                outList.append(testCluster)

                # reset the counter
                numLoopsWithoutANewCluster = 0
                
                # remove root node from the root node list
                rootNodesRemaining.remove(rootNodeIndex)
                
                # remove each node index (including root node) 
                # from the possible branch nodes list
                for index in testCluster:
                    branchNodesRemaining.remove(index)
                
                if (len(rootNodesRemaining)==0) or (len(branchNodesRemaining)<min(validClusterLengths)):
                    # force an exit to the while loop as there are no points left to assign to clusters
                    numLoopsWithoutANewCluster=MaxNumAttempts
        else:
            # the randomly chosen node will never make a good root 
            # (not enough branch nodes within appropriate distance)
            # so remove it from the remaininRootNodeList.
            # but it could still be a useful branch
            rootNodesRemaining.remove(rootNodeIndex)
            if (len(rootNodesRemaining)==0):
                # force an exit to the while loop as there are no points left to assign 
                numLoopsWithoutANewCluster=MaxNumAttempts

        numLoopsWithoutANewCluster += 1
        
    return outList

def checkTestClusterAgainstList(testCluster, outList, xyzList, tubeRadius):
    # assume success
    retVal = True
     
    # loop through all the clusters in the list
    for cluster in outList:
        # check to see if the current cluster intersects anywhere with the new cluster
        retVal = doClustersIntersect(testCluster, cluster, xyzList, tubeRadius)
        
        # if it does then break 
        if retVal==False:
            break
    
    return retVal
    
def doClustersIntersect(cluster1, cluster2, xyzList, radius):

    # assume success
    goodCluster = True
    
    # first point in each cluster is common to all tubes in the cluster
    p1 = xyzList[cluster1[0]]
    p2 = xyzList[cluster2[0]]

    # loop through all combinations of end points until we find one that intersects    
    for endPointIndexC1 in cluster1[1:]:
        if goodCluster==True:
            for endPointIndexC2 in cluster2[1:]:
                goodCluster = doSpheroCylindersIntersect(p1, xyzList[endPointIndexC1], p2, xyzList[endPointIndexC2], radius)
                if goodCluster ==False:
                    break

    return goodCluster

def doSpheroCylindersIntersect(p1, q1, p2, q2, radius):
    # tests to see if two spherocylinders intersect.
    # the points are the centres of the two spheres at either end of the sphero cylinder.
    # p1 and q1 are the start and end of one spherocylinder, p2 and q2 the start and 
    # end of the second.    
    goodPos = True
    if cart.closestApproachTwoLineSegmentsSquared(p1, q1, p2, q2) < (2.0 * radius)**2:
        goodPos = False
        
    return goodPos


def pickTestCluster(rootNodeIndex, validIndexList, xyzList, maxTubeLength, numNodes):    

    outListIndex = [rootNodeIndex] 
    rootNode = xyzList[rootNodeIndex]
    
    # loop through the list and make note of the indices of at most numNodes 
    # nodes that are within maxTubeLength of the rootnode 
    for nodeIndex in validIndexList:
        if nodeIndex!=rootNodeIndex:
            if np.linalg.norm(rootNode - xyzList[nodeIndex]) < maxTubeLength:
                outListIndex.append(nodeIndex)  
        if len(outListIndex)==numNodes:
            break

    return outListIndex        


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
    backboneDict1['minDist'] = 0.5
    backboneDict1['bondLength'] = 0.5
    backboneDict1['Z1'] = 2
    backboneDict1['R1'] = 2
    backboneDict1['Z2'] = 1.5 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
    backboneDict1['R2'] = 10

    backboneDict2['filename'] = 'RandomPolymer.txt'  
    backboneDict2['mode'] = 'Polymer'
    backboneDict2['name'] = 'BlockB'
    backboneDict2['numMonomers'] = 6
    backboneDict2['monomerNames'] = backboneDict2['numMonomers'] * ['C']
    backboneDict2['alpha1'] = 40
    backboneDict2['alpha2'] = 50
    backboneDict2['beta1'] = 165
    backboneDict2['beta2'] = 185
    backboneDict2['minDist'] = 0.5
    backboneDict2['bondLength'] = 0.5
    backboneDict2['Z1'] = 2
    backboneDict2['R1'] = 2
    backboneDict2['Z2'] = 1.5 * backboneDict2['numMonomers'] * backboneDict2['bondLength'] + backboneDict2['Z1']
    backboneDict2['R2'] = 10

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
    brushDict1['bondLength'] = 1.0
    brushDict1['Z1'] = 2.0
    brushDict1['R1'] = 2.0
    brushDict1['Z2'] = brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
    brushDict1['R2'] = 10.0

    brushDict2['filename'] = 'RandomPolymer.txt'  
    brushDict2['mode'] = 'Polymer'
    brushDict2['name'] ='brush2'
    brushDict2['numMonomers'] = 4 
    brushDict2['monomerNames'] = ['Zn', 'O', 'Pb', 'N']#, 'B', 'H', 'S']
    brushDict2['phaseRange'] = 360.0 
    brushDict2['alpha1'] = 40
    brushDict2['alpha2'] = 50
    brushDict2['beta1'] = 155
    brushDict2['beta2'] = 175
    brushDict2['minDist'] = 0.5
    brushDict2['bondLength'] = 0.5
    brushDict2['Z1'] = 2
    brushDict2['R1'] = 50
    brushDict2['Z2'] = 1.5 * brushDict2['numMonomers'] * brushDict2['bondLength'] + brushDict2['Z1']
    brushDict2['R2'] = 50
              
    polymerBrushDict={}
    polymerBrushDict['backbones'] = [backboneDict1, backboneDict2]#, backboneDict1, backboneDict2]
    polymerBrushDict['brushes'] = [brushDict1, brushDict2]#, brushDict1, brushDict2]
    polymerBrushDict['connectors'] = [connectorDict12]#, connectorDict12, connectorDict12]
       
    numNodesInit = 100
    TubeRadius = 8
    nodeDist = 2 *(TubeRadius + backboneDict1['Z2'] - backboneDict1['Z1'] + backboneDict2['Z2'] - backboneDict2['Z1'])
    xSize = 1000  
    ySize = 1000
    zSize = 200
    MaxTubeLength = 500
    MaxNumAttempts = 1000
    ValidNodesPerCluster = [4]
    unimerBaseDist = 2
    packingFraction = 0.01    
    # generate the XYZVals packed in the outer cylinder
    WormNetwork = makeWormNetwork(numNodesInit, nodeDist, xSize, ySize, zSize, TubeRadius, MaxTubeLength, MaxNumAttempts, ValidNodesPerCluster, unimerBaseDist, packingFraction)
    fIO.saveXYZList(WormNetwork[0], WormNetwork[2], "wormBasePointNetwork.xyz")
    fIO.saveXYZ([ point + director for point, director in zip(WormNetwork[0], WormNetwork[1])], 'O'  ,"directorNetwork.xyz")

    from Projects.GeneralPolymerBrush import GeneralBlockPolymer as GBCP
 
    numUnimers = len(WormNetwork[0])
    UnimerStrands = []
    # generate unique polymer strands
    for basePos, director, index in zip(WormNetwork[0], WormNetwork[1], range(0, numUnimers)):
        print(index, " unimers out of: ", numUnimers)
        UnimerStrands+=[GBCP(polymerBrushDict)]
        
    # transform the strands to the worm network positions
    xyzValsList = []
    allNames = []
    for directorHat, pos, strand in zip(WormNetwork[1], WormNetwork[0], UnimerStrands):
        labRotation = rnd.uniform(0, 2 * np.pi)
        # add numpy points to python list 
        xyzValsList += [ pos for pos in coords.transformFromBlockFrameToLabFrame(directorHat, pos, labRotation, np.array([0.0, 0.0, 1.0]), strand[0][0], strand[0])]
        # add numpy array of strings to python list
        allNames += strand[1].tolist()

    fIO.saveXYZList(xyzValsList, allNames, "wormNetwork.xyz")
    
    print("example done")
