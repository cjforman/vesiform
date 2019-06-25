'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.constrainedPolymer import ConstrainedPolymerPackBBG as CPPBBG
from Library.VolumePackCuboid import VolumePackCuboidBBG as VPCBBG
import sys
import numpy as np
import random as rnd
import copy as cp
from Utilities import fileIO as fIO
from Utilities import coordSystems as coords
import networkx as nx

class BranchedPolymerPackBBG(BBG):
    '''
    Overloads the building block generator class, and makes use of the constrained polymer class,
    to generate a set of polymer chains between the nodes of arbitary graphs drawn 
    from a set of N points. Each distinct graph is called a cluster. 
    
    The user specifies a set containing the allowed numbers of nodes per cluster, and these 
    are distributed randomly. 
     
    Sets of 2,4 5 and 7 nodes gives rise to a 0, 1 2 or 3 branch points per cluster. 
    Each polymer that connects two nodes of a given graph is drawn using the constrained 
    polymer class. The rest of the network can be passed in as pointsToAvoid to avoid a self crossing.
    But this increases the cost of the calculation.
       
    '''
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        BBG.initialiseParameters(self) 
        
        self.particleName = self.getParam('particleName')

        self.CPPBBG = CPPBBG(self.paramFilename)
        self.VPCBBG = VPCBBG(self.paramFilename)
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for constrained polymer object")
            sys.exit()

    def generateBuildingBlock( self, 
                               numPoints,
                               x1,
                               x2,
                               y1,
                               y2,
                               z1,
                               z2,
                               networkMinDist,
                               neighborRadius,
                               monomerMinDist,
                               bondLength,
                               curlyFactor,
                               coneAngle,
                               minAngle,
                               threshold,
                               maxNeighbours,
                               angularRange=['None'],
                               pointsToAvoid=[],
                               visualiseEnvelope=(0,20),
                               envelopeList=["None"]):

        self.numPoints = numPoints
        self.bondLength = float(bondLength)
        self.networkMinDist = networkMinDist
        self.neighborRadius = neighborRadius
        self.monomerMinDist = monomerMinDist
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.z1 = z1
        self.z2 = z2
        self.coneAngle = coneAngle
        self.minAngle = minAngle
        self.threshold = threshold
        self.maxNeighbours = maxNeighbours
        self.curlyFactor = curlyFactor
        self.pointsToAvoid = pointsToAvoid
        self.visualiseEnvelope = visualiseEnvelope
        self.envelopeList = envelopeList
        self.angularRange = angularRange
        
        return BBG.generateBuildingBlock(self, numPoints, networkMinDist, envelopeList=envelopeList, visualiseEnvelope=visualiseEnvelope, pointsToAvoid=pointsToAvoid)

    def generateBuildingBlockDirector(self):
        return  np.array([0.0, 0.0, 1.0])
    
    def generateBuildingBlockRefPoint(self):
        return np.array( [ ( self.x1 + self.x2 ) / 2.0, 
                           ( self.y1 + self.y2 ) / 2.0, 
                           ( self.z1 + self.z2 ) / 2.0 ] )

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints

    def generateBuildingBlockXYZ(self):
        # pick n random points in the specified envelope within a cubic region of space.
        points = self.VPCBBG.generateBuildingBlock(self.numPoints, 
                                                   self.x1, 
                                                   self.x2, 
                                                   self.y1, 
                                                   self.y2, 
                                                   self.z1, 
                                                   self.z2, 
                                                   self.networkMinDist, 
                                                   envelopeList=self.envelopeList,
                                                   pointsToAvoid=self.envelopeList,
                                                   visualiseEnvelope=self.visualiseEnvelope)

        if self.dumpInterimFiles==1:
            fIO.saveXYZ(points.blockXYZVals, 'C', 'pointsInCube.xyz')

        # generate a random graph of connections between all the points,
        # such that no node has more than 3 connections and every node has at least one edge
        g = self.splitPointsIntoGraphs(points.blockXYZVals)
        
        if self.dumpInterimFiles==1:
            self.visualiseNetwork(g, points.blockXYZVals, self.monomerMinDist)

        return self.connectGraphWithConstrainedPolymer(g, points.blockXYZVals)

    def visualiseNetwork(self, g, points, radius):
        dZ = 0.1 * radius # ten particles per radius 
        xyzVals = []
        for edge in g.edges:
            node1 =  points[edge[0]]
            node2 =  points[edge[1]]
            director = node2 - node1
            length = np.linalg.norm(director)
            director = director/length
            numPoints = int(length/dZ)
            xyzVals += [  node1 + float(zIndex) * dZ * director for zIndex in range(0, numPoints)   ]
                
        fIO.saveXYZ(xyzVals, 'P', 'rawNetwork.xyz')        

    
    def splitPointsIntoGraphs(self, points):
        
        # nodes in graph are numbered 0 to numNodes - 1
        # which correspond to the index in points. 
        numNodes = len(points)
        
        # set up an empty graph
        g = nx.Graph()

        # add the nodes to the graph 
        for nodeNum in range(0, numNodes):
            g.add_node(nodeNum)

        params=[self.neighborRadius, self.coneAngle, self.minAngle, self.threshold, self.maxNeighbours]

        # Add the edges and return
        return self.add_edges(g, points, params, mode="natural")
    
    # switcher between edge adding algorithms    
    def add_edges(self, g, points, params, mode="natural"):
        
        if mode=="deterministic":
            outputGraph = self.add_edges_deterministic(g)
    
        if mode=="natural":
            outputGraph = self.add_edges_natural(g, points, params[0], params[1], params[2], params[3], params[4])

        # find singleton nodes
        nodesToRemove=[]
        for node in outputGraph.nodes():
            if sum(1 for _ in g.neighbors(node))==0:
                nodesToRemove += [node]

        # remove singleton nodes
        for node in nodesToRemove:
            outputGraph.remove_node(node)

        return outputGraph

    # randomly picks two points less than radius apart and uses that as the axis for a cone of interior angle coneangle
    # Then picks point inside that cone, such that no points subtend an angle less than minAngle.
    # keep going until threshold percentage of points have at least one edge
    def add_edges_natural(self, g, points, neighborRadius, coneAngle, minAngle, threshold, maxEdgesPerNode):
    
        numPoints = len(points)

        numFreeNodes = self.countFreeNodes(g)
        minNumFreeNodes = int((1.0 - threshold) * float(numPoints))    
                
        while ( numFreeNodes > minNumFreeNodes):  

            # pick a node at random
            node1 = rnd.randint(0, numPoints-1)
            
            # only pick a node if it has less than maxEdgesPerNode
            if sum(1 for _ in g.neighbors(node1)) < maxEdgesPerNode:
                
                # pick a second node at random
                node2 = rnd.randint(0, numPoints-1)
                
                # compute director between vector if node1 and node 2 are not the same node
                if not node1==node2:
                    director = points[node2] - points[node1]
                    length = np.linalg.norm(director)
                    directorHat = director/length 
                
                # only pick second node if it not node 1, has less than maxEdgesPerNode and is within radius of the first node 
                if (not node1==node2) and (sum(1 for _ in g.neighbors(node2)) < maxEdgesPerNode) and length < neighborRadius:
                    
                    # add edge between the two root nodes
                    g.add_edge(node1, node2)

                    # begin a list of points in the current chain that we are making 
                    xyzList = [ points[node1], points[node2] ]

                    # generate a list of points that are in the allowed cone                    
                    coneList = self.pickPointsInCone(points[node1], coneAngle, directorHat, points, xyzList)

                    # loop through the cone list checking each point in order to see
                    # a) it's less than radius from last point on list
                    # b) forms an angle greater than minAngle with previous 2 points on list
                    # c) has less edges than maxEdgesPerNode
                    lastNode = node2
                    for node in coneList:
                        if (np.linalg.norm(points[node] - xyzList[-1]) < neighborRadius ): 
                            if coords.bondAngle(xyzList[-2], xyzList[-1], points[node]) > minAngle:
                                if (sum(1 for _ in g.neighbors(node2)) < maxEdgesPerNode):
                                    g.add_edge(lastNode, node)
                                    xyzList += points[node]
                                    lastNode = node
            numFreeNodes = self.countFreeNodes(g)

        return g
                    
    # returns a list of indices of points in pointsList that are within a cone
    # defined by apex, director and cone angle.
    # Sorts the list of indices by the distance from the apex point along the cone axis
    def pickPointsInCone(self, apex, coneAngle, director, pointsList, excludedPoints):
        indexList = []
        tanAngle = np.tan(coneAngle/2)
        for index, point in enumerate(pointsList):
            
            if not True in [ np.linalg.norm(point - testPoint)<1e-6 for testPoint in excludedPoints ]: # ignore the two points that were used to define the cone 
                z = np.dot(point - apex, director)
                y = np.linalg.norm(point - apex - z*director)
                if y/z <= tanAngle and z > float(0.0) :
                    indexList += [(index, z)]

        indexList = sorted(indexList, key=lambda a:a[1])
        return [ item[0] for item in indexList]
                    
            
    def countFreeNodes(self, g):
        return sum( 1 for node in g.nodes if sum(1 for _ in g.neighbors(node))==0 )
 
    
    def add_edges_deterministic(self, g):
        # loop through the nodes once and add one two or 3 edges to the node.
        for node in g.nodes():
            
            # randomly choose how many edges this node should have
            numEdgesToAdd = rnd.choice([1,2,3])
            
            # count number of current edges on the node
            numCurrentEdges = sum(1 for _ in g.neighbors(node))
            
            # reduce number of edges 
            numEdgesToAdd  = numEdgesToAdd - numCurrentEdges
            
            # only add edges to this node if we still need to 
            if numEdgesToAdd>0:
            
                # finds, or creates the right number of nodes that can 
                # still be connected to (excluding the current node).
                freeNodes = self.findFreeNodes(g, [node], numEdgesToAdd)
                
                # connect the chosen nodes to the current node 
                for node2 in freeNodes:
                    # add an edge between the node and the selected free nodes
                    g.add_edge(node, node2)
            
        return g

    def findFreeNodes(self, g, protectedNodes, targetNumFreeNodes):

        freeNodes = []

        # make sure the requested number of nodes is smaller than the number of nodes
        # we are allowed to look at (total # of nodes - # of protected nodes)
        if targetNumFreeNodes <= ( len(g) - len(protectedNodes) ): 
            
            # first look for nodes with no neighbours, then 1, then 2.
            for targetNumEdges in [0, 1, 2]:
                # loop through the nodes (backwards!) and look for nodes with targetNumEdges
                # backwards encourages a more broken up network with several components 
                for node in range(len(g)-1, 0, -1):
                    # make sure the current node is not on the protected list 
                    if not node in protectedNodes:
                        # if the current node has the targetNumEdges then add it to the free 
                        # node list.
                        if (sum(1 for _ in g.neighbors(node))==targetNumEdges):
                            freeNodes+=[node]
        
                    # break the inner loop if we already have enough nodes                    
                    if len(freeNodes)==targetNumFreeNodes:
                        break
                # break the outer loop if we already have enough nodes                    
                if len(freeNodes)==targetNumFreeNodes:
                    break
    
            # if we got through the 0, 1 and 2 nodes and there aren't enough free nodes
            # we need to disconnect some edges.   
            # If we are in this state then we know that all the nodes 
            # have 3 edges and there are enough non-protected nodes to generate the free_list
            if (len(freeNodes)<targetNumFreeNodes)==True:
                # loop through unprotected nodes and remove an edge
                # add the current node to the free nodes.
                # The second node is also free, but don't want to do an insertion.
                # might make for a dull networks.So this will leave 3 nodes open - will help
                # with later nodes anyway. 
                for node in g.nodes():
                    if not node in protectedNodes:
                        for node2 in g.neighbours(node):
                            g.remove_edge(node, node2)
                            break # only take the first neighbor
                    if len(freeNodes)<targetNumFreeNodes:
                        freeNodes+=node
                    if len(freeNodes)>=targetNumFreeNodes:
                        break;
        else:
            print("Impossible to find enough free nodes to connect to.")
            exit()

        if len(freeNodes)!=targetNumFreeNodes:
            print("Warning: numfreeNodes not equal to target number but run out of things to try")

        return freeNodes


    # checks to see how many nodes have no neighbours.
    def numZeroNodes(self, g):
        numZeroNodes = 0
        for node in g.nodes():
            if sum(1 for _ in g.neighbors(node))==0:
                numZeroNodes+=1
        return numZeroNodes
        
    def connectGraphWithConstrainedPolymer(self, g, points):
        '''
        # we now traverse the graph generating a polymer for each edge in the network.
        # the nodes are the indices of the A and B points for each network, and the distance between
        # the points determines the minmum length of polymer, with extra length allowed by the curly factor.
        # We supply the growing list of points as points to avoid for subsequent chains.
        # Might make some impossible situations, so choose parameters carefully.
        '''
        
        # create an array to store output
        outputXYZ = []  
        
        # specify number of crank moves to try to fold polymer in to zone
        numCrankMoves = 0
        
        # loop through the edges in g and generate a self avoiding random polymer
        # that connects the specified points in the space
        for edge in g.edges():
            pointA = points[edge[0]]
            pointB = points[edge[1]]
            
            # generate local points to avoid 
            
            # copy externally supplied points
            l_pointsToAvoid = cp.copy(self.pointsToAvoid) 
            
            # add existing output from other polymers providing they are more than minDist away from pointA and pointB
            for point in outputXYZ:
                if np.linalg.norm(point - pointA) > self.minDist and np.linalg.norm(point - pointB) > self.monomerMinDist:   
                    l_pointsToAvoid += [point] 
            
            # add the other node points in the network providing they are far enough away from pointA and pointB
            for point in points:
                if np.linalg.norm(point - pointA) > self.minDist and np.linalg.norm(point - pointB) > self.monomerMinDist:   
                    l_pointsToAvoid += [point] 

            # compute the number of vertices we will need between the two points.        
            numPoints = int(self.curlyFactor * np.linalg.norm(pointA - pointB)/self.bondLength)
        
            polymerBB = self.CPPBBG.generateBuildingBlock(  numPoints, 
                                                            pointA, 
                                                            pointB, 
                                                            self.monomerMinDist, 
                                                            self.bondLength, 
                                                            numCrankMoves, 
                                                            l_pointsToAvoid, 
                                                            visualiseEnvelope=self.visualiseEnvelope, 
                                                            envelopeList=self.envelopeList,
                                                            angularRange=self.angularRange)
            
            # add the newly minted points to the list of output points. 
            outputXYZ += [ point for point in polymerBB.blockXYZVals]

        return outputXYZ

if __name__ == '__main__':
   
    # create the generator
    BPPBBG = BranchedPolymerPackBBG('branchedPolymerNetwork.txt')
    numPoints = 100
    x1 = -100
    x2 = 100
    y1 = -100
    y2 = 100
    z1 = -100
    z2 = 100
    networkMinDist = 15.0
    neighborRadius = networkMinDist * 5 
    monomerMinDist = 1.0
    bondLength = 2.0
    curlyFactor = 1.4
    pointsToAvoid=[]
    visualiseEnvelope=(0,20)
    envelopeList=["None"]
    angularRange=[0.0, 10, 165, 175]
    # envelopeList=["outersphere 12.0", "innersphere 10.0"],
    # pointsToAvoid=pointsToAvoid
    coneAngle = 45 * np.pi/180
    minAngle = 165 * np.pi/180
    threshold = 0.7
    maxNeighbours = 3
    
        
    # generate a set of points and build a constrained polymer pack network between them
    BranchedPolymerPackBB = BPPBBG.generateBuildingBlock(  numPoints, 
                                                           x1, x2,
                                                           y1, y2,
                                                           z1, z2, 
                                                           networkMinDist,
                                                           neighborRadius,
                                                           monomerMinDist, 
                                                           bondLength,
                                                           curlyFactor,
                                                           coneAngle,
                                                           minAngle,
                                                           threshold,
                                                           maxNeighbours)

    fIO.saveXYZ(BranchedPolymerPackBB.blockXYZVals, 'C', "branchedPolymerNetwork.xyz")
    print("branched Polymer Done")