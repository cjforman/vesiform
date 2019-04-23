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
                               monmerMinDist,
                               bondLength,
                               curlyFactor,
                               angularRange=['None'],
                               pointsToAvoid=[],
                               visualiseEnvelope=(0,20),
                               envelopeList=["None"]):

        self.numPoints = numPoints
        self.bondLength = float(bondLength)
        self.networkMinDist = networkMinDist
        self.monomerMinDist = monomerMinDist
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.z1 = z1
        self.z2 = z2
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

        return self.connectGraphWithConstrainedPolymer(g, points.blockXYZVals)

    # uses the sizeSet information to breaks the points list up
    # into Ni graphs of size Mi where 0<i<len(sizeSet).
    # sizeSet is a list of tuples [(N1, M1), (N2, M2)...]
    # such that Ni is number of graphs of that size and M2 is number of points in that
    # graph.  
    def splitPointsIntoGraphs(self, points):
        
        # nodes in graph are numbered 0 to numNodes - 1
        # which correspond to the index in points. 
        numNodes = len(points)
        
        # set up an empty graph
        g = nx.Graph()

        # add the nodes to the graph 
        for nodeNum in range(0, numNodes):
            g.add_node(nodeNum)

        # Add the edges.
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
    numPoints = 30
    x1 = -100
    x2 = 100
    y1 = -100
    y2 = 100
    z1 = -10
    z2 = 10
    networkMinDist = 30.0
    monomerMinDist = 1.0
    bondLength = 2.0
    curlyFactor = 1.4
    pointsToAvoid=[]
    visualiseEnvelope=(0,20)
    envelopeList=["None"]
    angularRange=[0.0, 10, 165, 175]
    # envelopeList=["outersphere 12.0", "innersphere 10.0"],
    # pointsToAvoid=pointsToAvoid
        
    # generate a set of points and build a constrained polymer pack network between them
    BranchedPolymerPackBB = BPPBBG.generateBuildingBlock(  numPoints, 
                                                           x1, x2,
                                                           y1, y2,
                                                           z1, z2, 
                                                           networkMinDist,
                                                           monomerMinDist, 
                                                           bondLength,
                                                           curlyFactor)

    fIO.saveXYZ(BranchedPolymerPackBB.blockXYZVals, 'C', "branchedPolymerNetwork.xyz")
    print("branched Polymer Done")