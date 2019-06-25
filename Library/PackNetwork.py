'''
Created on 30 Apr 2019

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import random as rnd
import copy as cp
import numpy as np
import networkx as nx
from Utilities import fileIO as fIO
from Utilities import cartesian as cart


class PackNetworkBBG(NBB):
    '''
    Over loads the NPackBB object to give a packing on a spatial locus which is 
    defined as a function G of a network and set of spatial points where no point in the 
    locus has a euclidean distance less than minDist from the others. A set of parameters 
    must also be supplied which are used by the function of the network to map the network
    to a spatial locus.
    
    i.e. the set of points defined by G(network, points, params).  
    
    The input points are geometric spatial points, and the network defines 
    which edges between the points are connected.  The two spatial points connected by an 
    edge form a straight line in 3-space. The network of lines connecting the points is considered
    to be the object that is being packed.  It is up to the user
    to overload the function "checkPointStatus" which tests whether a test point is 
    on the surface, inside or outside G(network, points, params). 
    
    The user can specify (via the packType parameter) which set of points to pack: either 'surface', 'inside' or 'outside'.    
    
    Only one kind of network functional is defined, but the user can overload the network function as they see fit.
    
    WormlikeNetwork: This takes the nodes of the graph as points in 3Space and the edges as lines between those points
    to define a skeleton of straight lines. G(network, point, params) is the set of all points that are precisely wormRadius or 
    closer to any part of the network connected lines.  Any singleton points that are not connected are ignored. Only edge connected points
    concern us. 
    
    It is easy to check if a point belongs to this set by computing the distance of a test point from every line segment connecting 
    each pair of connected nodes. Each pair of connected nodes contributes a spherocylinder (swept volume of radius wormRadius) to the 
    locus of G.

    If the point is within any of these swept volumes it is inside the set. Being inside the volume of one swept volume supercedes
    being on the surface or outside of another.  
    
    The user can specify using the packType parameter whether they want to return a set of points which is "inside" G(netwrok, points, params),
    "outside" or on the "surface"  of G. ( Surface is defined as a thin shell of thickness 2 * epsilon).
    
    All the usual rules such as envelopes and points to avoid work as usual. Thus one can define the G(Network, points, params) and 
    whatever other volumes may wish to be defined in the envelopes or points to avoid points.
    
    The default reference point of the network is the centre of mass of the *input* points.
    The default director is the z-axis.
    
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for NetworkPack object")
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, minDist, graph, points, params, packType, epsilon, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200,'envelope.xyz'), showBlockDirector=False):
        self.graph = graph
        self.points = points
        self.params = params
        self.packType = packType # inside, outside or surface
        self.epsilon = epsilon # tolerance on equality between two real values.  
        self.xMin, self.xMax, self.yMin, self.yMax, self.zMin, self.zMax = self.boundingRegion()
        BuildingBlock = NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector)
        directorStatus = self.computeDirectors()
        
        # extract and normalise the directors
        directors = [ dirn[1]/np.linalg.norm(dirn[1]) for dirn in directorStatus]
        return BuildingBlock, directors
    
    def generateBuildingBlockDirector(self):
        return np.array([0.0, 0.0, 1.0]) # z-axis.

    def generateBuildingBlockRefPoint(self):
        # default reference point is the mean of the input geometry. Can always redefine. 
        return np.mean(self.points, 0)
    
    def pickFirstPoints(self):
        self.nList = [self.pickRandomPointInDefinedSpace()]
        self.nAttempts = [0]
        return True
        
    def pickRandomPointInDefinedSpace(self):
        
        inBounds = False
        
        # keep looping until we pick a point within the bounding volume and determine if it is surface, inside or outside as desired.
        
        while inBounds==False:
            # pick a random point within the bounding cuboid of G(network, points, params) (computed at instantiation).
            pos = np.array([ rnd.uniform(self.xMin, self.xMax), rnd.uniform(self.yMin, self.yMax), rnd.uniform(self.zMin, self.zMax) ]) 
        
            # determine if the point is inside, outside or on the surface of the network. In process constructs the director 
            packType = self.checkPointStatus(pos)
            if packType==self.packType:
                inBounds = True
        
        return pos

    # returns the status and the director of the list of points going into the building block wrt to the network Graph member variable
    def computeDirectors(self):
        return [ self.checkPointStatus(point, returnDirector=True) for point in self.buildingBlockXYZ ] 
    
    def checkPointStatus(self, pos, returnDirector=False):
        # returns the status of the packtype of the particle with respect to G(network, points, params)
        # either 'surface', 'inside' or 'outside'.  
        
        # TO DO: Overide this function with your own functional. Must return a string "inside", "outside" or "surface", and implement
        # return director switch - allows you to optionally return the director of a point whose status you have checked)
        # outside points will always return the director as pointing to the closest point on the last edge in the network.

        # assume the point is outside G(network, p#oints, params)
        packType = 'outside'

        surfList = []
        
        # loop through the edges of the graph
        for edge in self.graph.edges:

            # don't check edges that self-connect nodes to themselves.            
            if edge[0]!=edge[1]:
            
                # compute the closest approach of the test point to the line segment between the points referred to by 
                # the nodes of the graph.  
                
                # if 
                (distSq, director) = cart.closestApproachPointToLineSegmentSquared(pos, self.points[edge[0]], self.points[edge[1]], returnVec=True)
                dist = np.sqrt(distSq)
                
                # if the dist to the line segment is within self.epsilon of parameter[0] - the worm radius - 
                # then test point is on the surface of this segment of the graph.  
                if np.abs(dist - self.params[0]) < self.epsilon:
                    # point is on surface of this segment, but could be inside any other segment so keep on trucking. Don't break
                    # the loop. 
                    packType='surface'
                    surfList += [(dist, cp.copy(director))] # generate a list of all the edges for which this point is a surf point
                else:
                    if np.abs(dist - self.params[0] - self.params[1]) < self.epsilon:
                        # point is on outer surface of a bilayer but could also be inside any other segment so keep on trucking. Don't break
                        # the loop. 
                        packType='bilayer'
                        surfList += [(dist, cp.copy(director))] # generate a list of all the edges for which this point is a surf point
                    else:
                        if np.abs(dist - self.params[0] + self.params[1]) < self.epsilon:
                            # point is on inner surface of a bilayer but could also be inside any other segment so keep on trucking. Don't break
                            # the loop. 
                            packType='bilayer'
                            surfList += [(dist, cp.copy(-1 * director))] # generate a list of all the edges for which this point is a surf point
                        else:
                            # if dist is < wormRadius then point belongs to the volume, so keep it. Break Now. Doesn't matter about
                            # any other point. 
                            if (dist < self.params[0]):
                                packType='inside'
                                break
                
        output = packType
        if returnDirector:
            # if we get to the end, point is a surface type and returnDirector is set then choose the director closest to the true surface
            if packType in ['surface', 'bilayer']:
                try:
                    minD = np.inf
                    for directorPair in surfList:
                        if directorPair[0]<minD:
                            minD = directorPair[0]
                            director = directorPair[1]
                        
                    #director = min(surfList)[1] # finds the director associated with the smallest dist.
                except ValueError:
                    print("Caught some weird bug") 
                      
            output = (packType, director)
        return output
        
  
    def boundingRegion(self):
        # returns outer limits on a cuboid shaped region of space which is guaranteed to 
        # be bigger than G(network, points, params).  Thinking of G(network, points, params)
        # as a closed level set. But user may want to think of G(network, points, params) as an open surface,
        # in which case it may run off to infinity even for a finite network. In that case, user must still define
        # region on interest in which to compute G(network, points, params). 
        #
        # TO DO Redefine bound for you purposes. 
        #
        # Here we find the connected nodes in the network with lowest and highest X, Y and Z coords, and increase them away from zero 
        # by param[0] - the worm radius.  We ignore singleton points in the points list that are not in the connected network. 
        # All points further away than are outside G(network. points, params)
        
        # set up the initial boundary as infinite box. 
        xMin = np.inf
        xMax = -np.inf
        yMin = np.inf
        yMax = -np.inf
        zMin = np.inf
        zMax = -np.inf
        
        wormRad = self.params[0]
        
        for edge in self.graph.edges():
            x1, y1, z1 = self.points[edge[0]]
            x2, y2, z2 = self.points[edge[1]]
            if  x1 - wormRad < xMin:
                xMin = x1 - wormRad 
            if  x2 - wormRad < xMin:
                xMin = x2 - wormRad 
            if  x1 + wormRad > xMax:
                xMax = x1 + wormRad 
            if  x2 + wormRad > xMax:
                xMax = x2 + wormRad 

            if  y1 - wormRad < yMin:
                yMin = y1- wormRad 
            if  y2 - wormRad < yMin:
                yMin = y2- wormRad 
            if  y1 + wormRad > yMax:
                yMax = y1 + wormRad 
            if  y2 + wormRad > yMax:
                yMax = y2 + wormRad   
                 
            if  z1 - wormRad < zMin:
                zMin = z1 - wormRad 
            if  z2 - wormRad < zMin:
                zMin = z2 - wormRad 
            if  z1 + wormRad > zMax:
                zMax = z1 + wormRad 
            if  z2 + wormRad > zMax:
                zMax = z2 + wormRad 
        
        # export the boundaries with a worm radius margin for kicks.
        return xMin - wormRad, xMax + wormRad, yMin - wormRad, yMax + wormRad, zMin- wormRad, zMax + wormRad 
        
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # create the NPack object.
    NetworkPackBBG = PackNetworkBBG('networkPack.txt')

    numPoints = 1000
    minDist = 2
    points = [ np.array([rnd.uniform(-40, 40), rnd.uniform(-40, 40), rnd.uniform(-40, 40)]) for _ in range(0, 20) ]
    graph = nx.Graph()
    graph.add_nodes_from( range(0, 20) )
    node1List = [ rnd.choice( range( 0, 20) ) for _ in range( 0, 10 ) ]
    node2List = [ rnd.choice( range( 0, 20) ) for _ in range( 0, 10 ) ]
    [ graph.add_edge(a, b) for a, b in zip(node1List, node2List) ] 
    params = [5]
    packType = 'surface'
    epsilon = 0.1

    # dump the network to file to help judge the output
    fIO.visualiseNetwork(graph, points, 'rawNetwork.xyz', particleName='P', dZ=1)
       
    # generate the building block
    NetworkPackBB, directors = NetworkPackBBG.generateBuildingBlock(numPoints, minDist, graph, points, params, packType, epsilon)
    fIO.saveXYZ(NetworkPackBB.blockXYZVals, 'C', "networkPack.xyz")
    
    dirPosns = [ pos + director for (pos, director) in zip(NetworkPackBB.blockXYZVals, directors) ]
    
    fIO.saveXYZ(dirPosns, 'O', "networkPackDirectorPoints.xyz")
    
    print("NetworkPack.xyz Example written")