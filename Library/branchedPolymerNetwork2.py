'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
from Library.VolumePackCuboid import VolumePackCuboidBBG as VPCBBG
from Builder.BuildingBlock import BuildingBlock
import sys
import numpy as np
import random as rnd
import copy as cp
import networkx as nx
from Utilities import fileIO as fIO
from Utilities import coordSystems as coords
from Utilities import cartesian as cart

class BranchedPolymerPackBBG(BBG):
    '''
    Overloads the building block generator class, and makes use of the random polymer pack generator,
    to generate a set of interweaving polymer chains. The directors of the first generation of chains can be 
    correlated within a specified deviation to a global director.
    
    The first generation of polymer chains can then be used as a seed points for branching polymer chains 
    which have a specified density and minimum interior branch angle.  Each of the subsequent unimers branches
    may then be used as a seed points for additional unimers up to a specified number of generations of branching. 
    
    Every unimer that is generated will have a envelope which is a cone. The same shape cone will be used 
    for all unimers - they are related molecules afterall - but the orientation will be selected based on the degree
    of correlation selected.  Prior to computing each new unimer any of the over all set of points that are in the cone
    will be isolated as points to Avoid for that unimer. In this way no unimer thread should intersect any other. 
    
    The self avoidance can be turned off. It turns out the polymers miss each other most the time anyway. Can live with 
    the odd clash, if we're only making cartoons.   
    
    A networkx graph of the polymer network is created whose nodes are the indices of the points in output xyzlist
    and whose edges are the points in that list which are connected by a bond. 
    
    There are many additional parameters: 
    the alpha and beta ranges for the internal angles of the polymer.
    the networkMinDist is the distance between the starting points of the polymer chains
    MonomerMinDist is as close as two chains can come. 
    BondLength is the distance between successive unimers along the chain. Bondlength must be longer than minDist
    or the chain self rejects. (sorry).
    Minpoints per unimer is the smallest number of points allowed on a chain. typically the chain length is 
    determined by the distance to the outer box, and then 50 % of the time the polymers are 
    randomly truncated by a randomfactor between half and 1. MaxNumAttempts is a bailout counter for some
    while loops which repeat various unimer constructions. Oftimes the starting point is just in a bad position. 
     
    Uggh. I'm bored. figure out the rest of the params yourself. They're pretty self explanatory.
    
    
    '''
    def __init__(self, filename):
        BBG.__init__(self, filename)

    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        BBG.initialiseParameters(self) 
        
        self.particleName = self.getParam('particleName')

        self.RPPBBG = RPPBBG(self.paramFilename)
        self.VPCBBG = VPCBBG(self.paramFilename)
        
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for branched polymer network object")
            sys.exit()

    def generateBuildingBlock( self, 
                               numGen1Unimers,
                               x1,
                               x2,
                               y1,
                               y2,
                               z1,
                               z2,
                               networkMinDist,
                               alpha1, 
                               alpha2,
                               beta1,
                               beta2,
                               monomerMinDist,
                               bondLength,
                               minNumPointsUnimer,
                               maxNumAttempts,
                               selfAvoid,
                               coneAngle,
                               directorMaxPitch,
                               globalDirector,
                               numBranchGenerations,
                               branchDensity,
                               minBranchAngle,
                               angularRange=['None'],
                               pointsToAvoid=[],
                               visualiseEnvelope=(0,20),
                               envelopeList=["None"],
                               SpaceCurveTransform=True):

        numPoints = 0 # computed on the fly in this case
        self.networkGraph = nx.Graph() 
        self.numGen1Unimers = numGen1Unimers 
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.z1 = z1
        self.z2 = z2
        self.networkMinDist = networkMinDist  
        self.monomerMinDist = monomerMinDist
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.beta1 = beta1
        self.beta2 = beta2
        self.bondLength = bondLength
        self.minNumPointsUnimer = minNumPointsUnimer
        self.maxNumAttempts = maxNumAttempts
        self.selfAvoid = selfAvoid
        self.coneAngle = coneAngle * np.pi/180.0
        self.directorMaxPitch = directorMaxPitch * np.pi/180.0
        self.globalDirector = globalDirector
        self.numBranchGenerations = numBranchGenerations
        self.branchDensity = branchDensity
        self.minBranchAngle = minBranchAngle * np.pi/180.0
        self.pointsToAvoid = pointsToAvoid
        self.visualiseEnvelope = visualiseEnvelope
        self.envelopeList = envelopeList + ['cuboid ' + str(x1) + ' ' + str(x2) + ' ' + str(y1) + ' ' + str(y2) + ' ' + str(z1) + ' ' + str(z2)]
        self.angularRange = angularRange
        self.blockRefPoint = self.generateBuildingBlockRefPoint()
        self.SpaceCurveTransform=SpaceCurveTransform
         
        return BBG.generateBuildingBlock(self, numPoints, networkMinDist, envelopeList=self.envelopeList, visualiseEnvelope=visualiseEnvelope, pointsToAvoid=pointsToAvoid)

    def generateBuildingBlockDirector(self):
        return  np.array([0.0, 0.0, 1.0])
    
    def generateBuildingBlockRefPoint(self):
        return np.array( [ ( self.x1 + self.x2 ) / 2.0, 
                           ( self.y1 + self.y2 ) / 2.0, 
                           ( self.z1 + self.z2 ) / 2.0 ] )

    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints
    
    # over load the final function to get a building block, so that the networkGraph can be passed back as well. 
    def getBuildingBlock(self):
        # access function that returns a building block based on the backbone
        return (BuildingBlock(self.buildingBlockXYZ, 
                             self.blockNames, 
                             self.blockConnectors, 
                             self.blockRefPoint,
                             self.blockDirectorHat), self.networkGraph)

    def generateBuildingBlockXYZ(self):
        # first pick n random points in the specified envelope within a cubic region of space.
        # These are the unimer start points. Leave a margin at the edge to try to eliminate dodgy scenarios
        startPoints = self.VPCBBG.generateBuildingBlock(   self.numGen1Unimers, 
                                                           self.x1 + self.minNumPointsUnimer * self.bondLength, 
                                                           self.x2 - self.minNumPointsUnimer * self.bondLength, 
                                                           self.y1 + self.minNumPointsUnimer * self.bondLength, 
                                                           self.y2 - self.minNumPointsUnimer * self.bondLength, 
                                                           self.z1 + self.minNumPointsUnimer * self.bondLength, 
                                                           self.z2 - self.minNumPointsUnimer * self.bondLength, 
                                                           self.networkMinDist, 
                                                           envelopeList=self.envelopeList,
                                                           pointsToAvoid=self.pointsToAvoid)

        if self.dumpInterimFiles==1:
            fIO.saveXYZ(startPoints.blockXYZVals, 'C', 'pointsInCube.xyz')

                
        # for each point generate a direction relative to the globalDirector in the allowed pitch range.
        # azimuthally can choose anything, but elevation wise we are restrained to be within directorMaxPitch 
        # angle of the globalDirector. Choose a direction that yields at least minNumPointsUnimer
        gen1UnimerDirectors = [ self.genDirector(point, self.minNumPointsUnimer, self.globalDirector, self.directorMaxPitch, self.maxNumAttempts)  for point in startPoints.blockXYZVals ] 
        
        # set up output arrays
        xyzPoints = []
        UnimerIndices = []
        unimersToBranchFrom = []

        currentUnimer = 0
        # Now generate the gen1 Unimers
        for point, director in zip(startPoints.blockXYZVals, gen1UnimerDirectors):
            print("Generating Unimer: ", currentUnimer + 1, " of ", self.numGen1Unimers)
            
            # generate the next list of xyz values
            if self.selfAvoid:
                # director[0] is the vector, director[1] is the number of points
                unimerXYZVals  = self. generateUnimer(point, director[0], director[1],  xyzPoints + self.pointsToAvoid, self.maxNumAttempts)
            else:
                unimerXYZVals  = self. generateUnimer(point, director[0], director[1],  self.pointsToAvoid, self.maxNumAttempts)
            
            # make a tuple to note the start and end index of the new unimer in the points array.
            # add that tuple to a list of indices
            UnimerIndices += [ (len(xyzPoints), len(xyzPoints) + len(unimerXYZVals)) ]
            
            # add the xyz points of the unimer to the growing list of xyzValues
            xyzPoints += cp.copy(unimerXYZVals) 

            # make a note of the index of the latest tuple pair added.
            # Each generation only adds to the last worms created 
            unimersToBranchFrom += [len(UnimerIndices)-1]

            currentUnimer += 1
        
        print("First Generation Unimers Complete.")
        
        if self.dumpInterimFiles==1:
            fIO.saveXYZ(xyzPoints, 'Ca', 'FirstGenUnimers.xyz')

        # loop through the number of generations we want to have        
        for genNum in range(0, self.numBranchGenerations):
            print("Now creating unimers for generation ", genNum + 1, " of ", self.numBranchGenerations)

            # create a new generation of unimers that will branch off the specified unimers.
            # return a full list of all the points, the indices of each unimer in points
            # and a list of new unimers from which we can branch from on the next generation 
            xyzPoints, UnimerIndices, unimersToBranchFrom = self.generateNextUnimerGeneration(xyzPoints, UnimerIndices, unimersToBranchFrom)

        # we have a list of tuples which contain the start and end indices of each unimer.
        # Add this information to the graph network to complete the graphs.
        # The nodes that connect the daughter and parent unimers together are already in there. 
        self.addNodesAndEdgesToGraphFromUnimerIndices(UnimerIndices)
            
        return xyzPoints
    
    # Each entry in unimer indices indicates a connected set of nodes.
    # Add those nodes to the graph and then edges between them.
    def addNodesAndEdgesToGraphFromUnimerIndices(self, UnimerIndices):
        for segment in UnimerIndices:
            
            # add the first node of the segment if it is not there already
            if not self.networkGraph.has_node(segment[0]):
                self.networkGraph.add_node(segment[0])
            
            # for each subsequent node, add it and connect it to the previous node in the segment
            for node in range(segment[0] + 1, segment[1] ):
                if not self.networkGraph.has_node(node):
                    self.networkGraph.add_node(node)
                
                # add the edge connecting the previous node to the current node
                self.networkGraph.add_edge(node - 1, node)

    # generates a director for a given point. Makes sure that there is space for at least N bond lengths 
    # from the point to the edge of the box for that director. Also ensure the picked director is within 
    # angle radians of gDirector.
    def genDirector(self, point, minNumPoints, gDirector, angle, maxNumAttempts, within=True):
        
        numPoints = 0
        maxNumPoints = 0
        numLoops = 0
        # loop until we exceed the minimum number of points requested or we have tried 10 times
        while numPoints<minNumPoints and numLoops<maxNumAttempts:
            if numLoops>0:
                print("Finding alternative director. Attempt: ", numLoops + 1, "of ", maxNumAttempts)
            if within:
                # for choosing first generation unimers (must be within angle of gDirector)
                theta, phi = coords.pickRandomPointOnUnitSphereInAngRange(np.pi/2 - angle, np.pi/2, -np.pi, np.pi)
            else:
                # for choosing nth generation unimers (must be between angle and 2 angle of gdirector)
                theta, phi = coords.pickRandomPointOnUnitSphereInAngRange(np.pi/2 - 2 * angle, np.pi/2 - angle, -np.pi, np.pi)
            
            # convert theta, phi pair to xyz vectors
            directorHat = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi]))
            
            #  rotate director to take into account orientation of globalDirector
            directorHatRot = coords.transformFromBlockFrameToLabFrame( gDirector, np.array([0.0, 0.0, 0.0]), 0.0, 
                                                                       np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), [directorHat])[0]
    
            # calculate distance along director vector to the outer edge of the box.
            dist2Box = self.dist2Box(point, directorHatRot)
                
            # compute number of points we can fit in there assuming straight line to edge of box
            numPoints = int(dist2Box/self.bondLength)
            
            # make a note of the direction that gives the largest number of points so far
            if numPoints>maxNumPoints:
                maxDirector=cp.copy(directorHatRot)
                maxNumPoints = numPoints
            
            if numPoints<minNumPoints:
                print("low num points")
            
            numLoops+=1

        # if we have exited the loop then check to see if we have succeeded.
        # if we have failed then return the best director we encountered.
        if numPoints<minNumPoints:
            print("Warning: Director choice resulted in short polymer")
            numPoints = maxNumPoints
            directorHatRot = maxDirector

        return (directorHatRot, numPoints)

    def generateNextUnimerGeneration(self, points, indices, unimersToGrowFrom):
        # The points for the entire network are contained in flat list of points. 
        # Each unimer is identified by a tuple of indices referring to the start and end points in the array.
        # The list of tuples is stored in indices. 
        # Unimers to Grow from is a list of unimers i.e. their position in the index list, from which to grow
        # a new unimer.
        
        # create an array to keep track of the index pairs of each new unimer we add them.
        # This will form the list to grow from for the next generation.
        newUnimersToGrowFrom = []        
        
        # loop through the unimers to grow from.
        for unimerNum, unimer in enumerate(unimersToGrowFrom):
            
                print("Generating Unimer: ", unimerNum + 1, " of ", len(unimersToGrowFrom))
            
                # get the indices of the current unimer from which we are going to branch
                startIndex = indices[unimer][0]
                endIndex = indices[unimer][1] 
                length = endIndex - startIndex

                # compute how many branch points to add               
                numBranchPoints = int(float(length) * self.branchDensity)
                
                # generate a list of indexes where the branches will happen
                seedPoints = []
                for _ in range(0, numBranchPoints):
                    newSeedPoint = rnd.randint(startIndex + 2, endIndex - 2) # don't use the first few points in a unimer - makes for crappy looking networks. 
                    while newSeedPoint in seedPoints: # loop until we don't pick the same index multiple times
                        newSeedPoint = rnd.randint(startIndex + 2, endIndex - 2)
                    seedPoints += [ newSeedPoint ]
            
                # loop through the selected seed points 
                for seedPoint in seedPoints:
                    # generate a director that is not within branch angle of the local 
                    # seed unimer
                    director = self.genDirector(points[seedPoint], 
                                                self.minNumPointsUnimer, 
                                                points[seedPoint + 1] - points[seedPoint], 
                                                self.minBranchAngle, self.maxNumAttempts, within=False)

                    
                    if self.selfAvoid:
                        newUnimerXYZ = self.generateUnimer(points[seedPoint], director[0], director[1], points + self.pointsToAvoid, self.maxNumAttempts)
                    else:    
                        newUnimerXYZ = self.generateUnimer(points[seedPoint], director[0], director[1], self.pointsToAvoid, self.maxNumAttempts)

                    # Add the current new unimer's point indices to the list of unimers start and stop indices - we subtract 1
                    # because the seed point is already in the parent unimer. We are only adding the new points.
                    indices += [ (len(points), len(points) + len(newUnimerXYZ) - 1) ]

                    # At this point this is the only time in the program we know the index values of 
                    # both the parent and daughter unimers and where they are connected.
                    # so add this information to the graph.
                    
                    # add the current seed point as a node if it isn't already there
                    if not self.networkGraph.has_node(seedPoint):
                        self.networkGraph.add_node(seedPoint)  
                    
                    # add the start index of the most recent unimer added to the unimer index list
                    if not self.networkGraph.has_node(indices[-1][0]):
                        self.networkGraph.add_node(indices[-1][0])
                    
                    # connect those two nodes in the graph
                    self.networkGraph.add_edge(seedPoint, indices[-1][0])

                    # add the xyz points of the unimer to the growing list of xyzValues
                    # only add the new points - don't add the seed point which would be redundant.
                    points += cp.copy(newUnimerXYZ[1:]) 
                    
                    # add the index number of the latest unimer to the newUnimersToGrowFrom list for the next gen
                    # subtract 1 because it is zero based, where the length of the array is not.   
                    newUnimersToGrowFrom += [len(indices) - 1]    

        return points, indices, newUnimersToGrowFrom

    # generates a single unimer which starts at point point, in direction of director
    # avoiding all the points in xyzList, and uses global params in the class to
    # define the unimer details.     
    def generateUnimer(self, point, director, numPoints, pointsToAvoid, maxNumAttempts):

        # half the time, scale the number of points by a random factor between 0.5 and 1.0
        # otherwise all the unimers will leave the box
        if rnd.uniform(0.0, 1.0) > 1.0:
            numPoints *= rnd.uniform(0.5, 1.0)
        
        # make the numPoints an integer
        numPoints = int(numPoints)

        # develop the frustum definition in the block frame of the unimer 
        z1 = -self.bondLength
        r1 = self.bondLength
        z2 = numPoints * self.bondLength
        r2 = z2 * np.tan(self.coneAngle/2.0)
        
        localEnvelopeList = ["frustum " + str(z1) + ' ' + str(r1) + ' ' + str(z2) + ' ' + str(r2)]  
        
        if numPoints==0:
            print("Zero Length Polymer")
        numPointsOutside = numPoints
        numLoops = 0
        # loop until we create a unimer that satisfies constraints. Can set maxNumAttempts depending on how bothered we are about constraints.
        while numPointsOutside > 0 and numLoops<maxNumAttempts:
            if numLoops>0:
                print("Re-generating Unimer - exceeds global envelope. Attempt: ", numLoops + 1, "of ", maxNumAttempts)
            # generate the unimer in it's own block frame starting at origin. Don't worry about outer envelope or other points at this stage
            unimer = self.RPPBBG.generateBuildingBlock( numPoints, 
                                                        np.array([0.0, 0.0, 0.0]),
                                                        self.alpha1, 
                                                        self.alpha2,
                                                        self.beta1,
                                                        self.beta2,
                                                        self.monomerMinDist,
                                                        self.bondLength, 
                                                        pointsToAvoid=[],
                                                        envelopeList=localEnvelopeList,
                                                        SpaceCurveTransform=self.SpaceCurveTransform)
                                                        #visualiseEnvelope=(100000,200,'envelopeUnimer.xyz'))

            if self.dumpInterimFiles==1:
                fIO.saveXYZ(unimer.blockXYZVals, 'Ca', 'preRotUnimer.xyz')
    
            # now rotate and translate the unimer xyzVals to their final point
            rotXYZVals =  coords.transformFromBlockFrameToLabFrame(director, point, 0.0, 
                                                                  np.array([0.0, 0.0, 1.0]), 
                                                                  np.array([0.0, 0.0, 0.0]),
                                                                  unimer.blockXYZVals)  
            
            # Now check that the unimer is entirely within the outer box, and any global specified user envelope,
            # Also make sure that the unimer is networkMinDist points away from existing pointsToAvoid.
            # Won't self check the unimer against itself only the points to avoid. Allows to have two minDists.
            # One for the monomer distance along the backbone, and one for the distance between unimers (2 * worm tube radius).
            # enables the unimer bondlength to be small which makes for smoother unimer chains.
            numPointsOutside = len(self.checkEnvelope(rotXYZVals, ignorePTA=False, pointsToAvoid=pointsToAvoid))
            
            # increment loop counter
            numLoops += 1
        
        return rotXYZVals
    
#   def choosePointsToAvoid(self, points, Z1, R1, Z2, R2):
#       # go through points and keep any points that are inside the frustum
#       outPoints = []
#       for point in points:
#           if coords.checkPointInFrustum(point, Z1, R1, Z2, R2, 0): 
#               outPoints += [ point ]
#       return outPoints
            
    def dist2Box(self, point, director):
        # Find the distance from the point to the outer box in the direction of director
        
        # compute the line parameter where the given line intersects each of the six planes of the box.
        paramVals = []
        paramVals += [ cart.distanceOfLineToAPlane(point, director, np.array([ self.x1, 0, 0]), np.array([1.0, 0.0, 0.0]) ) ]
        paramVals += [ cart.distanceOfLineToAPlane(point, director, np.array([ self.x2, 0, 0]), np.array([1.0, 0.0, 0.0]) ) ]
        paramVals += [ cart.distanceOfLineToAPlane(point, director, np.array([ 0, self.y1, 0]), np.array([0.0, 1.0, 0.0]) ) ]
        paramVals += [ cart.distanceOfLineToAPlane(point, director, np.array([ 0, self.y2, 0]), np.array([0.0, 1.0, 0.0]) ) ]
        paramVals += [ cart.distanceOfLineToAPlane(point, director, np.array([ 0, 0, self.z1]), np.array([0.0, 0.0, 1.0]) ) ]
        paramVals += [ cart.distanceOfLineToAPlane(point, director, np.array([ 0, 0, self.z2]), np.array([0.0, 0.0, 1.0]) ) ]

        # choose the smallest +ve distance. (selects the plane in the director of director)
        # when the line is parallel with the plane we're not bothered about that plane.
        # when the line goes through the corner of the box or two of the planes, it's the same distance 
        # so we don't care. Automatically handles these cases.
        return min([a for a in paramVals if a > 0.0]) 
        
if __name__ == '__main__':
   
    # create the generator
    BPPBBG = BranchedPolymerPackBBG('branchedPolymerNetwork.txt')
    numGen1Unimers = 4
    x1 = -100
    x2 = 100
    y1 = -100
    y2 = 100
    z1 = -100
    z2 = 100    
    networkMinDist = 50.0
    bondLength = 1.0
    alpha1 = 45
    alpha2 = 65
    beta1 = 150 #170
    beta2 = 170 #180

    minNumPointsUnimer = 5
    maxNumAttempts = 100
    
    coneAngle = 45 
    directorMaxPitch = 20 
    globalDirector = np.array([1.0, 1.0, 0.0])
    
    numBranchGenerations = 1
    minBranchAngle = 25 
    branchDensity = 0.05

    pointsToAvoid=[]
    visualiseEnvelope=(0,20)
    envelopeList=["None"]
    selfAvoid = False
    
    # generate a set of points and build a constrained polymer pack network between them. 
    # Returns a list of points and a networkx graph showing connections between those points
    BranchedPolymerPackBB, networkGraph = BPPBBG.generateBuildingBlock(  numGen1Unimers, 
                                                                           x1, x2,
                                                                           y1, y2,
                                                                           z1, z2, 
                                                                           networkMinDist,
                                                                           alpha1, alpha2,
                                                                           beta1, beta2,
                                                                           0.9 * bondLength, 
                                                                           bondLength,
                                                                           minNumPointsUnimer,
                                                                           maxNumAttempts,
                                                                           selfAvoid,
                                                                           coneAngle,
                                                                           directorMaxPitch,
                                                                           globalDirector,
                                                                           numBranchGenerations,
                                                                           branchDensity,
                                                                           minBranchAngle,
                                                                           visualiseEnvelope=(10000, x2 - x1, 'outerBoxEnv.xyz'),
                                                                           SpaceCurveTransform=False)

    fIO.saveXYZ(BranchedPolymerPackBB.blockXYZVals, 'C', "branchedPolymerNetwork.xyz")
    print("branched Polymer Done")