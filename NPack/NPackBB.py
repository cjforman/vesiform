'''
Created on 14 Dec 2017

@author: chris
'''
import sys
import numpy as np
import random as rnd
import Utilities.cartesian as cart
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG

class NPackBB(BBG):
    ''' This class inherits the building block generator class and adds to it the ability 
    to generate random structures using the core NPack algorithm and checking routine. 
    The generateBuildingXYZ function is then able to randomly select a list of numpy xyz vectors 
    that satisfy the applied constraints. Some example classes are then provided which uses this capability
    showing how to use the class - 3DSpacePack, Flat, Cylinder, Sphere, Ellipsoids and Chain.
    
    Static values for a given class of object are provided in a textfile with examples in the 
    library. Parameter handling is automatic since the keyProc base class is inherited via the 
    BuildingBlockGenerator class. Just add parameters to the initialise parameters function. This
    system automatically lists necessary parameters that are not included in the text file at 
    run time. 
    
    Dynamic values such as numPoints, startPosition, global orientation are are provided at
    run time as for all BuildingBlockGenerator objects
    
    NPack Parameters that must be defined in the params file are:
    
    '''
    
    def __init__(self, paramFilename):
        # initialise the parameter dictionary
        BBG.__init__(self, paramFilename)
  
    def initialiseParameters(self):
        # copy across necessary parameters common to all possible structures.
        self.maxNumAttempts = self.getParam('maxNumAttempts')
        self.maxNumOperations = self.getParam('maxNumOperations')
        self.maxNumPrunesBeforeChop = self.getParam('maxNumPrunesBeforeChop')
        self.maxNumChops = self.getParam('maxNumChops')
        self.maxLivesPerNode = self.getParam('maxLivesPerNode')
        self.maxLivesPerNewNode = self.getParam('maxLivesPerNewNode') 
        self.particleName = self.getParam('particleName')
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')
        self.verbose = self.getParam('verbose')
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for NPack object")
            sys.exit()
    
    def generateBuildingBlock(self, numPoints, minDist, envelopeList='[None]', pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False, defaultBlockRefPoint=None):
        self.numPoints = numPoints
        self.minDist = minDist
        return BBG.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=showBlockDirector, defaultBlockRefPoint=defaultBlockRefPoint)
    
    def generateBuildingBlockXYZ(self):
        # this is the main new bit for this class relative to the old NPack class.
        # Calls the generateNConstrainedRandomPoints function to populate the
        # BuildingBlockXYZ variable allowing the building block generator sub class
        # to process the resulting xyz lists.
        # This is inefficient since the same object now has two copies of the same list.
        # (four in fact because there is BuildingBlockXYZ, NList, longestNList and BuildingBlockFinalXYZ)
        # Without redesigning the whole thing can't get round this so tough luck.
        self.generateNConstrainedRandomPoints()
        return self.getNList()

    def generateBuildingBlockDirector(self):
        return np.array([0, 0, 1])
    
    def generateBuildingBlockRefPoint(self):
        return cart.getCentreOfMass(self.buildingBlockXYZ)
    
    def generateBuildingBlockNames(self):
        return [self.particleName] * self.numPoints 
    
    def generateNConstrainedRandomPoints(self):
        ''' Returns True if a list of numPoints numpy vectors are randomly found that obey the user supplied conditions. 
        Returns False if the number of attempts to find a list exceeds maxNumAttempts
        '''
        print("Starting packing process")
        
        # initialise random number generator
        rnd.seed()

        # initialise flags and counters
        processSucceeded = False   # Assume failure. This is the return variable
        numAttempts = 0            # number of times the algorithm has tried to succeed
        
        # keep looping until the number of attempts exceeds the max number allowed or the process succeeds 
        while numAttempts < self.maxNumAttempts and processSucceeded == False:
            
            # increment the number of attempts.
            numAttempts += 1 

            # attempt to seed the list
            listSeeded = self.pickFirstPoints()
            
            # initialise longestNList
            self.longestNList = self.nList
            
            # if the list seeded successfully then attempt to find N points from the current seed
            if listSeeded:
                # set the number of seed points - this may be varied by the users overloaded function
                self.numSeedPoints = len(self.nList)
                
                # given the current seed, attempt to find a list of n points that obeys the conditions.
                # Reports on the success or failure of the endeavour.
                processSucceeded = self.attemptToFindNPointsFromCurrentSeed()
            
        if processSucceeded == True:
            print("Packing Process Successful")
        else:
            print("Process failed maxNumAttempts times. Consider relaxing constraints or increasing maxNumAttempts.")
            
        return processSucceeded

    def attemptToFindNPointsFromCurrentSeed(self):
        ''' Returns True if a seedlist can be grown into an array of numPoint random numpy vectors that 
        the obey the user overloaded conditions. 
        
        Returns false if 
            1) the current list can neither be added to nor pruned from (tried to delete root node)
        or  2) the total number of attempts to operate on the list exceeds a certain threshold.
        or  3) the number of times a list was chopped exceeds maxNumChops (means pruning was not leading to additional growth)
        
        Number 2) is a catch all to ensure that we don't get stuck in an infinite loop that the user never 
        foresaw with their conditions.  Can set this number super high if we want to give a larger chance 
        for success conditions to emerge.  Down to the patience of the user!
        
        '''
        processFailed = False # assume attempt will be successful 
        
        # initialise a dynamic array that keeps track of how many times the
        # algorithm looped while each node was the end point of the list
        self.nAttempts = [0] * self.numSeedPoints
            
        # set the counter for the number of times we had to perform a chop
        numChops = 0
            
        # initialise longest length record 
        longestLengthYet = self.numSeedPoints

        # highest number of prunes before growth occured
        largestPruneBeforeGrowthYet = 0
                
        # initialise counter for num of prunes that have occurred since the last time a new max length was encountered
        numTimesPrunedSinceLastMaxLengthGrowth = 0
                
        # initialise counter for counting how many operations we have attempted on the list in total
        numOperations = 0
        
        # Each loop either adds a node or prunes a node (or group of nodes) from the list. 
        # We loop until:
        #    1) success condition arises (list hits target length)
        #    2) processFailed = True
        
        # Process fails when:
        #    1) we can neither add nor prune a node (typically a seed is jammed up against the edge of an excluded region)  
        #    2) the total number of attempts to operate on the list exceeds a threshold. (fail safe to avoid infinite loops)
        #    3) we had to perform too many chops. Chops occur when many small prunes failed to achieve list growth. so problem is further down the list.
        #    4) an attempt was made to chop a short list (short = 2 * seed array length). Chopping short lists suggests that the seed is bad.
        
        while len(self.nList) < self.numPoints and processFailed==False: 
        
            # each loop constitutes one operation
            numOperations += 1
            
            if (numOperations % 100) == 0:
                print("num Operations: ", numOperations)
        
            if self.verbose==2:
                print("current Length: ", len(self.nList), ", numOperations: ", numOperations, ", numTimesPruned: ", numTimesPrunedSinceLastMaxLengthGrowth,  ", nAttempts: ", self.nAttempts)
        
            # if we hit the maximum number of operations then bail out the loop
            if numOperations>=self.maxNumOperations:
                print("Max Number of allowed operations exceeded; current seed has failed.")
                processFailed = True
        
            if not processFailed:
                # can only attempt to add a node to the current node up to maxLivesPerNode times. 
                nodeAdded = False                        
                while self.nAttempts[-1]<self.maxLivesPerNode and nodeAdded == False: 
                    
                    # add one to the number of times we tried to add a node to the current end node
                    self.nAttempts[-1] += 1
    
                    # attempt to add a node - returns true if we successfully find a node
                    nodeAdded = self.attemptToAddANode()
    
                # check to see if we added a node.
                if nodeAdded == True:
                    
                    # check to see if the list is the longest yet produced.
                    if len(self.nList)>longestLengthYet:
                        
                        # make back up of list
                        self.longestNList = self.nList[:]
                        
                        # Record the longest length so far.
                        longestLengthYet = len(self.nList)
                        
                        # reset the numTimesPrunedSinceLastMaxLengthGrowth 
                        numTimesPrunedSinceLastMaxLengthGrowth = 0
    
                        # send this rare but informative report to keep user aware that stuff is going on.
                        print("New Longest Length: ", longestLengthYet, " out of ", self.numPoints)
                
                else: # unable to add a node
                    
                    # attempt to prune the current list
                    nodesPruned = self.pruneNodes()
                                
                    # check to see if we successfully pruned the tree
                    if nodesPruned:
    
                        # increment count of number of times pruned 
                        numTimesPrunedSinceLastMaxLengthGrowth += 1
                        
                        if numTimesPrunedSinceLastMaxLengthGrowth>largestPruneBeforeGrowthYet:
                            largestPruneBeforeGrowthYet = numTimesPrunedSinceLastMaxLengthGrowth
                            print("New Largest Prune: ", largestPruneBeforeGrowthYet, "chop occurs at: ", self.maxNumPrunesBeforeChop, " Num Chops: ", numChops)  
    
                        # if we have successfully pruned the current list without causing new growth for ages then perform a chop.
                        # This has the effect of moving the list a long way in search space 
                        if numTimesPrunedSinceLastMaxLengthGrowth >= self.maxNumPrunesBeforeChop:
                            
                            # only perform a chop if the list is long enough to warrant it
                            if len(self.nList) > 2 * self.numSeedPoints:
                                numChops += 1
                                
                                if numChops>=self.maxNumChops:
                                    print("Max Num Chops reached. Current Seed has failed.")
                                    processFailed = True
    
                                if not processFailed:
                                    # pick a chop point which is a random point in the list between 
                                    # the seedArray and half way down the list at the current length
                                    chopPoint = rnd.randint(self.numSeedPoints, np.floor(len(self.nList)/2))
                                
                                    del self.nList[chopPoint:]
                                    del self.nAttempts[chopPoint:]
                                    
                                    # send this rare but informative message by way of a heart beat to the user
                                    print("Num Chops: ", numChops, " out of: ", self.maxNumChops)
                                    
                                    # reset the pruning counters
                                    numTimesPrunedSinceLastMaxLengthGrowth = 0
                                    largestPruneBeforeGrowthYet = 0
    
                            else:
                                print("Attempt to perform a chop on short list. Current seed has failed.")
                                processFailed = True
                    else:
                        # both nodedAdded and nodePruned were False
                        print("Unable to either add or a prune a node therefore current seed has failed.")
                        processFailed = True

        if processFailed:
            print("Current Attempt Has failed.")
            print("longest Length Reached: ", longestLengthYet, " number of operations: ", numOperations, " numChops: ", numChops, " highest number of Prunes: ", largestPruneBeforeGrowthYet) 
                        
        # need to report success so report inverse of failure
        return not processFailed

    def attemptToAddANode(self):                        
        ''' Returns True if a good node was found for the current list.
        Returns False if we had to try more than maxLivesPerNewNode times.  
        
        Each time a new node is picked in the user defined space, 
        it is checked to make sure it is within the user defined bounds within that space.
        
        If the node is in bounds in the user space then it is converted into a xyz vector in the lab space.   
        The new xyz vector is then checked against all the other members of the list.
        
        if the point is still good at this point then we add it to the list and create a new entry in the 
        numAttempts array with a value of zero.
        
        return whether or not we found a good point. '''
        
        goodPoint = False # assume failure

        numLoops = 0 # count number of times we try to find a new node 
                            
        # loop while trying to find a good new point maxLivesPerNewNode time
        while ((goodPoint==False) and (numLoops < self.maxLivesPerNewNode)): 
            
            # increment counter
            numLoops += 1
            
            # pick a new point in the defined sub space. 
            # Result is written into self.newPoint if successful
            
            newPoint = self.pickRandomPointInDefinedSpace()
            
            # check the new Point is within the acceptable part of the defined space.
            goodPoint = self.checkPointInBounds(newPoint)
                            
            # if we found a good new point then perform a global check - 
            # only perform a global check if the point is in bounds.               
            if goodPoint==True:
                                  
                # convert from point in user space to XYZ using user defined function.
                # useful for generating chains
                newPointInXYZ = self.convertPointToXYZ(newPoint)
            
                # check to see if the new point satisfies the global constraints
                goodPoint = self.checkPointAgainstList(newPointInXYZ)
                
        # if we found a good point that meets all the criterion then 
        # add it to the list
        if goodPoint==True:
            # append the latest point to the list
            self.nList.append(newPointInXYZ)
                            
            # the new point is now the last entry on the chain and has never been added to. 
            # Initialise the last position counter array for the newest point to zero
            self.nAttempts.append(0)

        return goodPoint


    def pruneNodes(self):

        # Starting from the last node, which only gets called for pruning when it's age hits self.maxLivesPerNode,
        # delete all the adjacent nodes whose age has reached self.maxLivesPerNode.
        
        pruneNodesSuccess = False # default return value assumes failure was encountered 
        continueToPrune = True # flag to indicate that pruning should continue 
        numNodesPruned = 0
         
        # keep pruning until told to stop
        while (continueToPrune):
            
            # if the end point is old then perform pruning operations.
            if self.nAttempts[-1] >= self.maxLivesPerNode:
                if len(self.nList)==self.numSeedPoints:
                    print("Attempt to delete the seed points. Current Seed is not viable.")
                    pruneNodesSuccess = False
                    continueToPrune = False    
                else:
                    pruneNodesSuccess = True
                    del self.nList[-1]
                    del self.nAttempts[-1]
                    numNodesPruned += 1    
            else:
                # end point is no longer old enough to prune
                continueToPrune = False
         
        return pruneNodesSuccess
        
    def pickFirstPoints(self):
        ''' Returns true if the nList is successfully initialised so that the first few points in the list comply
        with the constraints. Returns false if some good first points cannot be found.'''
        # TO DO: Overload this function to initialise the list and the lives list
        self.nList = [np.array([rnd.uniform(-1.0, 1.0), rnd.uniform(-1.0, 1.0), rnd.uniform(-1.0, 1.0)])]
        self.nAttempts = [0]
        return True
    
    def pickRandomPointInDefinedSpace(self):
        ''' Returns true if a point is successfully randomly selected in the allowed space 
        and false otherwise. The point is stored in a member variable and it's format 
        is defined by the user but is wlays a numpy 3-array.'''
        # TO DO: Overload this function to randomly pick a point in the defined space
        # default is a point in cube between (-boxsize/2, boxsize/2) in x, y and z.
        boxsize = 1.3 * self.minDist * np.power(self.numPoints,1.0/3.0)
        
        return np.array([rnd.uniform(-boxsize/2.0, boxsize/2.0), 
                         rnd.uniform(-boxsize/2.0, boxsize/2.0), 
                         rnd.uniform(-boxsize/2.0, boxsize/2.0)]) 
        
    def checkPointAgainstList(self, pos):
        ''' Function returns true if the test position is greater than minDist distance from all the other points.
            Assumes that pos is not in the list already. Can over load to do any global checking you like.'''
        # assume position is good
        goodPos = True

        for zPos in self.nList:
            # if new position occurs within minDist of any of the existing positions in the list 
            # then throw a wobbly.
            if np.linalg.norm((zPos - pos)) <= self.minDist:
                goodPos = False
                if self.verbose==1: 
                    print("Packing violation")

        return goodPos

    def convertPointToXYZ(self, newPointInSpace):
        ''' function takes a point generated by the pickRandomPointInDefinedSpace function
        and converts it to the xyz position that it needs to be in for the given space. This function is 
        automatically called by the main algorithm and can be overridden.'''
        # TODO: overload this function for a given space if necessary.
        return newPointInSpace

    def getCurrentNList(self):
        ''' Returns the nList as it currently stands.'''
        return self.nList

    def getNList(self):
        ''' Returns the longest nList found.'''
        return self.longestNList
    
    
    
if __name__=="__main__":
    NBBGenerator =  NPackBB("NPackExample.txt")
    pointsToAvoid=[ np.array([float(x + 0.5), float(y + 0.5), 2.0]) for x in range(-5,5) for y in range(-5,5)]
    BB = NBBGenerator.generateBuildingBlock(200, 1, visualiseEnvelope=(3000,50), envelopeList=['frustum 20 20 -3 15'], pointsToAvoid=pointsToAvoid, showBlockDirector=True)
    BB.exportBBK("NPackExample")
        
        