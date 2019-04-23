'''
Created on 14 Dec 2017

@author: chris
'''
from NPack.NPackBB import NPackBB as NBB 
import sys
import random as rnd
import numpy as np
import itertools as it
from Utilities import fileIO as fIO


class PartialVolumePackCuboidBBG(NBB):
    '''
    Over loads the NPackBB object to give a non-uniform packing throughout a 3Space. 
    No point is less than minDist away from any other in the 
    cubic region of space, but the probability of finding a particle is higher nearer
    existing particles. The centre of the space is the origin. 
    '''
    def __init__(self, filename):
        NBB.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        NBB.initialiseParameters(self) 
        self.proximityScore = 0
        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for PlanePack object")
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, xRange1, xRange2, yRange1, yRange2, zRange1, zRange2, minDist, r0, proximityT, seedLength, envelopeList=['None'], pointsToAvoid=[], visualiseEnvelope=(0,200), showBlockDirector=False):
        self.xRange1 = xRange1 
        self.xRange2 = xRange2
        self.yRange1 = yRange1
        self.yRange2 = yRange2
        self.zRange1 = zRange1
        self.zRange2 = zRange2
        self.r0 = r0
        self.proximityScore = 0
        self.lenListForProximityScore = 0
        self.ProximityT = proximityT
        self.seedLength = seedLength
        return NBB.generateBuildingBlock(self, numPoints, minDist, envelopeList=envelopeList, pointsToAvoid=pointsToAvoid, visualiseEnvelope=visualiseEnvelope, showBlockDirector=False)
    
    def generateBuildingBlockDirector(self):
        return np.array([0,0,1])

    def generateBuildingBlockRefPoint(self):
        return np.array([(self.xRange1 + self.xRange2)/2, 
                         (self.yRange1 + self.yRange2)/2,
                         (self.zRange1 + self.zRange2)/2])
    
    def pickFirstPoints(self):
        from Library.VolumePackCuboid import VolumePackCuboidBBG as VPCBBG        
        
        BBGen = VPCBBG(self.paramFilename)
        # generate seedLength number of randomly disributed particles in the void.
        BB = BBGen.generateBuildingBlock(self.seedLength, self.xRange1, self.xRange2, self.yRange1, self.yRange2, self.zRange1, self.zRange2, self.minDist)
        self.nList = BB.blockXYZVals
        #self.nList = [np.array([rnd.uniform(self.xRange1, self.xRange2), rnd.uniform(self.yRange1, self.yRange2), rnd.uniform(self.zRange1, self.zRange2)])]
        self.nAttempts = self.seedLength * [0]
        self.lenListForProximityScore, self.proximityScore = self.computeEntireScore(self.nList, [self.r0])
        return True
        
    def computeEntireScore(self, particle_list, paramList):
        # compute the proximity score of the list
        score = 0
        
        pairs = it.combinations(particle_list, 2)
        numPairs = 0
        for pair in pairs:
            score += self.pairwiseScore(pair, paramList)
            numPairs += 1
            
        return len(particle_list), score
        
    def computePartialScore(self, pos, particle_list, paramList):
        # compute the contribution of the test particle to the proximity score against the rest of the list
        score = 0
        
        for lpos in particle_list:
            score += self.pairwiseScore((pos, lpos), paramList)
            
        return score
        
    def pairwiseScore(self, pair, paramList):
        # the bigger R is the larger the score. We aim to minimise score as we add particles. 
        r = np.linalg.norm(pair[1] - pair[0])/paramList[0]
        return r*r*r*r # score is essentially total distance between all the particles   
        
    def pickRandomPointInDefinedSpace(self):
        return np.array( [ rnd.uniform( self.xRange1, self.xRange2 ), 
                           rnd.uniform( self.yRange1, self.yRange2 ),
                           rnd.uniform( self.zRange1, self.zRange2 )] ) 

    def getParams(self):
        return self.params 
    
    def checkPointAgainstList(self, pos):
        ''' This overloaded function accepts new particles that move the system towards
        the ideal distribution and rejects, with a certain probability, particles that move the system away from the
        specified distribution.  Can literally use an energy score to compute the proximity of the particles with a stickiness 
        and range parameter. '''
        # assume position is good
        goodPos = True

        for zPos in self.nList:
            # if new position occurs within minDist of any of the existing positions in the list 
            # then throw a wobbly.
            if np.linalg.norm((zPos - pos)) <= self.minDist:
                goodPos = False
                if self.verbose==1: 
                    print("Packing violation - too close to another particle")
                    
        if goodPos:

            # compute the contribution of the current particle to the current score 
            partial_score = self.computePartialScore(pos, self.nList, [self.r0])

            # check that 2 or more particles were used to compute the current score otherwise accept the point.           
            if self.lenListForProximityScore > 1:
                
                # if the length of the list has changed since the last time a score was computed then recompute entire score
                if self.lenListForProximityScore != len(self.nList):
                    self.lenListForProximityScore, self.proximityScore = self.computeEntireScore(self.nList, [self.r0]) 
                
                # compute the previous score per pair - essentially the mean pairwise distance between the particles.  
                oldScorePerPair = self.proximityScore/ ( self.lenListForProximityScore * (self.lenListForProximityScore - 1) /2 ) 
                
                # compute the new score as it would be if we accept the test particle.
                newScorePerParticle = (self.proximityScore + partial_score)/ ((self.lenListForProximityScore + 1)  * self.lenListForProximityScore / 2)
                
                # compute the new means distance between the particles
                dMeanDistance = newScorePerParticle - oldScorePerPair # +ve if mean distance has grown
                
                prob = rnd.uniform(0,1) # pick a random number between 0 and 1
                
                # if dMeanDistance is +ve then new score is larger than old score 
                # the more +ve it is (the larger the increase in average distance) the smaller exp(-dR/T) becomes.
                # The randomly picked threshold between 0 and 1 corresponds to a mean distance.  
                # Reject the move if the actual increase is larger than the randomly picked mean distance.
                if np.exp( -dMeanDistance / self.ProximityT ) < prob :
                    goodPos = False
                    if self.verbose==1: 
                        print("New position increased mean distance between particles by too much.")                
                
                # if the dR is -ve (dist has gone down) then exp(-dR) is greater than 1 so we do nothing - leave goodPos as true and accept point.

        # if the pos was not rejected then remember the new score
        if goodPos:
            self.proximityScore += partial_score
            self.lenListForProximityScore += 1  
            # the actual point will be added to the list 
            # as part of the main algorithm, if it's rejected for some other reason, then we will recompute the score on in its  
            # entirety on the next go. otherwise, we just modify the mean slowly as we add a point. Bit of efficiency. 

        return goodPos
    
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    PlanePackBBG = PartialVolumePackCuboidBBG(filename)

    numPoints = 300
    centrePos = np.array([0, 0, 0])
    director= np.array([0, 0, 1])
    rotation = 0 
    xRange1 = -100
    xRange2 = 100
    yRange1 = -100
    yRange2 = 100
    zRange1= -100
    zRange2= 100
    minDist = .1
    r0 = .01
    proximityT = 100
    seedLength = 10
       
    # generate the building block
    SpacePackBB = PlanePackBBG.generateBuildingBlock(numPoints, xRange1, xRange2, yRange1, yRange2, zRange1, zRange2, minDist, r0, proximityT, seedLength)
    SpacePackBB.transformBBToLabFrame(director, centrePos, rotation)
    SpacePackBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    print("Done")
    