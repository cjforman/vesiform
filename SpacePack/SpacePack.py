'''
Created on 14 Dec 2017

@author: chris
'''
from Builder.BuildingBlockGenerator import BuildingBlockGenerator as BBG 
import sys
import numpy as np
import random as rnd
from Utilities import fileIO as fIO

class SpacePackBBG(BBG):
    '''
    A function that exploits the envelope framework to pack an arbitrary shaped space
    with a uniform distribution of objects. 
    '''
    def __init__(self, filename):
        BBG.__init__(self, filename)
        
    def initialiseParameters(self):
        # ensure parent initialisation takes place and core values are initialised 
        BBG.initialiseParameters(self) 
        
        self.maxNumPackMoves  = self.getParam("maxNumPackMoves")
        self.packT = self.getParam("packT")
        self.verbose = self.getParam("verbose")
        
        if self.noLoadErrors == False:            
            print "Critical Parameters are undefined for SpacePack object"
            sys.exit()        
    
    def generateBuildingBlock(self, numPoints, XRange, YRange, ZRange, minDist, visualiseEnvelope=(0,20), envelopeList=['None'], pointsToAvoid=[]):
        # XRange, YRange and ZRange are a pair of min and max points.
        
        # Note the number of points
        self.numPoints = numPoints
        self.minDist = minDist
        
        # specify the min and max of X
        self.XRange = XRange
        self.YRange = YRange
        self.ZRange = ZRange
        
        self.pointsToAvoid = pointsToAvoid
        
        return BBG.generateBuildingBlock(self, numPoints, minDist, visualiseEnvelope=visualiseEnvelope, envelopeList=envelopeList, pointsToAvoid = pointsToAvoid, )
    
    def generateBuildingBlockDirector(self):
        director = np.array([0.0, 0.0, 1.0])
        return director

    def generateBuildingBlockRefPoint(self):
        return np.array([0.0, 0.0, 0.0])
    
    def generateBuildingBlockXYZ(self):
        
        # generate a set of N random points in the outer cubic region.
        # At this stage don't care about anything except that there are the 
        # correct number of points in this region.
        curXYZVals = self.generateNRandPoints(self.numPoints)

        # make a copy of curXYZ Vals which is the best full set yet,
        bestXYZVals = curXYZVals[:]
        
        print "Checking Initial list for out of bounds points."
        # make a list of which points indices are in violation of the conditions 
        curPointsOutOfBounds = self.checkNewPointsAllTests( range(len(curXYZVals)), curXYZVals, curXYZVals, self.pointsToAvoid) 
            
        # make a note of the current number of points out of bounds.
        bestNumOutOfBoundPoints = len(curPointsOutOfBounds)
        curNumOutOfBoundPoints = len(curPointsOutOfBounds)

        # start a move counter        
        numMoves = 0

        # loop until there are no more out of bounds points or we have taken the max number of moves allowed by the user.        
        while bestNumOutOfBoundPoints > 0 and numMoves < self.maxNumPackMoves:
            
            # for each point that is Out Of Bounds generate a new point
            replacementVals = self.generateNRandPoints(curNumOutOfBoundPoints)

            # create a temporary list with the replaced points
            newXYZVals = curXYZVals[:]
            for index, replacementVal in zip(curPointsOutOfBounds, replacementVals):
                newXYZVals[index] = replacementVal
            
            # test the new points to see how they fare against the rest of the list 
            newPointsOutOfBounds = self.checkNewPointsAllTests(curPointsOutOfBounds, replacementVals, newXYZVals, self.pointsToAvoid) 
            
            # count the number of points in the new set that are out of bounds
            newNumOutOfBoundPoints = len(newPointsOutOfBounds)
            
            # if the new set of points is better than the current best set of points keep it as the new best set.
            if newNumOutOfBoundPoints < bestNumOutOfBoundPoints:
                bestXYZVals=newXYZVals[:]
                bestNumOutOfBoundPoints = newNumOutOfBoundPoints
                print "numMoves: ", numMoves, " bestNumOutOfBoundPoints: ", bestNumOutOfBoundPoints, "curNumOutOfBoundPoints:", curNumOutOfBoundPoints
                
            # if the new number of bad points is lower than the number associated with the current working copy 
            # then update the current working copy to the new copy
            if newNumOutOfBoundPoints < curNumOutOfBoundPoints:
                curXYZVals = newXYZVals[:]
                curNumOutOfBoundPoints = newNumOutOfBoundPoints
                curPointsOutOfBounds = newPointsOutOfBounds[:]
                # update the current copy
            else:
                # The new number of points outside the zone is higher than the current set.
                # We roll the dice. If the dice is lower than the threshold we revert back to 
                # the previous working set. Otherwise we keep this one.
                # The probability of keeping the new set as the current set shrinks exponentially 
                # with the difference on the number out of bounds vs the current best set
                if rnd.uniform(0.0, 1.0) < np.exp( (bestNumOutOfBoundPoints - newNumOutOfBoundPoints)/(self.packT * self.numPoints) ):
                    curXYZVals = newXYZVals[:]
                    curNumOutOfBoundPoints = newNumOutOfBoundPoints
                    curPointsOutOfBounds = newPointsOutOfBounds[:]                    
            
            numMoves += 1
            
            if numMoves % 20 == 0:
                print "step: ", numMoves, " bestNumOutOfBoundPoints: ", bestNumOutOfBoundPoints, "curNumOutOfBoundPoints:", curNumOutOfBoundPoints
             
        if bestNumOutOfBoundPoints > 0:
            print "Warning: Packing unable to pack all points into defined region." 
    
    
        print "step: ", numMoves, " bestNumOutOfBoundPoints: ", bestNumOutOfBoundPoints, "curNumOutOfBoundPoints:", curNumOutOfBoundPoints
        return bestXYZVals 
    
    def pickAPoint(self):
        return np.array([rnd.uniform(self.XRange[0], self.XRange[1]),
                         rnd.uniform(self.YRange[0], self.YRange[1]),
                         rnd.uniform(self.ZRange[0], self.ZRange[1])])
    
    def generateNRandPoints(self, N):
        return [ self.pickAPoint() for _ in range(N)]
    
    def checkPointAgainstDynamicList(self, index, posToTest, dynamicList):
        # assume point is good
        newPointGood = True 
        
        # set up array to export points made bad by proximity to this point
        badPoints = []

        for n, pos in enumerate(dynamicList):
            if not n == index:
                if np.linalg.norm(posToTest - pos) < self.minDist:
                    newPointGood=False
                    badPoints.append(n)
        
        # a flag saying whether or not it passed the test and a list of points it intersects with
        return newPointGood, badPoints

    def checkPointAgainstStaticList(self, posToTest, staticList):        
        # assume success
        goodPoint = True
        for pos in staticList:
            if np.linalg.norm(posToTest - pos) < self.minDist:
                goodPoint = False
                if self.verbose==1:
                    print "Points To Avoid Violation"
                break
        return goodPoint
        
    def checkNewPointsAllTests(self, pointIndices, newPoints, allPoints, pointsToAvoid):
        ''' Creates a list of all the points in XYZ that failed the tests.'''
        # set up output arrays        
        badPoints = []
        extraBadPoints = []
        
        for pointIndex, pointToTest in zip(pointIndices, newPoints):
            
            # perform spatial test
            pointGood = self.checkPointInBounds(pointToTest)

            # only perform static test if necessary            
            if pointGood:
                pointGood = self.checkPointAgainstStaticList(pointToTest, pointsToAvoid)
            
            # only perform dynamic test if necessary
            if pointGood:
                # test the point against all the other points except itself
                pointGood, extraBadPoints = self.checkPointAgainstDynamicList(pointIndex, pointToTest, allPoints)
            
            # if the point is bad then add the index to the array of bad points        
            if not pointGood:
                badPoints.append(pointIndex)
                
            # add any extraBadPoints to the badPoints array and eliminate any repeated indices
            badPoints = list(set(badPoints).union(set(extraBadPoints))) 
            
            # clear up the extra bad points list for the next go around
            extraBadPoints = []
        
        # return the index array of bad points and the positions themselves 
        return badPoints
         
    def getParams(self):
        return self.params 
        
if __name__ == '__main__':
    
    # get the file name from the command line
    filename = sys.argv[1]

    # create the NPack object.
    SpacePackBBG = SpacePackBBG(filename)

    numPoints = 1000
    centrePos = np.array([0.0, 0.0, 0.0])
    director= np.array([0.0, 0.0, 1.0])
    rotation = 0 
    xR = [-20, 20]
    yR = [-20, 20]
    zR = [ 0, 30]
    minDist = .2
    
    envelopeList = ['endsphere 20.0 10.0']
    envelopeList.append('frustum 20.0 10.0 -10.0 0.0')
       
    # generate the building block
    SpacePackBB = SpacePackBBG.generateBuildingBlock(numPoints, xR, yR, zR, minDist, envelopeList = envelopeList, pointsToAvoid = [])
    SpacePackBB.transformBBToLabFrame(director, centrePos, rotation)
    # dump the list of values to file
    SpacePackBB.exportBBK(fIO.fileRootFromInfile(filename,'txt'))
    print "Space Pack Done"
    