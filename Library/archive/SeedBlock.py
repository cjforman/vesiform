'''
Created on 16 Jan 2018

@author: chris
'''
import sys
import numpy as np
import copy as cp
from Builder.BuildingBlock import BuildingBlock
import Library.Sphere as sphere
import Library.RandomPolymer as rPoly



class SeedBlock(BuildingBlock):
    # Distributes a list of building blocks around a list of points and orients each building block.
    # Creates a single seed building block which takes a list of seedpoints, a list of directors 
    # and a list of building blocks all the same length.
    # The reference point of each building block is placed at the corresponding seedPoint 
    # and the building block is rotated so that it's orientation vector is aligned with
    # the director associated with the seed point.
    # Each listBlock is assigned to the seedpoints in the same order as the lists
    # 
    # as the buildingblocks are created and oriented their connectors, xyzVals and names are 
    # copied into the child building block of the seedblock.
    #  
    # When the relevant access function is called a building block is created 
    # which uses the connectors, xyzVals, names, refPoint and orientation of the child BuildingBlock 

    def __init__(self, seedPoints, directors, BBList, orientation, refPoint):
        # seedpoints is an array of n np 3-arrays
        # direcotrs is an array of n np 3-arrays
        # BBList is a list of n buildingBlocks
        # orientation is an np array defining orientation of convolutedSeedBlock
        # refPoint is an np array defining ref Point of convolutedSeedBlock
        
        if len(seedPoints)==len(directors)==len(BBList):
            
            # call the Building Block constructor with dummy variables. to set up the object
            BuildingBlock.__init__(self, np.array([]), [], np.array([[]]), refPoint, orientation)

            # add the seedBlock variables            
            self.seedPoints = seedPoints
            self.directors = directors
            self.BBList = BBList
            self.convoluteSeedBlock()
        else:
            print "seedpoints, directors and building block list are unequal lengths"
            sys.exit()

        return self.getSeedBlock()

    def convoluteSeedBlock(self):
        
        for seedPoint, curBB, curDirector in zip(self.seedPoints, self.BBList, self.directors):
            
            # don't modify the original building block - just copy it.
            tempBB = cp.copy(curBB)
            
            # Align the building block so it's refPoint is now at seedPoint
            tempBB.translateRefPointToTarget(seedPoint)
            
            # Rotate the building block about it's refPoint so it's orientation is parallel with director 
            tempBB.orientToDirector(curDirector)
            
            # get the index base Point
            indexOffset = self.countAtoms() 
            
            # add the reorientated atoms to the seedPoint Building Block  
            self.addAtoms(tempBB.getAtoms())
            
            # add the curBB connectors to the seedBlock connectors, adding the indexOffset
            self.addConnectors(tempBB.getAllConnectionIndices(), indexOffset)
        
    def getSeedBlock(self):
        # returns a building block containing the convoluted seed block
        return BuildingBlock(self.xyzVals, self.atomNames, self.connectors, self.refPoint, self.orientation)
        

if __name__=="__main__":
        
        # generate the seed points
        sphereGenerator = sphere.SphereNBB("sphere.txt")
        maxNumAttempts = 10
        curNumAttempts = 0
        while (not sphereGenerator.generateNConstrainedRandomPoints() and curNumAttempts < maxNumAttempts):
            curNumAttempts += 1
        
        seedPointList = sphereGenerator.getNList()
        # the directors of a sphere are unit vectors pointing in the direction of each point.
        directors = [ seedPoint/np.linalg.norm(seedPoint) for seedPoint in seedPointList]
        
        # generator the polymers for placing at each point in the sphere
        polymerGenerator = rPoly.RandomPolymerPackNBB("unimer.txt")

        # count the number of points in the seedPoints
        numPoints = len(seedPointList)

        # set up array to store the building blocks
        bbList = []
        
        # loop to generate the building blocks for the unimers
        curBB = 0
        while (curBB < numPoints):
            print curBB, " out of ", numPoints
            # generate a new polymer
            if polymerGenerator.generateNConstrainedRandomPoints():
                polymerXYZList = polymerGenerator.getNList()
                unimerLen = polymerGenerator.numPoints
                bbList.append( BuildingBlock( polymerXYZList,
                                              ['C'] * unimerLen,
                                              np.array([[2, 1, 0], [unimerLen -3, unimerLen - 2, unimerLen - 1]]),
                                              polymerXYZList[0], 
                                              polymerXYZList[-1] - polymerXYZList[0]) )
                curBB += 1
            else:
                print "re-doing that one" 
                curBB -= 1 
        
        sBlock = SeedBlock(seedPointList, directors, bbList, np.array([0,0,1]), np.array([0,0,1]))
        
        # perform the convolution of each unimer
        sBlock.convoluteSeedBlock()
        
        # get the overall seedBlock now just a regular BuildingBlock
        sBlockAll = sBlock.getSeedBlock()
        
        sBlockAll.exportBBK("vesicle")
        
        print "done"
        