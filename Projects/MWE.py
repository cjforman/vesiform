import numpy as np
from Builder.BuildingBlockGenerator import BuildingBlockGenerator


BBG = BuildingBlockGenerator('MWEParameters.txt')
numAtoms = 30
startPosition = np.array([0.0, 0.0, 0.0])
orientationVector = np.array([0.0, 0.0, 1.0])
rotation = 45 * np.pi/180
minDist = 1.0
testBuildBlock = BBG.generateBuildingBlock(numAtoms, minDist)
testBuildBlock.transformBBToLabFrame(orientationVector, startPosition, rotation) 
testBuildBlock .exportBBK('example')
print "example done"