import numpy as np
import copy as cp
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG

import Utilities.fileIO as fIO

SphereBBG = SPSBBG('SurfacePackSphere.txt')

bondlength = 1.5

Polymer1Generator = RPBBG('RandomPolymer.txt')

numA = 35
numB = 35
numMonomersPerPolymer = numA + numB
numPolymersPerSphere = 88
FMaxRadius = 15
FMinRadius = 3
FZ1 = 50
FZ2 = 10
alpha1= 50
alpha2 = 90
beta1= 80
beta2 = 170
AtomicMinDist = 1.0

centerPos = np.array([0.0, 0.0, 0.0])

# generate the XYZVals in the packed spaced
Polymer1SphereBB = SphereBBG.generateBuildingBlock(numPolymersPerSphere, FZ2, -90, 90, -180, 180, FMinRadius)
Polymer1SphereBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
Polymer1SphereBB.exportBBK("sphereBasePoints")
Polymer1SpherePoints = Polymer1SphereBB.blockXYZVals 

envelopeList = ['frustum ' + str(FZ1) + ' ' + str(FMaxRadius) + ' ' +str(FZ2-AtomicMinDist) + ' ' +str(FMinRadius)]

polymerStartPoint = np.array([0.0, 0.0, FZ2])
polymer1Strands = [ Polymer1Generator.generateBuildingBlock(numMonomersPerPolymer, 
                                                polymerStartPoint,
                                                alpha1,
                                                alpha2,
                                                beta1, 
                                                beta2,
                                                AtomicMinDist,
                                                bondlength,
                                                envelopeList=envelopeList,
                                                visualiseEnvelope=(0, 100)) for _ in range(numPolymersPerSphere) ]

[  strand.setBlockRefPoint(polymerStartPoint) for strand in polymer1Strands ]

names =  ['O'] * numA 
names = np.concatenate( (names, ['C'] * numB), 0 ) 

for strand in polymer1Strands:
    strand.blockAtomNames = names[:]

directors = [ (pos - centerPos)/np.linalg.norm(pos - centerPos) for pos in Polymer1SpherePoints] 

AllNames = []

curStrand = 0
for director, pos, strand in zip(directors, Polymer1SpherePoints, polymer1Strands):
    strand.transformBBToLabFrame(director, pos, 0.0)
    if curStrand==0:
        xyzVals = strand.blockXYZVals
        allNames = strand.blockAtomNames
    else:
        xyzVals = np.concatenate((xyzVals, strand.blockXYZVals), 0)
        allNames = np.concatenate( (allNames, strand.blockAtomNames), 0)
    curStrand += 1

fIO.saveXYZList(xyzVals, allNames, "polymer1.xyz")


print "example done"