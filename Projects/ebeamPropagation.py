import numpy as np
import copy as cp
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
from Library.XYZBuildingBlock import XYZBuildingBlockGenerator as XYZBBG
import Utilities.fileIO as fIO

SphereBBG = SPSBBG('SurfacePackSphere.txt')
Polymer1Generator = RPBBG('RandomPolymer.txt')
XYZGenerator = XYZBBG()

numPEG = 10
numHPMA = 0
numTriS = 11
PEGAlpha1= 40.0
PEGAlpha2 = 80.0
PEGBeta1= 130.0
PEGBeta2 = 165.0
HPMAAlpha1 = 20.0
HPMAAlpha2 = 40.0
HPMABeta1= 150.0
HPMABeta2 = 165.0
TriSAlpha1 = 40.0
TriSAlpha2 = 80.0
TriSBeta1= 130.0
TriSBeta2 = 165.0

numPolymersPerSphere = 340
FactorA = 2.5
FactorB = 1.0
AtomicMinDist = 1.0
bondLength = 1.5
FMaxRadius = 8.0
FMinRadius = 1.5
FZ1 = 80
FZ2 = 10
AtomicMinDistA = FactorA * AtomicMinDist
bondLengthA = FactorA * bondLength
FMaxRadiusA = FactorA * FMaxRadius
FMinRadiusA = FactorA * FMinRadius
FZ1A = FactorA * FZ1
FZ2A = FactorA * FZ2

AtomicMinDistB = FactorB * AtomicMinDist
bondLengthB = FactorB * bondLength
FMaxRadiusB = FactorB * FMaxRadius
FMinRadiusB = FactorB * FMinRadius
FZ1B = FactorB * FZ1
FZ2B = FactorB * FZ2




centerPos = np.array([0.0, 0.0, 0.0])

# generate the XYZVals in the packed spaced
SphereBB = SphereBBG.generateBuildingBlock(numPolymersPerSphere, FZ2A, -90, 90, -180, 180, FMinRadiusA)
SphereBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
SphereBB.exportBBK("sphereBasePoints")
SpherePoints = SphereBB.blockXYZVals 

envelopeListA = ['frustum ' + str(FZ1A) + ' ' + str(FMaxRadiusA) + ' ' +str(FZ2A - AtomicMinDistA) + ' ' +str(FMinRadiusA)]
envelopeListB = ['frustum ' + str(FZ1B) + ' ' + str(FMaxRadiusB) + ' ' +str(FZ2B - AtomicMinDistB) + ' ' +str(FMinRadiusB)]

polymerStartPointA = np.array([0.0, 0.0, FZ2A + AtomicMinDistA])
polymerStartPointB = np.array([0.0, 0.0, FZ2B + AtomicMinDistB])

PegStrands = [ Polymer1Generator.generateBuildingBlock(numPEG, 
                                                       polymerStartPointA,
                                                       PEGAlpha1,
                                                       PEGAlpha2,
                                                       PEGBeta1, 
                                                       PEGBeta2,
                                                       AtomicMinDistA,
                                                       bondLengthA,
                                                       envelopeList=envelopeListA,
                                                       visualiseEnvelope=(0, 100)) for _ in range(numPolymersPerSphere) ]

[  strand.setBlockRefPoint(polymerStartPointA) for strand in PegStrands ]

if numHPMA>0:
    HPMAStrands = [ Polymer1Generator.generateBuildingBlock(   numHPMA, 
                                                               polymerStartPointA,
                                                               HPMAAlpha1,
                                                               HPMAAlpha2,
                                                               HPMABeta1, 
                                                               HPMABeta2,
                                                               AtomicMinDistA,
                                                               bondLengthA,
                                                               envelopeList=envelopeListA,
                                                               visualiseEnvelope=(0, 100)) for _ in range(numPolymersPerSphere) ]
    
    [  strand.setBlockRefPoint(polymerStartPointA) for strand in HPMAStrands ]

TriSStrands = [ Polymer1Generator.generateBuildingBlock(   numTriS, 
                                                           polymerStartPointB,
                                                           TriSAlpha1,
                                                           TriSAlpha2,
                                                           TriSBeta1, 
                                                           TriSBeta2,
                                                           AtomicMinDistB,
                                                           bondLengthB,
                                                           envelopeList=envelopeListB,
                                                           visualiseEnvelope=(0, 100)) for _ in range(numPolymersPerSphere) ]

[  strand.setBlockRefPoint(polymerStartPointB) for strand in TriSStrands ]

sTNB = np.identity(3)
S1 = np.array([0.0, 0.0, 0.0])
C1 = np.array([bondLengthB, 0.0, 0.0])
S2 = np.array([bondLengthB * (1.0 + np.cos( (180.0 - 110.0) * np.pi/180.0 )), bondLengthB * np.sin( (180.0 - 110.0) * np.pi/180.0 ), 0.0])
S3 = np.array([bondLengthB * (1.0 + np.cos( (110.0/2.0) * np.pi/180.0 )), -1.0 * bondLengthB * np.sin( ( 110.0 / 2.0) * np.pi/180.0), 0.0])
fIO.saveXYZList([S1, C1, S2, S3], ['S', 'C', 'S', 'S'], 'S3C.xyz')
S3C = XYZGenerator.generateBuildingBlock('S3C.xyz')

names =  ['N'] * numPEG 
if numHPMA>0:
    names = np.concatenate( (names, ['O'] * numHPMA), 0 ) 
names = np.concatenate( (names, ['S', 'C', 'S', 'S']), 0 )
names = np.concatenate( (names, ['C'] * numTriS), 0 ) 
  
strandNum = 0
if numHPMA>0:
    for PEG, HPMA, TriS in zip(PegStrands, HPMAStrands, TriSStrands):
        strand, staple1 = PEG.addBuildingBlock(HPMA, 1, 0, bondLengthA, -57.0, 116.0, 180.0, 122.0, -47.0, polymerStartPointA, np.array([0.0, 0.0, 1.0]))
        strand, staple2 = strand.addBuildingBlock(S3C, 1, 0, 0.8 * bondLengthA, -57.0, 116.0, 47.0, 160.0, -47.0, polymerStartPointA, np.array([0.0, 0.0, 1.0]))
        strand, staple3 = strand.addBuildingBlock(TriS, 1, 0, bondLengthB, 157.0, 160.0, 47.0, 122.0, -47.0, polymerStartPointA, np.array([0.0, 0.0, 1.0]))
        strand.blockAtomNames = names[:]
        strand.exportBBK("strand_" + str(numHPMA) + "_" + str(strandNum))
        if strandNum==0:
            strandList = [cp.deepcopy(strand)]
        else:
            strandList.append(cp.deepcopy(strand))
        strandNum += 1
else:
    for PEG, TriS in zip(PegStrands, TriSStrands):
        strand, staple1 = PEG.addBuildingBlock(S3C, 1, 0, 0.8 * bondLengthA, -57.0, 116.0, 47.0, 160.0, -47.0, polymerStartPointA, np.array([0.0, 0.0, 1.0]))
        strand, staple3 = strand.addBuildingBlock(TriS, 1, 0, bondLengthB, 157.0, 160.0, 47.0, 122.0, -47.0, polymerStartPointA, np.array([0.0, 0.0, 1.0]))
        strand.blockAtomNames = names[:]
        strand.exportBBK("strand_" + str(numHPMA) + "_" + str(strandNum))
        if strandNum==0:
            strandList = [cp.deepcopy(strand)]
        else:
            strandList.append(cp.deepcopy(strand))
        strandNum += 1

directors = [ (pos - centerPos)/np.linalg.norm(pos - centerPos) for pos in SpherePoints] 

AllNames = []

curStrand = 0
for director, pos, strand in zip(directors, SpherePoints, strandList):
    strand.transformBBToLabFrame(director, pos, 0.0)
    if curStrand==0:
        xyzVals = strand.blockXYZVals
        allNames = strand.blockAtomNames
    else:
        xyzVals = np.concatenate((xyzVals, strand.blockXYZVals), 0)
        allNames = np.concatenate( (allNames, strand.blockAtomNames), 0)
    curStrand += 1

fIO.saveXYZList(xyzVals, allNames, "sphere_" + str(numHPMA) + ".xyz")

print "ebeam done"