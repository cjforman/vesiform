import numpy as np
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO

SphereBBG = SPSBBG('SurfacePackSphere.txt')
PolymerGenerator = RPBBG('RandomPolymer.txt')

numA = 30
numB = 40
numMonomersPerPolymer = numA + numB
bondLength = 1.5
maxPolymerLength = float(numMonomersPerPolymer) * bondLength
minRadius = 50
midLayerRadius = minRadius + maxPolymerLength
minTheta = -90
maxTheta = 90
minPhi = -135
maxPhi = 135

numPolymersPerOuterSphere = 1100
outerSphereFrustumOuterZ = midLayerRadius + bondLength/2.0 + maxPolymerLength  
outerSphereFrustumInnerZ = midLayerRadius + bondLength/2.0
outerSphereFrustumOuterRadius = np.sqrt(0.75 * 4.0  * np.power(outerSphereFrustumOuterZ, 2.0) / numPolymersPerOuterSphere)  
outerSphereFrustumInnerRadius = np.sqrt(0.75 * 4.0  * np.power(outerSphereFrustumInnerZ, 2.0) / numPolymersPerOuterSphere)
outerPolymerDihedralMin = 50
outerPolymerDihedralMax = 90
outerPolymerBondAngleMin = 100
outerPolymerBondAngleMax = 170

numPolymersPerInnerSphere = 500
innerSphereFrustumOuterZ = midLayerRadius - bondLength/2.0 
innerSphereFrustumInnerZ = midLayerRadius - bondLength/2.0 - maxPolymerLength
innerSphereFrustumOuterRadius = np.sqrt(0.75 * 4.0  * np.power(innerSphereFrustumOuterZ, 2.0) / numPolymersPerInnerSphere)
innerSphereFrustumInnerRadius = np.sqrt(0.75 * 4.0  * np.power(innerSphereFrustumInnerZ, 2.0) / numPolymersPerInnerSphere)
innerPolymerDihedralMin = 60
innerPolymerDihedralMax = 80
innerPolymerBondAngleMin = 130
innerPolymerBondAngleMax = 170

print "Outer: numPolymers:", numPolymersPerOuterSphere, "outerZ:", outerSphereFrustumOuterZ, "outerR", outerSphereFrustumOuterRadius, "innerZ:", outerSphereFrustumInnerZ, "innerR", outerSphereFrustumInnerRadius    
print "Inner: numPolymers:", numPolymersPerInnerSphere, "outerZ:", innerSphereFrustumOuterZ, "outerR", innerSphereFrustumOuterRadius, "innerZ:", innerSphereFrustumInnerZ, "innerR", innerSphereFrustumInnerRadius

atomicMinDist = 1.0

centerPos = np.array([0.0, 0.0, 0.0])

# generate the outer sphere root positions at the mid layer radius and their directors pointing out from the center 
PolymerSphereOuterBB = SphereBBG.generateBuildingBlock(numPolymersPerOuterSphere, outerSphereFrustumInnerZ, minTheta, maxTheta, minPhi, maxPhi, outerSphereFrustumInnerRadius)
PolymerSphereOuterBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
PolymerSphereOuterBB.exportBBK("OuterSphereBasePoints")
PolymerSphereOuterPoints = PolymerSphereOuterBB.blockXYZVals 
PolymerSphereOuterPointsPolar = [ coords.XYZ2SphericalPolar(pos) for pos in PolymerSphereOuterBB.blockXYZVals]
outerDirectorsHat = [ coords.sphericalPolar2XYZ(np.array([1.0, pos[1], pos[2]]) ) for pos in PolymerSphereOuterPointsPolar]

# generate the inner sphere root positions at the mid layer radius and their directors pointing in towards the center
PolymerSphereInnerBB = SphereBBG.generateBuildingBlock(numPolymersPerInnerSphere, innerSphereFrustumOuterZ, minTheta, maxTheta, minPhi, maxPhi, innerSphereFrustumOuterRadius)
PolymerSphereInnerBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
PolymerSphereInnerBB.exportBBK("OuterSphereBasePoints")
PolymerSphereInnerPoints = PolymerSphereInnerBB.blockXYZVals 
PolymerSphereInnerPointsPolar = [ coords.XYZ2SphericalPolar(pos) for pos in PolymerSphereInnerBB.blockXYZVals]
innerDirectorsHat = [ -1.0 * coords.sphericalPolar2XYZ(np.array([1.0, pos[1], pos[2]]) ) for pos in PolymerSphereInnerPointsPolar]

print "Outer: numPolymers:", numPolymersPerOuterSphere, "outerZ:", outerSphereFrustumOuterZ, "outerR", outerSphereFrustumOuterRadius, "innerZ:", outerSphereFrustumInnerZ, "innerR", outerSphereFrustumInnerRadius    
print "Inner: numPolymers:", numPolymersPerInnerSphere, "outerZ:", innerSphereFrustumOuterZ, "outerR", innerSphereFrustumOuterRadius, "innerZ:", innerSphereFrustumInnerZ, "innerR", innerSphereFrustumInnerRadius


# generate the outer polymers
outerEnvelopeList = ['frustum ' + str(outerSphereFrustumOuterZ) + ' ' + str(outerSphereFrustumOuterRadius) + ' ' +str(outerSphereFrustumInnerZ - atomicMinDist) + ' ' +str(outerSphereFrustumInnerRadius)]
outerPolymerStartPoint = np.array([0.0, 0.0, outerSphereFrustumInnerZ])
outerPolymerStrands = []
names =  ['O'] * numA # red Hydrophobic
names = np.concatenate( (names, ['C'] * numB), 0 ) # blue hydrophilic (later monomers blue are towards outer shell) 
for strandNum in range(numPolymersPerOuterSphere): 
    print "starting strand: ", strandNum
    envelopeSize = 0
    if strandNum == 0:
        envelopeSize = 1000000
    strand = PolymerGenerator.generateBuildingBlock( numMonomersPerPolymer, 
                                                     outerPolymerStartPoint,
                                                     outerPolymerDihedralMin,
                                                     outerPolymerDihedralMax,
                                                     outerPolymerBondAngleMin, 
                                                     outerPolymerBondAngleMin,
                                                     atomicMinDist,
                                                     bondLength,
                                                     envelopeList=outerEnvelopeList,
                                                     visualiseEnvelope=(envelopeSize, 400, 'outerEnvelope.xyz'))
    strand.setBlockRefPoint(outerPolymerStartPoint)
    strand.blockAtomNames = names[:]
    strand.exportBBK("outer_strand_" + str(strandNum))
    outerPolymerStrands.append(strand)

# generate the inner polymers
innerSphereFrustumInnerZInverted = innerSphereFrustumOuterZ - innerSphereFrustumInnerZ + innerSphereFrustumOuterZ 
innerEnvelopeList = ['frustum ' + str(innerSphereFrustumOuterZ - atomicMinDist) + ' ' + str(innerSphereFrustumOuterRadius) + ' ' +str(innerSphereFrustumInnerZInverted) + ' ' +str(innerSphereFrustumInnerRadius)]
innerPolymerStartPoint = np.array([0.0, 0.0, innerSphereFrustumOuterZ])
innerPolymerStrands = []
for strandNum in range(numPolymersPerInnerSphere): 
    print "starting strand: ", strandNum
    envelopeSize = 0
    if strandNum == 0:
        envelopeSize = 1000000
    strand = PolymerGenerator.generateBuildingBlock( numMonomersPerPolymer, 
                                                     innerPolymerStartPoint,
                                                     innerPolymerDihedralMin,
                                                     innerPolymerDihedralMax,
                                                     innerPolymerBondAngleMin, 
                                                     innerPolymerBondAngleMin,
                                                     atomicMinDist,
                                                     bondLength,
                                                     envelopeList=innerEnvelopeList,
                                                     visualiseEnvelope=(envelopeSize, 400, 'innerEnvelope.xyz'))
    strand.setBlockRefPoint(innerPolymerStartPoint)
    strand.blockAtomNames = names[:]
    strand.exportBBK("inner_strand_" + str(strandNum))
    innerPolymerStrands.append(strand)

AllNames = []

curStrand = 0
for director, pos, strand in zip(outerDirectorsHat, PolymerSphereOuterPoints, outerPolymerStrands):
    strand.transformBBToLabFrame(director, pos, 0.0)
    if curStrand==0:
        xyzVals = strand.blockXYZVals
        allNames = strand.blockAtomNames
    else:
        xyzVals = np.concatenate((xyzVals, strand.blockXYZVals), 0)
        allNames = np.concatenate( (allNames, strand.blockAtomNames), 0)
    curStrand += 1

for director, pos, strand in zip(innerDirectorsHat, PolymerSphereInnerPoints, innerPolymerStrands):
    strand.transformBBToLabFrame(director, pos, 0.0)
    xyzVals = np.concatenate((xyzVals, strand.blockXYZVals), 0)
    allNames = np.concatenate( (allNames, strand.blockAtomNames), 0)
    curStrand += 1

fIO.saveXYZList(xyzVals, allNames, "Vesicle.xyz")

print "classic vesicle example done"