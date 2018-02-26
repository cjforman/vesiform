import numpy as np
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.VolumePackEllipsoid import VolumePackEllipsoidBBG as VPEBBG
import Utilities.fileIO as fIO
peptideBBG = PBG('betastrand.txt')
sphereBBG = SPSBBG('SurfacePackSphere.txt')
VolSphereBBG = VPEBBG('VolumePackEllipsoid.txt')

numResidues = 5
minDist = 1.0

numStrands = 100
centerPos = np.array([-0,-0, -0])
director = np.array([0, 0, 1])
rotation = 0
radius1 = 22
radius2 = 12
radius3 = 8
theta1 = -90
theta2 = 90
phi1 = -135
phi2 = 135
minDist1 = 3.0
minDist2 = 1.0
minDist3 = 1.0

# generate the XYZVals in the packed spaced
SphereBB = sphereBBG.generateBuildingBlock(numStrands, radius1, theta1, theta2, phi1, phi2, minDist1)
SphereBB.transformBBToLabFrame(director, centerPos, rotation)
SphereBB.exportBBK("sphereBasePoints")
SpherePoints = SphereBB.blockXYZVals


SphereBB1 = sphereBBG.generateBuildingBlock(int(2.5*numStrands), radius2, theta1, theta2, phi1, phi2, minDist2)
SphereBB1.transformBBToLabFrame(director, centerPos, rotation)
SphereBB1.exportBBK("sphereBasePoints1")

numSpherePoints = 200
volPointBB = VolSphereBBG.generateBuildingBlock(numSpherePoints, 
                                                radius3, 
                                                radius3,
                                                radius3, 
                                                minDist3)
volPointBB.exportBBK("VolSpherePoints.xyz")


peptideStrands = [ peptideBBG.generateBuildingBlock(numResidues) for _ in range(numStrands) ]
directors = [ (pos - centerPos)/np.linalg.norm(pos - centerPos) for pos in SpherePoints] 

vesicleBBS = [ strand.transformBBToLabFrame(director, pos, 0.0) for director, pos, strand in zip(directors, SpherePoints, peptideStrands)]

xyzVals = peptideStrands[0].blockXYZVals
names =  peptideStrands[0].blockAtomNames
for strand in peptideStrands[1:]:
    xyzVals = np.concatenate( (xyzVals, strand.blockXYZVals), 0)
    names = np.concatenate( (names, strand.blockAtomNames), 0)

fIO.saveXYZList(xyzVals, names, "peptideVesicle.xyz")


print "example done"