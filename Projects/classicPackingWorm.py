import numpy as np
from Library.SurfacePackCylinder import SurfacePackCylinderBBG as SPCBBG
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO

CylinderBBG = SPCBBG('SurfacePackCylinder.txt')

bondlength = 1.5

Polymer1Generator = RPBBG('RandomPolymer.txt')

numA = 30
numB = 40
numMonomersPerPolymer = numA + numB
numPolymersPerCylinder =  430 #430
CZ1 = 0
CZ2 = 150
FMaxRadius = 10
FMinRadius = 2
FZ1 = 80
FZ2 = 3
alpha1= 60
alpha2 = 80
beta1= 150
beta2 = 170
AtomicMinDist = 1.0

centerPos = np.array([0.0, 0.0, 0.0])

# generate the XYZVals packed in the outer cylinder
Polymer1CylinderBB = CylinderBBG.generateBuildingBlock(numPolymersPerCylinder, FZ2, FZ2, CZ1, CZ2, -180, 180, FMinRadius)
Polymer1CylinderBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
Polymer1CylinderBB.exportBBK("cylinderBasePoints")
Polymer1CylinderPoints = Polymer1CylinderBB.blockXYZVals 

# take the xyzpoint of each cylinder point and convert it to cylindrical coords
Polymer1CylinderPointsCylCoords = [ coords.XYZ2Cyl(pos) for pos in Polymer1CylinderPoints ] 

# find the radial director vector at each point
directors = [ coords.cyl2XYZ(np.array([pos[0], pos[1], 0.0])) for pos in Polymer1CylinderPointsCylCoords ]
directorsHat = [ director/np.linalg.norm(director) for director in directors]

# compute the base point of each vector on the inner cylinder
basePoints = [ coords.cyl2XYZ(np.array([FZ2, pos[1], pos[2]])) for pos in Polymer1CylinderPointsCylCoords ]

# set up the frustum to start on the inner surface
envelopeList = ['frustum ' + str(FZ1) + ' ' + str(FMaxRadius) + ' ' + str(FZ2 - AtomicMinDist) + ' ' + str(FMinRadius)]
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
                                                visualiseEnvelope=(0, 100)) for _ in range(numPolymersPerCylinder) ]

[  strand.setBlockRefPoint(polymerStartPoint) for strand in polymer1Strands ]

names =  ['O'] * numA 
names = np.concatenate( (names, ['C'] * numB), 0 ) 

strandNum = 0
for strand in polymer1Strands:
    strand.blockAtomNames = names[:]
    strand.exportBBK("strand" + str(strandNum))
    strandNum += 1

AllNames = []

curStrand = 0
for directorHat, pos, strand in zip(directorsHat, basePoints, polymer1Strands):
    strand.transformBBToLabFrame(directorHat, pos, 0.0)
    if curStrand==0:
        xyzVals = strand.blockXYZVals
        allNames = strand.blockAtomNames
    else:
        xyzVals = np.concatenate((xyzVals, strand.blockXYZVals), 0)
        allNames = np.concatenate( (allNames, strand.blockAtomNames), 0)
    curStrand += 1

fIO.saveXYZList(xyzVals, allNames, "polymerCylinder.xyz")


print "example done"