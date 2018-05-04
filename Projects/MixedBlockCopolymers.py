import numpy as np
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
import Utilities.fileIO as fIO


def makeSphereAndPolymers(numA, numB, numPolymersPerSphere, FMaxRadius, FMinRadius, FZ1, FZ2, alpha1, alpha2, beta1, beta2 , AtomicMinDist, bondLength, filename):
    SphereBBG = SPSBBG('SurfacePackSphere.txt')
    Polymer1Generator = RPBBG('RandomPolymer.txt')

    centerPos = np.array([0.0, 0.0, 0.0])

    # generate the XYZVals in the packed spaced
    Polymer1SphereBB = SphereBBG.generateBuildingBlock(numPolymersPerSphere, FZ2, -90, 90, -180, 180, FMinRadius)
    Polymer1SphereBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
    Polymer1SphereBB.exportBBK(filename + "_sphereBasePoints")
    Polymer1SpherePoints = Polymer1SphereBB.blockXYZVals 

    envelopeList = ['frustum ' + str(FZ1) + ' ' + str(FMaxRadius) + ' ' +str(FZ2 - AtomicMinDist) + ' ' +str(FMinRadius)]

    polymerStartPoint = np.array([0.0, 0.0, FZ2])
    polymer1Strands = [ Polymer1Generator.generateBuildingBlock(numMonomersPerPolymer, 
                                                                polymerStartPoint,
                                                                alpha1,
                                                                alpha2,
                                                                beta1, 
                                                                beta2,
                                                                AtomicMinDist,
                                                                bondLength,
                                                                envelopeList=envelopeList,
                                                                visualiseEnvelope=(0, 100)) for _ in range(numPolymersPerSphere) ]

    [  strand.setBlockRefPoint(polymerStartPoint) for strand in polymer1Strands ]

    names =  ['O'] * numA 
    names = np.concatenate( (names, ['C'] * numB), 0 ) 
    
    strandNum = 0
    for strand in polymer1Strands:
        strand.blockAtomNames = names[:]
        strand.exportBBK(filename + "_strand" + str(strandNum))
        strandNum += 1
    
    directors = [ (pos - centerPos)/np.linalg.norm(pos - centerPos) for pos in Polymer1SpherePoints] 
    
    allNames = []
    
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
    
    fIO.saveXYZList(xyzVals, allNames, filename + "_sphere.xyz")


if __name__=="__main__":

    
    bondLength = 1.5
    numA = 15
    numB = 135
    numMonomersPerPolymer = numA + numB
    numPolymersPerSphere = 180 
    FMaxRadius = 40
    FMinRadius = 20
    FZ1 = 500
    FZ2 = 100
    alpha1= 40
    alpha2 = 80
    beta1= 130
    beta2 = 180
    AtomicMinDist = 1.0

    makeSphereAndPolymers(numA, numB, numPolymersPerSphere, FMaxRadius, FMinRadius, FZ1, FZ2, alpha1, alpha2, beta1, beta2 , AtomicMinDist, bondLength, "A_15_B_135")
    makeSphereAndPolymers(numB, numA, numPolymersPerSphere, FMaxRadius, FMinRadius, FZ1, FZ2, alpha1, alpha2, beta1, beta2 , AtomicMinDist, bondLength, "A_135_B_15")


print "example done"