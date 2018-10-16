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
    AtomicMinDist = 1.0
    numA = 15
    numB = 135
    numMonomersPerPolymer = numA + numB
    numPolymersPerSphere_A = 50 
    FMaxRadius_A = 20
    FMinRadius_A = 2
    FZ1_A = 500
    FZ2_A = 5
    alpha1_A= 40
    alpha2_A = 60
    beta1_A = 165
    beta2_A = 185
    

    # makeSphereAndPolymers(numA, numB, numPolymersPerSphere_A, FMaxRadius_A, FMinRadius_A, FZ1_A, FZ2_A, alpha1_A, alpha2_A, beta1_A, beta2_A, AtomicMinDist, bondLength, "A_15_B_135")
    
    
    numPolymersPerSphere_B = 73
    FMaxRadius_B = 20
    FMinRadius_B = 4
    FZ1_B = 150
    FZ2_B = 12
    alpha1_B= 40
    alpha2_B = 50
    beta1_B = 145
    beta2_B = 155
    
    makeSphereAndPolymers(numB, numA, numPolymersPerSphere_B, FMaxRadius_B, FMinRadius_B, FZ1_B, FZ2_B, alpha1_B, alpha2_B, beta1_B, beta2_B , AtomicMinDist, bondLength, "A_135_B_15")


    print "example done"