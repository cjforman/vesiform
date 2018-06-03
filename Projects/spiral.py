import numpy as np
# from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
from Library.constrainedPolymer import ConstrainedPolymerPackBBG as CPBBG
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO


def Spiral(numArms, numParticlesArm, armEndRad,  particleName, polyGen, spiralAngularDisplacement, bondLength, alpha1, alpha2, beta1, beta2, AtomicMinDist, filename):

    # Sets up an N pointed star where the initial director at point A points out to a point, 
    # and the end of the polymer is attracted to a pointB which is an angular displacement away from the 
    # initial director point. Can be quite a large angle. Generates a sort of a spiral.
    # Don't bother with envelopes and self-intersection checking. Choose your angular ranges carefully.

    directors = [  coords.sphericalPolar2XYZ(np.array([1.0, 0.0, phi])) for phi in np.linspace(0, 2*np.pi, numArms, endpoint=False) ] 
    pointsB = [  coords.sphericalPolar2XYZ(np.array([armEndRad, 0.0, phi + spiralAngularDisplacement])) for phi in np.linspace(0, 2*np.pi, numArms, endpoint=False) ] 

    numCrankMoves= 0

    polymerStartPoint = np.array([0.0, 0.0, 0.0])
    strands = [ polyGen.generateBuildingBlock(numParticlesArm, 
                                              polymerStartPoint,
                                              pointB,
                                              AtomicMinDist,
                                              bondLength,
                                              numCrankMoves,
                                              visualiseEnvelope=(0, 100),
                                              angularRange = [alpha1, alpha2, beta1, beta2],
                                              startDirector = director) for pointB, director in zip(pointsB, directors) ]

    [ strand.setBlockRefPoint(polymerStartPoint) for strand in strands ]

    names =  [particleName] * numParticlesArm 
    
    strandNum = 0
    for strand in strands:
        strand.blockAtomNames = names[:]
        strandNum += 1


    curStrand = 0
    for director, strand in zip(directors, strands):
        strand.transformBBToLabFrame(director, polymerStartPoint, 0.0)
        if curStrand==0:
            xyzVals = strand.blockXYZVals
            allNames = strand.blockAtomNames
        else:
            for xyzVal in strand.blockXYZVals:
                xyzVals.append(xyzVal)
            for name in strand.blockAtomNames:
                allNames.append(name)
        curStrand += 1
    
    fIO.saveXYZList(xyzVals, allNames, filename)
    
    return xyzVals, allNames

if __name__=="__main__":
    
    polyGen = CPBBG("ConstrainedPolymer.txt")
    
    numArms = 5
    numParticlesArm = 40
    bondLength = 1.5
    armEndRad = 0.6 * numParticlesArm * bondLength 
    particleName = 'C'
    spiralAngularDisplacement = 60 * np.pi/180
    alpha1 = 45 * np.pi/180
    alpha2 = 55 * np.pi/180
    beta1 = 170 * np.pi/180
    beta2 = 180 * np.pi/180
    atomicMinDist = 1.0
    filename = "spiral.xyz"
    Spiral(numArms, numParticlesArm, armEndRad, particleName, polyGen, spiralAngularDisplacement, bondLength, alpha1, alpha2, beta1, beta2, atomicMinDist, filename)    

print "example done"