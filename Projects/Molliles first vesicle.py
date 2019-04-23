import numpy as np
import copy as cp
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO


if __name__=="__main__":
    SphereBBG = SPSBBG('SurfacePackSphere.txt')
    PolymerGenerator = RPBBG('RandomPolymer.txt')

    numPoints  = 100
    radius = 20
    theta1 = -90
    theta2 = 90
    phi1 = -180
    phi2 = 180
    minDist = 1

    baseSPhereBB = SphereBBG.generateBuildingBlock(numPoints, radius, theta1, theta2, phi1, phi2, minDist)
    fIO.saveXYZ(baseSPhereBB.blockXYZVals, 'C', "mollieShere.xyz")


    numMonomers  = 100
    pointA = np.array([0.0, 0.0, 0.0]) 
    
    alpha1 = 40
    alpha2 = 50
    beta1 = 110 
    beta2 = 130
    minDist = 1.0
    bondLength = 2.0
    
    polymers = [ PolymerGenerator.generateBuildingBlock(numMonomers, 
                                                        pointA, 
                                                        alpha1, 
                                                        alpha2, 
                                                        beta1, 
                                                        beta2, 
                                                        minDist, 
                                                        bondLength) for _ in range(0, numPoints)]

    outputXYZ = []
    # hello
    for unimer, basePoint in zip( polymers, baseSPhereBB.blockXYZVals ):
        
        director = basePoint / np.linalg.norm(basePoint) 
        
        newUnimerXYZ = coords.transformFromBlockFrameToLabFrame( director, 
                                                                 basePoint, 
                                                                 0.0, 
                                                                 np.array([0.0, 0.0, 1.0]), 
                                                                 pointA, 
                                                                 unimer.blockXYZVals)
        
        outputXYZ += cp.copy(newUnimerXYZ)
            
    fIO.saveXYZ(outputXYZ, 'C', "unimers.xyz")
    
    
    print("Example done")




