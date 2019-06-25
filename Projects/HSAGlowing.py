import numpy as np
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO

from Library.XYZBuildingBlock import XYZBuildingBlockGenerator as XYZBBG
from Library.VolumePackCuboid import VolumePackCuboidBBG as VPCBBG



def makeHSADistribution(num0, num1, num2, num3, num4, num5, x, y, z, minDist):

    # create the cubic distribution
    cubeGen = VPCBBG('VolumePackCuboid.txt')

    # create the NPack object.
    numPoints = num0 + num1 + num2 + num3 + num4 + num5 
    
    # generate the building block
    SpacePackBB = cubeGen.generateBuildingBlock(numPoints, -x/2, x/2, -y/2, y/2, -z/2, z/2, minDist)

    # dump the base points to file    
    fIO.saveXYZ(SpacePackBB.blockXYZVals, 'CA', 'points.xyz')
    
    # generate sufficient random orientiation vectors
    thetaPhiArray = [coords.pickRandomPointOnUnitSphere() for _ in range(numPoints)]
    directors =[ coords.sphericalPolar2XYZ(np.array([1.0, angs[0], angs[1]])) for angs in thetaPhiArray ]
    
    xyzVals0, names0 = genHsaPoints('hsa_apo.xyz', SpacePackBB.blockXYZVals[0:num0], directors[0:num1] )
    xyzVals1, names1 = genHsaPoints('hsa_with1.xyz', SpacePackBB.blockXYZVals[num0:num0 + num1], directors[num0:num0 + num1])
    xyzVals2, names2 = genHsaPoints('hsa_with2.xyz', SpacePackBB.blockXYZVals[num0 + num1:num0 + num1 + num2], directors[num0 + num1:num0 + num1 + num2])
    xyzVals3, names3 = genHsaPoints('hsa_with3.xyz', SpacePackBB.blockXYZVals[num0 + num1 + num2:num0 + num1 + num2 + num3], directors[num0 + num1 + num2:num0 + num1 + num2 + num3])
    xyzVals4, names4 = genHsaPoints('hsa_with4.xyz', SpacePackBB.blockXYZVals[num0 + num1 + num2 + num3:num0 + num1 + num2 + num3 + num4], directors[num0 + num1 + num2 + num3:num0 + num1 + num2 + num3 + num4])
    xyzVals5, names5 = genHsaPoints('hsa_with5.xyz', SpacePackBB.blockXYZVals[num0 + num1 + num2 + num3 + num4:], directors[num0 + num1 + num2 + num3 + num4:])
    
    fIO.saveXYZList(xyzVals0, names0, 'hsa0All.xyz')
    fIO.saveXYZList(xyzVals1, names1, 'hsa1All.xyz')
    fIO.saveXYZList(xyzVals2, names2, 'hsa2All.xyz')
    fIO.saveXYZList(xyzVals3, names3, 'hsa3All.xyz')
    fIO.saveXYZList(xyzVals4, names4, 'hsa4All.xyz')
    fIO.saveXYZList(xyzVals5, names5, 'hsa5All.xyz')


def genHsaPoints(filename, xyzVals, directors):
    
    print("working on file:", filename)
    
    # create the pdb object generator
    XYZObjGen = XYZBBG()
    
    # create the hsa object
    hsa = XYZObjGen.generateBuildingBlock(filename) 

    xyzValsOut = []
    xyzNames = []
    
    # minVals = np.min(hsa.blockXYZVals, 0)
    # maxVals = np.max(hsa.blockXYZVals, 0)
    
    for pos, dirn in zip(xyzVals, directors):
        xyzValsOut += coords.transformFromBlockFrameToLabFrame(dirn, pos, 0, hsa.blockDirectorHat, hsa.blockRefPoint, hsa.blockXYZVals)
        newNames = [ 'W' if name=='W' else 'C' for name in hsa.blockAtomNames ]
        xyzNames += newNames
    
    return xyzValsOut, xyzNames


if __name__=="__main__":

    num0 = 16 
    num1 = 16 
    num2 = 16
    num3 = 16
    num4 = 16 
    num5 = 16 
    x = 2500 
    y = 2400 
    z = 2900 
    minDist = 250

    makeHSADistribution(num0, num1, num2, num3, num4, num5, x, y, z, minDist)    
    print("xyz building block done")