import numpy as np
import random as rnd
import copy as cp
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
from Library.constrainedPolymer import ConstrainedPolymerPackBBG as CPBBG
#from Library.VolumePackSquareBasedPyramid import VolumePackSquareBasedPyramidBBG as VPSBPBBG
from Library.VolumePackCuboid import VolumePackCuboidBBG as VPCBBG
from Projects.spiral  import Spiral
#from Projects.Kagome import computeKagomeLattice as CKL
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords


if __name__=="__main__":

    # create the NPack and polymer objects.
    CuboidPackBBG = VPCBBG("VolumePackCuboid.txt")
    ConstrainedPolyGen = CPBBG("ConstrainedPolymer.txt")
    PolyGen = RPBBG("RandomPolymer.txt")
    
    BoxX = 250
    BoxY = 250
    BoxZ =  70
    
    # Spiral Params
    num3Spirals = 2
    num4Spirals = 2
    num5Spirals = 2
    numStraightPolymers = 4
    numHelicalPolymers = 4
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
    

    # create the positions of all the spirals using the points to avoid facility
    spirals3XYZBB = CuboidPackBBG.generateBuildingBlock(num3Spirals, -BoxX/2.0, BoxX/2.0, -BoxY/2.0, BoxY/2.0, -BoxZ/2.0, BoxZ/2.0, 1.3 * armEndRad)
    pointsToAvoid = cp.copy(spirals3XYZBB.blockXYZVals)
    
    spirals4XYZBB = CuboidPackBBG.generateBuildingBlock(num4Spirals, -BoxX/2.0, BoxX/2.0, -BoxY/2.0, BoxY/2.0, -BoxZ/2.0, BoxZ/2.0, 1.3 * armEndRad, pointsToAvoid = pointsToAvoid)
    for point in  spirals4XYZBB.blockXYZVals:
        pointsToAvoid.append(cp.copy(point))
    
    spirals5XYZBB = CuboidPackBBG.generateBuildingBlock(num5Spirals, -BoxX/2.0, BoxX/2.0, -BoxY/2.0, BoxY/2.0, -BoxZ/2.0, BoxZ/2.0, 1.3 * armEndRad, pointsToAvoid = pointsToAvoid)
    for point in  spirals5XYZBB.blockXYZVals:
        pointsToAvoid.append(cp.copy(point))    
        
    # create the positions of all the straight polymers
    straightPolymers = CuboidPackBBG.generateBuildingBlock(numStraightPolymers, -BoxX/2.0, BoxX/2.0, -BoxY/2.0, BoxY/2.0, -BoxZ/2.0, BoxZ/2.0, 1.3 * armEndRad, pointsToAvoid = pointsToAvoid)
    for point in  straightPolymers.blockXYZVals:
        pointsToAvoid.append(cp.copy(point))    
    
    # create the positions of all the Helical polymers
    helicalPolymers = CuboidPackBBG.generateBuildingBlock(numHelicalPolymers, -BoxX/2.0, BoxX/2.0, -BoxY/2.0, BoxY/2.0, -BoxZ/2.0, BoxZ/2.0, 1.3 * armEndRad, pointsToAvoid = pointsToAvoid)
    
    # generate the 3 spirals and build up a list of xyzVals and names as we go for all the elements of the image
    spiral = 0
    for spiralPos in spirals3XYZBB.blockXYZVals:    
        filename = "uncut_spiral_3_" + str(spiral) + ".xyz"
        curXYZVals, curNamesAll = Spiral(3, numParticlesArm, armEndRad, particleName, ConstrainedPolyGen, spiralAngularDisplacement, bondLength, alpha1, alpha2, beta1, beta2, atomicMinDist, filename)
        theta, phi = coords.pickRandomPointOnUnitSphere()
        director = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi]))
        curXYZValsNew = coords.transformFromBlockFrameToLabFrame(director, spiralPos, rnd.uniform(0, 2* np.pi), np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), curXYZVals)
        if spiral == 0:
            xyzValsAll = cp.copy(curXYZValsNew) 
            namesAll = cp.copy(curNamesAll)
        else:
            for xyzVal in curXYZValsNew:
                xyzValsAll.append(cp.copy(xyzVal))
            for name in curNamesAll:
                namesAll.append(cp.copy(name))            
        spiral += 1
    
    # generate the 4 spirals
    for spiralPos in spirals4XYZBB.blockXYZVals:    
        filename = "uncut_spiral_4_" + str(spiral) + ".xyz"
        curXYZVals, curNamesAll = Spiral(4, numParticlesArm, armEndRad, particleName, ConstrainedPolyGen, spiralAngularDisplacement, bondLength, alpha1, alpha2, beta1, beta2, atomicMinDist, filename)
        theta, phi = coords.pickRandomPointOnUnitSphere()
        director = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi]))
        curXYZValsNew = coords.transformFromBlockFrameToLabFrame(director, spiralPos, rnd.uniform(0, 2* np.pi), np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), curXYZVals)
        for xyzVal in curXYZValsNew:
            xyzValsAll.append(cp.copy(xyzVal))
        for name in curNamesAll:
            namesAll.append(cp.copy(name)) 

    # generate the 5 spirals
    for spiralPos in spirals5XYZBB.blockXYZVals:    
        filename = "uncut_spiral_5_" + str(spiral) + ".xyz"
        curXYZVals, curNamesAll = Spiral(5, numParticlesArm, armEndRad, particleName, ConstrainedPolyGen, spiralAngularDisplacement, bondLength, alpha1, alpha2, beta1, beta2, atomicMinDist, filename)
        theta, phi = coords.pickRandomPointOnUnitSphere()
        director = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi]))
        curXYZValsNew = coords.transformFromBlockFrameToLabFrame(director, spiralPos, rnd.uniform(0, 2* np.pi), np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]), curXYZVals)
        for xyzVal in curXYZValsNew:
            xyzValsAll.append(cp.copy(xyzVal))
        for name in curNamesAll:
            namesAll.append(cp.copy(name)) 

    startPointPolymer = np.array([0.0, 0.0, 0.0])

    # generate the straight polymers
    for polyPos in straightPolymers.blockXYZVals:
        curPoly = PolyGen.generateBuildingBlock(numParticlesArm, startPointPolymer, alpha1*180/np.pi, alpha2*180/np.pi, beta1*180/np.pi, beta2*180/np.pi, atomicMinDist, bondLength)
        theta, phi = coords.pickRandomPointOnUnitSphere()
        director = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi]))
        curPoly.transformBBToLabFrame(director, polyPos, rnd.uniform(0, 2* np.pi))
        for xyzVal in curPoly.blockXYZVals:
            xyzValsAll.append(cp.copy(xyzVal))
        for name in curPoly.blockAtomNames:
            namesAll.append(cp.copy(name)) 

    # generate the helical polymers
    beta1Helical = 110
    beta2Helical = 120
    for helixPos in helicalPolymers.blockXYZVals:
        curPoly = PolyGen.generateBuildingBlock(numParticlesArm, startPointPolymer, alpha1 * 180/np.pi, alpha2 * 180/np.pi, beta1Helical, beta2Helical, atomicMinDist, bondLength)
        theta, phi = coords.pickRandomPointOnUnitSphere()
        director = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi]))
        curPoly.transformBBToLabFrame(director, helixPos, rnd.uniform(0, 2* np.pi))
        for xyzVal in curPoly.blockXYZVals:
            xyzValsAll.append(cp.copy(xyzVal))
        for name in curPoly.blockAtomNames:
            namesAll.append(cp.copy(name))  

    fIO.saveXYZList(xyzValsAll, namesAll, "uncut_autodermis.xyz")

print "autoDermis done"