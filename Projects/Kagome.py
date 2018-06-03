import numpy as np
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
#from Library.VolumePackSquareBasedPyramid import VolumePackSquareBasedPyramidBBG as VPSBPBBG
from Library.VolumePackCuboidSCP import VolumePackCuboidSCParticlesBBG as VPCSCPBBG
import Utilities.fileIO as fIO


def computeKagomeLattice(NumAX, NumAY, NumAZ, KagomeAngle, AX, AY, AZ, BX, BY, BZ):
    # initialise the arrays
    APosList = []
    BPosList = []
    ADirList = []
    BDirList = []
    ADir = np.array([0.0, 0.0, 1.0])
    BDir = np.array([BZ * np.sin(KagomeAngle), 0.0, BZ * np.cos(KagomeAngle)])

    # values that never change (relative intervals between various points in the lattice)    
    BDiag = np.sqrt(np.power((BX),2) + np.power((BZ),2))
    BAngle = np.arctan2(BZ, BX)
    X_AA= AX + BX * np.cos(KagomeAngle)
    Y_AA = AY + BY 
    Z_AA = AZ + BZ * np.cos(KagomeAngle)
    X_AB = AX/2.0 + BDiag * np.sin(KagomeAngle + BAngle)/2.0
    Y_AB = (AY + BY)/2.0
    Z_AB = AZ/2.0 + BDiag * np.cos(KagomeAngle + BAngle)/2.0
    DX = BZ * np.sin(KagomeAngle)
    DZ = BX * np.sin(KagomeAngle)

    # Generate List of X, y, Z positions of the A cubes and B cubes in separate lists
    # At the same time generate a list of the directions of each cube.    
    for XCube in range(NumAX):
        for YCube in range(NumAY):
            for ZCube in range(NumAZ):
            
                if len(APosList)==0:
                    APosList = [ np.array( [0.0, 0.0, 0.0] ) ]
                    BPosList = [ np.array( [X_AB, Y_AB, Z_AB] ) ]
                    ADirList = [ ADir ]
                    BDirList = [ BDir ]
                    
                else:
                    APosList.append( [XCube * X_AA + ZCube * DX, YCube * Y_AA, ZCube * Z_AA - XCube * DZ ] )
                    BPosList.append( [XCube * X_AA + ZCube * DX + X_AB, YCube * Y_AA + Y_AB, ZCube * Z_AA - XCube * DZ + Z_AB ] )
                    ADirList.append(ADir)
                    BDirList.append(BDir)

    return APosList, BPosList, ADirList, BDirList  


# generate xyz positions and names for all the monomers in the all the polymers in all the cubes of a given type.
def polyCubes(PosList, 
              DirList, 
              numPolymersPerCube,
              X,
              Y,
              Z,
              thetad1,
              thetad2,
              phid1,
              phid2,
              polymerCylinderDiam,
              polymerCylinderLength,
              numSpheresPerParticle, 
              numMonomers,
              Alpha1,
              Alpha2,
              Beta1, 
              Beta2,
              PolyGen,
              filename):

    # compute some universal vals
    envelopeList = ['frustum ' + str(polymerCylinderLength) + ' ' + str(polymerCylinderDiam/2.0) + ' ' +str(-polymerCylinderDiam) + ' ' +str(polymerCylinderDiam)]
    polymerStartPoint = np.array([0.0, 0.0, polymerCylinderDiam])
    bondLength = polymerCylinderLength/numMonomers 
    
    curStrand = 0
    for Pos, Dirn in zip(PosList, DirList):
        # generate a building block containing the positions and orientations of the polymer envelopes.
        Cuboid_BB = CuboidPackSCBBG.generateBuildingBlock(numPolymersPerCube, X, Y, Z, thetad1, thetad2, phid1, phid2, polymerCylinderDiam, polymerCylinderLength, numSpheresPerParticle)
        
        # transform the positions and director of the building block to the current cube in the lattice. 
        Cuboid_BB.transformBBToLabFrame(Dirn, Pos, 0.0)

        # extract the start position and director for each polymer from the building block
        PolyStartPoints = Cuboid_BB.blockXYZVals[0:-1:2] 
        PolyDirectors = [ (PosB - PosA)/np.linalg.norm(PosB - PosA) for PosA, PosB in zip(Cuboid_BB.blockXYZVals[0:-1:2], Cuboid_BB.blockXYZVals[1:-1:2])] 

        # compute a polymer for each position, director pair in the building block. build the list across all the cubes.        
        for director, pos in zip(PolyDirectors, PolyStartPoints):
            print curStrand, " out of ", len(PolyStartPoints) * len(PosList), " polymers of type ", filename 
            # compute a polymer that fits in a simple envelope with the given polymer parameters
            Polymer = PolyGen.generateBuildingBlock( numMonomers, 
                                                      polymerStartPoint,
                                                      Alpha1,
                                                      Alpha2,
                                                      Beta1, 
                                                      Beta2,
                                                      polymerCylinderDiam,
                                                      bondLength,
                                                      envelopeList=envelopeList,
                                                      visualiseEnvelope=(0, 100))    
            Polymer.setBlockRefPoint(polymerStartPoint)
            # transform the polymer to the correct positions and orientation 
            Polymer.transformBBToLabFrame(director, pos, 0.0)
            
            # append the xyzvals and names of each sphere to a mega list
            if curStrand==0:
                xyzVals = Polymer.blockXYZVals
                allNames = Polymer.blockAtomNames
            else:
                xyzVals = np.concatenate((xyzVals, Polymer.blockXYZVals), 0)
                allNames = np.concatenate( (allNames, Polymer.blockAtomNames), 0)
            curStrand += 1
        
    fIO.saveXYZList(xyzVals, allNames, filename)


def dumpCubes(PosList, X, Y, Z, KagomeAngle, filename, particleName):
    posOut = []
    curPos = 0
    R = np.matrix([[np.cos(KagomeAngle), 0, - np.sin(KagomeAngle)],
                  [0 , 1, 0],
                  [np.sin(KagomeAngle), 0,  np.cos(KagomeAngle)] ])
    
    Diag = np.sqrt(np.power(X, 2) + np.power(Z, 2))/2.0
    BAngle = np.arctan2(X, Z)
    
    P1 = np.array([ Diag * np.sin(BAngle), 0,  Diag * np.cos(BAngle)])
    P2 = np.array([ Diag * np.cos(BAngle), 0, -Diag * np.sin(BAngle)])
    P3 = np.array([-Diag * np.sin(BAngle), 0, -Diag * np.cos(BAngle)])
    P4 = np.array([-Diag * np.cos(BAngle), 0,  Diag * np.sin(BAngle)])
    
    P1R = P1.dot(R).A1
    P2R = P2.dot(R).A1
    P3R = P3.dot(R).A1
    P4R = P4.dot(R).A1
    
    for pos in PosList:
        if curPos == 0:
            posOut =     [pos + P1R + np.array([0, Y/2, 0])]
            posOut.append(pos + P2R + np.array([0, Y/2, 0]))
            posOut.append(pos + P3R + np.array([0, Y/2, 0]))
            posOut.append(pos + P4R + np.array([0, Y/2, 0]))
            posOut.append(pos + P1R - np.array([0, Y/2, 0]))
            posOut.append(pos + P2R - np.array([0, Y/2, 0]))
            posOut.append(pos + P3R - np.array([0, Y/2, 0]))
            posOut.append(pos + P4R - np.array([0, Y/2, 0]))
        else:
            posOut.append(pos + P1R + np.array([0, Y/2, 0]))
            posOut.append(pos + P2R + np.array([0, Y/2, 0]))
            posOut.append(pos + P3R + np.array([0, Y/2, 0]))
            posOut.append(pos + P4R + np.array([0, Y/2, 0]))
            posOut.append(pos + P1R - np.array([0, Y/2, 0]))
            posOut.append(pos + P2R - np.array([0, Y/2, 0]))
            posOut.append(pos + P3R - np.array([0, Y/2, 0]))
            posOut.append(pos + P4R - np.array([0, Y/2, 0]))            
        curPos += 1

    fIO.saveXYZ(posOut, particleName, filename)

if __name__=="__main__":

    # create the NPack object.
    CuboidPackSCBBG = VPCSCPBBG("VolumePackCuboidSC.txt")
    PolymerBBG = RPBBG("RandomPolymer.txt")

    # lattice parameters
    NumAX = 10
    NumAY = 10
    NumAZ = 3
    KagomeAngle = 30.0 * np.pi/180.0
    AX = 15.0
    AY = 15.0
    AZ = 15.0
    BX = 8.0
    BY = 8.0
    BZ = 8.0

    print "Computing Lattice Positions"

    # compute positions and directions of each cube in the super array
    APosList, BPosList, ADirList, BDirList = computeKagomeLattice(NumAX, NumAY, NumAZ, KagomeAngle, AX, AY, AZ, BX, BY, BZ)
    
    dumpCubes(APosList, AX, AY, AZ, 0.0, 'ACubeCorners.xyz', 'H')
    dumpCubes(BPosList, BX, BY, BZ, KagomeAngle, 'BCubeCorners.xyz', 'O')
    fIO.saveXYZ(APosList, 'Ca', 'ACubeCenters.xyz')
    fIO.saveXYZ(BPosList, 'O', 'BCubeCenters.xyz')
    
    # cube internal parameters    
    numPolymersPerCubeA = 120
    thetad1A = -5.0
    thetad2A = 5.0
    phid1A = -10.0
    phid2A = 10.0
    polymerCylinderDiamA = 2.0
    polymerCylinderLengthA = 0.8 * AX
    numSpheresPerParticleA = 2

    numPolymersPerCubeB = 200
    thetad1B = -75.0
    thetad2B = -65.0
    phid1B = -10.0
    phid2B = 10.0
    polymerCylinderDiamB = 1.0
    polymerCylinderLengthB = BZ/2.0
    numSpheresPerParticleB = 2
    
    # set up parameters for the A type polymers    
    Alpha1A = 40.0
    Alpha2A = 60.0
    Beta1A = 110.0
    Beta2A = 120.0
    numMonomersA = 30
    
    # set up parameters for the B type polymers    
    Alpha1B = 0.0
    Alpha2B = 90.0
    Beta1B = 170.0
    Beta2B = 180.0
    numMonomersB = 15
    
    # generate a file for the A pos list
    PolymerBBG.particleName = 'C'
    polyCubes( APosList, 
               ADirList, 
               numPolymersPerCubeA,
               AX,
               AY,
               AZ,
               thetad1A,
               thetad2A,
               phid1A,
               phid2A,
               polymerCylinderDiamA,
               polymerCylinderLengthA,
               numSpheresPerParticleA, 
               numMonomersA,
               Alpha1A,
               Alpha2A,
               Beta1A, 
               Beta2A,
               PolymerBBG,
               'polyCubeA.xyz')
    
    # generate the B poly cubes
    PolymerBBG.particleName = 'O'
    polyCubes( BPosList, 
               BDirList, 
               numPolymersPerCubeB,
               BX,
               BY,
               BZ,
               thetad1B,
               thetad2B,
               phid1B,
               phid2B,
               polymerCylinderDiamB,
               polymerCylinderLengthB,
               numSpheresPerParticleB, 
               numMonomersB,
               Alpha1B,
               Alpha2B,
               Beta1B, 
               Beta2B,
               PolymerBBG,
               'polyCubeB.xyz')

print "Kagome done"