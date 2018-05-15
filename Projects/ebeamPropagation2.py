import numpy as np
import copy as cp
import random as rnd
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
#from Library.VolumePackSquareBasedPyramid import VolumePackSquareBasedPyramidBBG as VPSBPBBG
from Library.VolumePackCuboid import VolumePackCuboidBBG as VPCBBG
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords


def eBeamSphere(numPolymersPerSphere, numPEG, numHPMA, FZ1, FMinRadius1, FZ2, FMinRadius2, AtomicMinDist, bondLength, filename, SphereBBG, PolyGen, centerPos):

    # generate the root points in the spherical inner shell
    SphereBB = SphereBBG.generateBuildingBlock(numPolymersPerSphere, FZ2, -90, 90, -180, 180, FMinRadius2)
    SphereBB.transformBBToLabFrame(np.array([0.0, 0.0, 1.0]), centerPos, 0.0)
    SpherePoints = SphereBB.blockXYZVals 

    directors = [ (pos - centerPos)/np.linalg.norm(pos - centerPos) for pos in SpherePoints] 

    allNames = []
    
    curStrand = 0
    for director, pos in zip(directors, SpherePoints):
        print curStrand, " out of ", len(SpherePoints)
        strand = generateEBeamPolymer(numPEG, numHPMA, PolyGen, FZ1, FMinRadius1, FZ2, FMinRadius2, AtomicMinDist, bondLength)
        strand.transformBBToLabFrame(director, pos, 0.0)
        if curStrand==0:
            xyzVals = strand.blockXYZVals
            allNames = strand.blockAtomNames
        else:
            xyzVals = np.concatenate((xyzVals, strand.blockXYZVals), 0)
            allNames = np.concatenate( (allNames, strand.blockAtomNames), 0)
        curStrand += 1
    
    fIO.saveXYZList(xyzVals, allNames, filename)
    

def generateEBeamPolymer(numPEG, numHPMA, PolyGen, FZ1, FRadius1, FZ2, FRadius2, AtomicMinDist, bondLength):
    

    if numHPMA==0:
        numPEG += 1
    else:
        numHPMA += 1

    PEGAlpha1= 40.0
    PEGAlpha2 = 80.0
    PEGBeta1= 130.0
    PEGBeta2 = 165.0
    HPMAAlpha1 = -15.0 # 10.0
    HPMAAlpha2 = 45.0 # 40.0
    HPMABeta1= 115 # 120.0
    HPMABeta2 = 125 # 165.0

    envelopeList = ['frustum ' + str(FZ1) + ' ' + str(FRadius1) + ' ' +str(FZ2 - AtomicMinDist) + ' ' +str(FRadius2)]

    polymerStartPoint = np.array([0.0, 0.0, FZ2 + AtomicMinDist])

    if numHPMA>0:
        HPMAStrand = PolyGen.generateBuildingBlock(   numHPMA, 
                                                      polymerStartPoint,
                                                      HPMAAlpha1,
                                                      HPMAAlpha2,
                                                      HPMABeta1, 
                                                      HPMABeta2,
                                                      AtomicMinDist,
                                                      bondLength,
                                                      envelopeList=envelopeList,
                                                      visualiseEnvelope=(0, 100))    
        HPMAStrand.setBlockRefPoint(polymerStartPoint)


    PEGStrand = PolyGen.generateBuildingBlock(  numPEG, 
                                                polymerStartPoint,
                                                PEGAlpha1,
                                                PEGAlpha2,
                                                PEGBeta1, 
                                                PEGBeta2,
                                                AtomicMinDist,
                                                bondLength,
                                                envelopeList=envelopeList,
                                                visualiseEnvelope=(0, 100))

    PEGStrand.setBlockRefPoint(polymerStartPoint)

    if numHPMA>0:
        names = ['O'] * numHPMA
        names = np.concatenate( (names, ['N'] * numPEG), 0 )
    else:
        names = ['N'] * numPEG  
    names[0] = 'S'

    retStrand = PEGStrand
    if numHPMA>0:
        # retStrand, _ = HPMAStrand.addBuildingBlock(PEGStrand, 1, 0, bondLength, -57.0, 116.0, 180.0, 122.0, -47.0, polymerStartPoint, np.array([0.0, 0.0, 1.0]))
        retStrand = HPMAStrand
        # complete hack, just ignores connectors etc not using them anyway for this project later on.
        newXYZVals = [ HPMAStrand.blockXYZVals[-1] + np.array([0.0, 0.0, bondLength]) +  PegXYZ - PEGStrand.blockXYZVals[0] for PegXYZ in PEGStrand.blockXYZVals]
        retStrand.addAtoms((newXYZVals, PEGStrand.blockAtomNames))

    retStrand.blockAtomNames = names[:]

    return retStrand


if __name__=="__main__":

    PolyGen = RPBBG('RandomPolymer.txt')    
    SphereBBG = SPSBBG('SurfacePackSphere.txt')
    #PyramidBBG = VPSBPBBG('VolumePackPyramid.txt')
    CubeBBG = VPCBBG('VolumePackCuboid.txt')

    numPolymersPerSphere = 195
    AtomicMinDist = 1.0
    NumPEG = 60
    NumHPMA = 113
    bondLength = 1.5
    FRadius1 = 20.0
    FRadius2 = 2.0
    FZ1 = 150.0
    FZ2 = 10.0
    
    
    BaseX = 4000.0
    BaseY = 5000.0
    BaseZ = 1000.0
    Volume = BaseX * BaseY * BaseZ
    
    NoHPMAConcentration  = 2.5e-8
    NumPolysNoHPMA = int(np.floor((NoHPMAConcentration * Volume))) 
    NoHPMAMinDist = 0.8 * np.power(Volume/float(NumPolysNoHPMA), 1.0/3.0)
    
    WithHPMAConcentration  = 2.5e-8
    NumPolysWithHPMA = int(np.floor(WithHPMAConcentration * Volume))
    WithHPMAMinDist = 0.8 * np.power(Volume/float(NumPolysWithHPMA), 1.0/3.0)

    monomerConcentration  = 1e-7
    NumHPMAMonomers = int(np.floor(monomerConcentration * Volume))
    monomerMinDist = 0.5 * np.power(Volume/float(NumHPMAMonomers), 1.0/3.0)
    
    NumMicelles =  25
    micelleMinDist =  np.power(Volume/float(NumMicelles), 1.0/3.0)
    
    NumPegOnly = 5
    
    doMicelles = False
    doHPMAPolys = True
    doPEGPolys = True
    doMonomers = True
    doPegStumps = True
    
    if doMicelles:
        print "***** Constructing Micelles ****"
        # make the micelles - each is dumped as a separate file
        centrePoints = CubeBBG.generateBuildingBlock(NumMicelles, -BaseX/2.0, BaseX/2.0, -BaseY/2.0, BaseY/2.0, -BaseZ/2.0, BaseZ/2.0, micelleMinDist)
        micelleNum = 0
        for centrePoint in centrePoints.blockXYZVals:
            print "Constructing Micelle: ", micelleNum, " of ", NumMicelles
            eBeamSphere(numPolymersPerSphere, NumPEG, NumHPMA, FZ1, FRadius1, FZ2, FRadius2, AtomicMinDist, bondLength, 'Sphere_' + str(NumHPMA) + "_" + str(micelleNum) + ".xyz", SphereBBG, PolyGen, centrePoint)
            micelleNum += 1

    if doPEGPolys:
        print "***** Constructing polymers with no HMPA ****"
    
        # generate a bunch of points for the polymers that have no HPMA block and save them 
        PolysNoHPMAPoints = CubeBBG.generateBuildingBlock(NumPolysNoHPMA, -BaseX/2.0, BaseX/2.0, -BaseY/2.0, BaseY/2.0, -BaseZ/2.0, BaseZ/2.0, NoHPMAMinDist)
        polysNoHPMAXYZVals = []
        polysNoHPMANames = []
        curPoint = 0
        for point in PolysNoHPMAPoints.blockXYZVals:
            print curPoint, " of ", NumPolysNoHPMA 
            # produce said polymer and move it to where we want it
            polysNoHPMAStrand = generateEBeamPolymer(NumPEG, 0, PolyGen, FZ1, FRadius1, FZ2, FRadius2, AtomicMinDist, bondLength)
            theta, phi = coords.pickRandomPointOnUnitSphere()
            director = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi])) 
            rotation = rnd.uniform(0, 2 * np.pi)
            newXYZVals = coords.transformFromBlockFrameToLabFrame( director, 
                                                                   point, 
                                                                   rotation, 
                                                                   polysNoHPMAStrand.blockDirectorHat, 
                                                                   polysNoHPMAStrand.getCOM(), 
                                                                   polysNoHPMAStrand.blockXYZVals)
            if  curPoint==0:
                polysNoHPMAXYZVals = cp.copy(newXYZVals)
                polysNoHPMANames = cp.copy(polysNoHPMAStrand.blockAtomNames)
            else:
                polysNoHPMAXYZVals = np.concatenate( (polysNoHPMAXYZVals, cp.copy(newXYZVals)), 0 )
                polysNoHPMANames =  np.concatenate( (polysNoHPMANames, cp.copy(polysNoHPMAStrand.blockAtomNames)), 0)
            
            if curPoint<5:
                fIO.saveXYZList(polysNoHPMAStrand.blockXYZVals, polysNoHPMAStrand.blockAtomNames, "polysNoHPMAStrand_" + str(curPoint) +".xyz")
            curPoint += 1 
        
        fIO.saveXYZList(polysNoHPMAXYZVals, polysNoHPMANames, "polysNoHPMA.xyz")
    
    if doHPMAPolys:    
        print "***** Constructing polymers with HMPA ****"
        
        # generate a bunch of points for the polymers that do have HPMA block and save them 
        # PolysWithHPMAPoints = PyramidBBG.generateBuildingBlock(NumPolysWithHPMA, BaseX, BaseY, ZHeight, ZHeightApex, polymerMinDist)
        PolysWithHPMAPoints = CubeBBG.generateBuildingBlock(NumPolysWithHPMA, -BaseX/2.0, BaseX/2.0, -BaseY/2.0, BaseY/2.0, -BaseZ/2.0, BaseZ/2.0, WithHPMAMinDist)
         
        polysWithHPMAXYZVals = []
        polysWithHPMANames = []
        curPoint = 0
        for point in PolysWithHPMAPoints.blockXYZVals:
            print curPoint, " of ", NumPolysWithHPMA
            # produce said polymer and move it to where we want it
            polysWithHPMAStrand = generateEBeamPolymer(NumPEG, NumHPMA, PolyGen, FZ1, FRadius1, FZ2, FRadius2, AtomicMinDist, bondLength)
            theta, phi = coords.pickRandomPointOnUnitSphere()
            director = coords.sphericalPolar2XYZ(np.array([1.0, theta, phi])) 
            rotation = rnd.uniform(0, 2 * np.pi)
            newXYZVals = coords.transformFromBlockFrameToLabFrame(director, point, rotation, polysWithHPMAStrand.blockDirectorHat, polysWithHPMAStrand.getCOM(), polysWithHPMAStrand.blockXYZVals)
             
            if curPoint==0:
                polysWithHPMAXYZVals = cp.copy(newXYZVals)
                polysWithHPMANames = cp.copy(polysWithHPMAStrand.blockAtomNames)
            else:
                polysWithHPMAXYZVals = np.concatenate( (polysWithHPMAXYZVals, cp.copy(newXYZVals)), 0 )
                polysWithHPMANames =  np.concatenate( (polysWithHPMANames, cp.copy(polysWithHPMAStrand.blockAtomNames)), 0)

            if curPoint<5:
                fIO.saveXYZList(polysWithHPMAStrand.blockXYZVals, polysWithHPMAStrand.blockAtomNames, "polysWithHPMAStrand_" + str(curPoint) +".xyz")

            curPoint += 1
        
        fIO.saveXYZList(polysWithHPMAXYZVals, polysWithHPMANames, "polysWithHPMA.xyz")

    if doMonomers:    
        print "***** Constructing monomers ******"
        # drop a load of monomers to file
        # HPMAMonomerPoints = PyramidBBG.generateBuildingBlock(NumHPMAMonomers, BaseX, BaseY, ZHeight, ZHeightApex, monomerMinDist)
        HPMAMonomerPoints = CubeBBG.generateBuildingBlock(NumHPMAMonomers, -BaseX/2.0, BaseX/2.0, -BaseY/2.0, BaseY/2.0, -BaseZ/2.0, BaseZ/2.0, monomerMinDist)
        fIO.saveXYZ(HPMAMonomerPoints.blockXYZVals, 'O', "HPMAMonomers.xyz")

    if doPegStumps:
        print "***** Constructing PEG Only polymers ****"
        # generate the PEG Only polymers (no sulphur blob)
        for PegStump in range(NumPegOnly):
            PegStumpStrand = generateEBeamPolymer(NumPEG - 1, 0, PolyGen, FZ1, FRadius1, FZ2, FRadius2, AtomicMinDist, bondLength)
            fIO.saveXYZ(PegStumpStrand.blockXYZVals, 'N', "PegStump_" + str(PegStump) + ".xyz")

    fIO.saveXYZList( [np.array([0.0, 0.0, 0.0])], ['S'], "Suplhur.xyz")

print "ebeam done"