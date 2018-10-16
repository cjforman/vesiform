import numpy as np
import random as rnd
from Projects.PolymerBrush import polymerBrush as PB
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.VolumePackEllipsoid import VolumePackEllipsoidBBG as VPEBBG
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO

def Framboise(OuterRadius, InnerRadius, MicelleRadius, numOuter, numInner, numMicelles, numPointsMicelle, 
              minDistInner, minDistOuter, minDistIntraMicelle, minDistExtraMicelle,
              theta1, theta2, phi1, phi2, surfaceSphereFilename, volumeSphereFilename, micelleBrush, vesicleBrush):

    # defines an outer vesicle and series of inner micelles and populates each type of point in these lists with a suitably defined 
    # polymer brush.
    
    totalNumberOfBrushes = numMicelles * numPointsMicelle + numOuter + numInner
    
    
    # define the main outer vesicle points and directors
    vesicleData = vesicle(OuterRadius, InnerRadius, numOuter, numInner, minDistInner, minDistOuter, theta1, theta2, phi1, phi2, 
                          np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]), 0.0, surfaceSphereFilename, nameInner='C', nameOuter='C')


    # define centers of micelles 
    envelopeList = ['innersphere ' + str(InnerRadius + 0.7 * minDistExtraMicelle),
                    'outersphere ' + str(OuterRadius - 0.7 * minDistExtraMicelle)]
    
    spherePack = VPEBBG(volumeSphereFilename)
    micelleCentrePointsBB = spherePack.generateBuildingBlock(  numMicelles, 
                                                         OuterRadius, 
                                                         OuterRadius, 
                                                         OuterRadius, 
                                                         theta1, 
                                                         theta2, 
                                                         phi1, 
                                                         phi2, 
                                                         minDistExtraMicelle, 
                                                         envelopeList=envelopeList,
                                                         defaultBlockRefPoint = np.array([0.0, 0.0, 0.0]))

    # generate the points, directors for each micelle
    micelleData = []
    for micelleCentrePoint in micelleCentrePointsBB.blockXYZVals:
        newMicelle = micelle( MicelleRadius, numPointsMicelle, minDistIntraMicelle, 
                              -90, 90, -180, 180, 
                              micelleCentrePoint, np.array([0.0, 0.0, 1.0]), 0.0, surfaceSphereFilename, name='C') 
        
        micelleData.append(newMicelle)

    # generate polymerBrush for each point in each micelle and transfer to framboiseXYZ and framboiseNames 
    framboiseXYZ = []
    framboiseNames = []
    firstTime = True;
    i = 0
    for curMicelle in micelleData:
        for point, dirNorm in zip(curMicelle[0], curMicelle[2]):
            if firstTime:
                (brushDir, brushRefPoint, brushXyzVals, brushNames) = PB(micelleBrush)
                newXyzVals = coords.transformFromBlockFrameToLabFrame(dirNorm, point, rnd.uniform(0, 2 * np.pi), brushDir, brushRefPoint, brushXyzVals)
                framboiseXYZ = newXyzVals
                framboiseNames = brushNames
                firstTime = False
            else:
                (brushDir, brushRefPoint, brushXyzVals, brushNames) = PB(micelleBrush)
                newXyzVals = coords.transformFromBlockFrameToLabFrame(dirNorm, point, rnd.uniform(0, 2 * np.pi), brushDir, brushRefPoint, brushXyzVals)
                framboiseXYZ += newXyzVals
                framboiseNames = np.concatenate( (framboiseNames, brushNames), 0 )
            i+=1
            if i%10==0:
                print "Completed Brush " + str(i) + " of " + str(totalNumberOfBrushes)

    # set the peptide brushnames to B (enables different colour to be used for centre of vesicles
    framboiseNames = [ 'B' if a=='C' else a for a in framboiseNames]

    # for each point in the vesicle create a polymer brush
    for point, dirNorm in zip(vesicleData[0], vesicleData[2]):
        (brushDir, brushRefPoint, brushXyzVals, brushNames) = PB(vesicleBrush)
        newXyzVals = coords.transformFromBlockFrameToLabFrame(dirNorm, point, rnd.uniform(0, 2 * np.pi), brushDir, brushRefPoint, brushXyzVals)
        framboiseXYZ += newXyzVals
        framboiseNames = np.concatenate( (framboiseNames, brushNames), 0 )
        i+=1
        if i%10==0:
            print "Completed Brush " + str(i) + " of " + str(totalNumberOfBrushes)

    return (framboiseXYZ,list(framboiseNames))

def vesicle(OuterRadius, InnerRadius, numOuter, numInner, minDistInner, minDistOuter, theta1, theta2, phi1, phi2, centrePos, zDirector, omega, filename, nameInner=None, nameOuter=None):
    # Defines a bunch of points in an inner sphere and an outer sphere, 
    # Returns the points on the spherical surface, the atomnames and a list of radial vectors corresponding to each point 
    # which is the direction from that point to or from the centre depending on inner or outer status.
    # The inner sphere has Radial points pointing away from centre 
    # The outer sphere has radial points point towards the centre.
    # No problem with making the inner outside the outer if that makes sense.
    # The points are in a domain defined in spherical polars allowing a section to be removed from the sphere for the classic vesicle look.  
    # The packing parameter on each surface is defined by minDistInner and minDistOuter which is the min allowed distance between points on the inner and outer respectively.
    # Translates the sphere to centrePos with the Z pole of sphere aligned with director and an angle omega rotation about that axis.
    # Filename defines the packing algorithm parameters. (How hard to try, name of atom in the structure etc).
    # Name can be overridden by the name parameters 
    SphereBBG = SPSBBG(filename)
    
    # generate the outer sphere
    outerSphereBB = SphereBBG.generateBuildingBlock(numOuter, OuterRadius, theta1, theta2, phi1, phi2, minDistOuter, visualiseEnvelope=(100000, 1000, 'sphereEnvelope.xyz'))
    outerSphereBB.transformBBToLabFrame(zDirector, centrePos, omega)
    outerSphereXYZ = outerSphereBB.blockXYZVals
    outerDirectors = [ centrePos - a for a in outerSphereXYZ]
    outerDirectorsNorm = [ a/np.linalg.norm(a) for a in outerDirectors]
    if nameOuter=='None':
        namesOuter = outerSphereBB.blockAtomNames
    else:
        namesOuter = [nameOuter] * numOuter
    
    # generate the inner sphere
    innerSphereBB = SphereBBG.generateBuildingBlock(numInner, InnerRadius, theta1, theta2, phi1, phi2, minDistInner)
    innerSphereBB.transformBBToLabFrame(zDirector, centrePos, omega)
    innerSphereXYZ = innerSphereBB.blockXYZVals
    innerDirectors = [ a - centrePos for a in innerSphereXYZ]
    innerDirectorsNorm = [ a/np.linalg.norm(a) for a in innerDirectors]
    
    if nameInner=='None':
        namesInner = innerSphereBB.blockAtomNames
    else:
        namesInner = [nameInner] * numInner

    vesXYZVals = outerSphereXYZ + innerSphereXYZ
    vesNames = namesOuter + namesInner     
    vesDirs = outerDirectorsNorm + innerDirectorsNorm     

    return (vesXYZVals, vesNames, vesDirs) 


def micelle(radius, num, minDist, theta1, theta2, phi1, phi2, centrePos, zDirector, omega, filename, name=None):
    # Defines a bunch of points in a sphere, 
    # Returns the points on the spherical surface, the atomnames and a list of radial vectors point from the surface to the centre. 
    # The points are in a domain defined in spherical polars allowing a section to be removed from the sphere for the classic look.  
    # The packing parameter on the surface is defined by minDist which is the min allowed distance between points.
    # Translates the micelle to centrePos with the Z pole of sphere aligned with director and an angle omega rotation about that axis.
    # Filename defines the packing algorithm parameters. (How hard to try, name of atom in the structure etc).
    # Name can be overridden by the name parameter 
    SphereBBG = SPSBBG(filename)

    # generate the inner sphere
    SphereBB = SphereBBG.generateBuildingBlock(num, radius, theta1, theta2, phi1, phi2, minDist)
    SphereBB.transformBBToLabFrame(zDirector, centrePos, omega)
    SphereXYZ = SphereBB.blockXYZVals
    Directors = [ a - centrePos for a in SphereXYZ]
    DirectorsNorm = [ a/np.linalg.norm(a) for a in Directors]
    
    if name=='None':
        names = SphereBB.blockAtomNames
    else:
        names = [name] * num

    return (SphereXYZ, names, DirectorsNorm) 


if __name__=="__main__":

    # first set up the brush polymer parameters
    brushDict={}
    # pack the brush dictionary
    brushDict['filenameBlock'] = 'RandomPolymer.txt'  
    brushDict['filenameBrush'] = 'RandomPolymer.txt' # use 'betastrand.txt' for mode peptide
    brushDict['mode'] = 'Polymer'
    brushDict['ABlock'] = {}
    brushDict['BBlock'] = {}
    brushDict['brushBlock'] = {}
    

    brushDict['ABlock']['num'] = 5 
    brushDict['ABlock']['alpha1'] = 40
    brushDict['ABlock']['alpha2'] = 50
    brushDict['ABlock']['beta1'] = 165
    brushDict['ABlock']['beta2'] = 185
    brushDict['ABlock']['minDist'] = 1.0
    brushDict['ABlock']['bondLength'] = 1.5
    brushDict['ABlock']['Z1'] = 2
    brushDict['ABlock']['R1'] = 20
    brushDict['ABlock']['Z2'] = brushDict['ABlock']['num'] * brushDict['ABlock']['bondLength'] + brushDict['ABlock']['Z1']
    brushDict['ABlock']['R2'] = 20


    brushDict['BBlock']['num'] = 40 
    brushDict['BBlock']['alpha1'] = 40.0
    brushDict['BBlock']['alpha2'] = 50.0
    brushDict['BBlock']['beta1'] = 145.0
    brushDict['BBlock']['beta2'] = 155.0
    brushDict['BBlock']['minDist'] = 1.0
    brushDict['BBlock']['bondLength'] = 1.5
    brushDict['BBlock']['Z1'] = 2
    brushDict['BBlock']['R1'] = 20
    brushDict['BBlock']['Z2'] = brushDict['BBlock']['num'] * brushDict['BBlock']['bondLength'] + brushDict['BBlock']['Z1']   
    brushDict['BBlock']['R2'] = 20
    
    brushDict['brushBlock']['num'] = 15
    brushDict['brushBlock']['alpha1'] = 40.0
    brushDict['brushBlock']['alpha2'] = 50.0
    brushDict['brushBlock']['beta1'] = 145.0
    brushDict['brushBlock']['beta2'] = 155.0
    brushDict['brushBlock']['minDist'] = 1.0
    brushDict['brushBlock']['bondLength'] = 1.5
    brushDict['brushBlock']['Z1'] = 0.0
    brushDict['brushBlock']['R1'] = 10.0
    brushDict['brushBlock']['Z2'] = brushDict['brushBlock']['num'] * brushDict['brushBlock']['bondLength'] + brushDict['brushBlock']['Z1']
    brushDict['brushBlock']['R2'] = 20.0
    
    
    
    # Set up the Framboise Parameters
    numMicelles = 65
    numPointsMicelle = 160
    MicelleRadius = 3.5 * brushDict['ABlock']['num'] * brushDict['ABlock']['bondLength']
    minDistIntraMicelle = 0.25 * brushDict['brushBlock']['num'] * brushDict['brushBlock']['bondLength']
    blockBLength = 0.8 * (brushDict['BBlock']['num'] - 1) * brushDict['BBlock']['bondLength']
    minDistExtraMicelle = 2.0 * blockBLength   +  MicelleRadius
    theta1 = -90.0
    theta2 = 90
    phi1 = -180.0
    phi2 = 90.0
    numInner = 2000
    InnerRadius = 300.0
    minDistInner = minDistIntraMicelle
    numOuter = 7000
    OuterRadius = InnerRadius + 4 * blockBLength
    minDistOuter = minDistIntraMicelle
    
    surfaceSphereFilename = "SurfacePackSphere.txt"
    volumeSphereFilename = "VolumePackEllipsoid.txt"
    micelleBrush = brushDict.copy()
    micelleBrush['mode'] = "PolymerRandom"
    vesicleBrush = brushDict.copy()
    vesicleBrush['mode'] = "PolymerRandom"

    # generate the Framboise
    (framboiseXYZVals, FramboiseNames) = Framboise(OuterRadius, InnerRadius, MicelleRadius, numOuter, numInner, numMicelles, numPointsMicelle, 
                                                   minDistInner, minDistOuter, minDistIntraMicelle, minDistExtraMicelle, 
                                                   theta1, theta2, phi1, phi2, surfaceSphereFilename, volumeSphereFilename, micelleBrush, vesicleBrush)
    
    fIO.saveXYZList(framboiseXYZVals, FramboiseNames, "Framboise.xyz")

    print "Framboise done"