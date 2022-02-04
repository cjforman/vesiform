import numpy as np
import random as rnd
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
from Library.VolumePackEllipsoid import VolumePackEllipsoidBBG as VPEBBG
import Utilities.coordSystems as coords  
import Utilities.fileIO as fIO

# builds up a random array of unimers in a spherical shell. 
def MelaninParticle(paramDict):

    # unpack the unimer dictionary
    filenameUnimer = paramDict['filenameUni']
    unimer = paramDict['unimer']
    unimerLength = unimer['unimerLength'] 
    unimerLengthRange = unimer['unimerLengthRange'] 
    alpha1 = unimer['alpha1']
    alpha2 = unimer['alpha2'] 
    beta1 = unimer['beta1'] 
    beta2 = unimer['beta2'] 
    minDist = unimer['minDist'] 
    bondLength = unimer['bondLength']
    Z1 = unimer['Z1'] 
    R1 = unimer['R1'] 
    Z2 = unimer['Z2'] 
    R2 = unimer['R2'] 

    # unpack the particle dictionary
    filenameMParticle = paramDict['filenameMP']    
    particle = paramDict['particle']
    numUnimers = particle['numUnimers'] 
    innerRadius = particle['innerRadius']
    outerRadius = particle['outerRadius']
    theta1 = particle['theta1']
    theta2 = particle['theta2']
    phi1 = particle['phi1']
    phi2 = particle['phi2']
    minDistParticle = particle['minDist'] 
    
    # construct a generator for the unimers
    unimerGenerator = RPPBBG(filenameUnimer)

    # construct a generator for the seedpoints
    seedpointGenerator = VPEBBG(filenameMParticle)
    
    # define the envelope where we want the unimers to live
    envelopeList = ['innersphere ' + str(innerRadius),
                    'outersphere ' + str(outerRadius)]

    seedPoints = seedpointGenerator.generateBuildingBlock(numUnimers, outerRadius, outerRadius, outerRadius, theta1, theta2, phi1, phi2, minDistParticle, envelopeList=envelopeList, defaultBlockRefPoint=np.array([0.0,0.0,0.0]))
    dirns = [ coords.pickRandomPointOnUnitSphere() for _ in range(numUnimers) ]
    seedDirectors = [ coords.sphericalPolar2XYZ(np.array([1.0, ang[0], ang[1]])) for ang in dirns ] 
    acceptedPoints= []
    polymerNames = []
    unimerCount = 0
    for seedPoint, seedDirector in zip(seedPoints.blockXYZVals, seedDirectors):

        # figure out the coords of the accepted Points in the reference frame of the block 
        rotPoints = coords.transformFromLabFrameToBlockFrame( np.array([0.0, 0.0, 1.0]), 
                                                              np.array([0.0, 0.0, 0.0]), 0.0, 
                                                              seedDirector, 
                                                              seedPoint,
                                                              acceptedPoints)
    
        # build a unimer avoiding the rotated accepted points
        unimer = unimerGenerator.generateBuildingBlock( rnd.choice(list(range(unimerLength - unimerLengthRange, unimerLength + unimerLengthRange + 1))),
                                                        np.array([0.0, 0.0, 0.0]),
                                                        alpha1,
                                                        alpha2,
                                                        beta1,
                                                        beta2, 
                                                        minDist,
                                                        bondLength,
                                                        pointsToAvoid = rotPoints)

            # add the new points to the rotPoints
        if len(rotPoints)==0:
            rotPoints = unimer.blockXYZVals
        else:
            rotPoints = np.concatenate((rotPoints, unimer.blockXYZVals), 0)

        if len(polymerNames)==0:            
            polymerNames = unimer.blockAtomNames
        else:
            polymerNames = np.concatenate((polymerNames, unimer.blockAtomNames), 0 )

        print("Growing list length:", len(polymerNames), " unimer ", unimerCount, " of ", numUnimers)

        # rotate the accepted points back to the lab frame
        acceptedPoints = coords.transformFromBlockFrameToLabFrame( np.array([0.0, 0.0, 1.0]), 
                                                                   np.array([0.0, 0.0, 0.0]), 0.0, 
                                                                   seedDirector, 
                                                                   seedPoint,
                                                                   rotPoints)
        
        unimerCount+=1
    
    return acceptedPoints, polymerNames

if __name__=="__main__":

    # pack the particle dictionary
    paramDict = {}
    paramDict['filenameUni'] = "RandomPolymer.txt"
    paramDict['filenameMP'] = "VolumePackEllipsoid.txt"
    unimer = {}
    unimer['unimerLength'] = 10
    unimer['unimerLengthRange'] = 3 
    unimer['alpha1'] = -20
    unimer['alpha2'] = 20
    unimer['beta1'] = 40 
    unimer['beta2'] = 60
    unimer['minDist'] = 1
    unimer['bondLength'] = 1
    unimer['Z1'] = 0
    unimer['R1'] = 3
    unimer['Z2'] = 10
    unimer['R2'] = 8
    particle = {}
    particle['numUnimers'] = 100
    particle['innerRadius'] = 50
    particle['outerRadius'] = 60
    particle['theta1'] = -90
    particle['theta2'] = 90
    particle['phi1'] = -90
    particle['phi2'] = 180
    particle['minDist'] = 10
    paramDict['unimer'] = unimer
    paramDict['particle'] = particle
      
    xyz, names = MelaninParticle(paramDict)
    fIO.saveXYZList(xyz, names, "melanin.xyz")

    print("example done")