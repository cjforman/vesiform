import numpy as np
from Library.randomPolymer import RandomPolymerPackBBG as RPBBG
import Utilities.fileIO as fIO


def makeBlockCopolymer(numA, numB, Z1_A, R1_A, Z2_A, R2_A, alpha1A, alpha2A, beta1A, beta2A, minDistA, bondLengthA, 
                                   Z1_B, R1_B, Z2_B, R2_B, alpha1B, alpha2B, beta1B, beta2B, minDistB, bondLengthB, filename):
    
    # create the generator
    RandomPolymerPackBBG = RPBBG(filename)
    
    # Z1 = envelopeParams[0]
    # Radius1 = envelopeParams[1]
    # Z2 = envelopeParams[2]
    # Radius2 = envelopeParams[3]
    
    pointA = np.array([ 0.0, 0.0, Z1_A])
    envelopeList = ['frustum ' + str(Z1_A) + ' ' + str(R1_A) + ' ' + str(Z2_A) + ' ' + str(R2_A)]
    
    
    rejectPolymer = True
    
    while rejectPolymer:
        # generate polymer A of appropriate length. 
        RandomPolymer_A = RandomPolymerPackBBG.generateBuildingBlock( numA + 1, 
                                                                      pointA,
                                                                      alpha1A,
                                                                      alpha2A,
                                                                      beta1A,
                                                                      beta2A,
                                                                      minDistA, 
                                                                      bondLengthA,
                                                                      envelopeList=envelopeList,
                                                                      visualiseEnvelope=(0, 100, 'blockCopolymer_envelope.xyz'))
        # check that Z coord of point 0 is below z coord of point numA - 1  otherwise reject polymer and choose another one.  
        if RandomPolymer_A.blockXYZVals[0][2] < RandomPolymer_A.blockXYZVals[-1][2]:
            rejectPolymer = False
        else:
            print "Polymer Rejected"



    pointB = np.array([ 0.0, 0.0, Z1_B])
    envelopeList = ['frustum ' + str(Z1_B) + ' ' + str(R1_B) + ' ' + str(Z2_B) + ' ' + str(R2_B)]
    
    # generate polymer B of appropriate length. 
    RandomPolymer_B = RandomPolymerPackBBG.generateBuildingBlock( numB, 
                                                                  pointB,
                                                                  alpha1B,
                                                                  alpha2B,
                                                                  beta1B,
                                                                  beta2B,
                                                                  minDistB, 
                                                                  bondLengthB,
                                                                  envelopeList=envelopeList,
                                                                  visualiseEnvelope=(0, 100, 'blockCopolymer_envelope.xyz'))
    
    xyzVals = RandomPolymer_A.blockXYZVals[0:-1]
    xyzVals = np.concatenate((xyzVals, [ a + RandomPolymer_A.blockXYZVals[-1] - RandomPolymer_B.blockXYZVals[0] for a in RandomPolymer_B.blockXYZVals] ), 0)
    names =  ['O'] * numA 
    names = np.concatenate( (names, ['P'] * numB), 0 ) 

    return (xyzVals, names)

if __name__=="__main__":

    filename = 'RandomPolymer.txt' 
    
    bondLengthA = 1.5
    minDistA = 1.0
    numA = 5
    Z1_A = 2
    R1_A = 20
    Z2_A = 500
    R2_A = 5
    alpha1A= 40
    alpha2A = 60
    beta1A = 165
    beta2A = 185
    
    numB = 40
    bondLengthB = 1.5
    minDistB = 1.0
    Z1_B = 2
    R1_B = 20
    Z2_B = 500
    R2_B = 5
    alpha1B= 40
    alpha2B = 50
    beta1B = 145
    beta2B = 155

    numPolymers = 20
    
    for i in range(numPolymers):
        xyzVals, names = makeBlockCopolymer(numA, numB, Z1_A, R1_A, Z2_A, R2_A, alpha1A, alpha2A, beta1A, beta2A, minDistA, bondLengthA, 
                                                        Z1_B, R1_B, Z2_B, R2_B, alpha1B, alpha2B, beta1B, beta2B, minDistB, bondLengthB, filename)

        fIO.saveXYZList(xyzVals, names, "blockCopolymer_" + str(i) + ".xyz")


    print "example done"