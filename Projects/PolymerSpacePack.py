import numpy as np
import random as rnd
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
import Utilities.coordSystems as coords 
import Utilities.fileIO as fIO

def PolymerSpacePack(polymerParams):

    # create a backbone generator object for the polymer
    backboneBBG = RPPBBG(polymerParams['inpFilename'])
            
    # generate the envelope containing the polymer
    envelopeList = ['frustum ' + str(polymerParams['Z1']) + ' ' + str(polymerParams['R1']) + ' ' + str(polymerParams['Z2']) + ' ' + str(polymerParams['R2'])]
    
    try:
        pointsToAvoid = polymerParams['PointsToAvoid']
    except KeyError:
        pointsToAvoid = []
    
    # generate polymerBB of appropriate length. 
    backboneBB = backboneBBG.generateBuildingBlock( polymerParams['numMonomers'],  
                                                    polymerParams['startPoint'],
                                                    polymerParams['alpha1'],
                                                    polymerParams['alpha2'],
                                                    polymerParams['beta1'],
                                                    polymerParams['beta2'],
                                                    polymerParams['minDist'],
                                                    polymerParams['bondLength'],
                                                    envelopeList=envelopeList,
                                                    visualiseEnvelope=(polymerParams['numEnv'], polymerParams['envScale'] , polymerParams['outFileroot'] + '_envelope.xyz'),
                                                    pointsToAvoid=pointsToAvoid)
        
    fIO.saveXYZ(backboneBB.blockXYZVals, "C", polymerParams['outFileroot']+'.xyz')

    
if __name__=="__main__":

    # packs a polymer made of a certain number of spheres and bondlength 
    # packed into different overall cylinders all of the same volume. 
    # Computes key geometric parameters for each polymer and plots the results in a graph against the input range of values
    import os
    import glob
    files = glob.glob("*.xyz")
    for f in files:
        os.remove(f)
    
    NumberOfParticlesInChain = 1000
    bondLength = 2.0
    particleRadius = 1.0
    particleVolume = (4 / 3) * np.pi * particleRadius**3 
    minDist = 2.0
    EmptySpacePerParticle = 100 * particleVolume 
    CylVolume =  NumberOfParticlesInChain *( EmptySpacePerParticle + particleVolume )  
    CylRadii = [10.0, 30.0, 50.0, 100.0, 200.0]
    Repeats = 1
    
    numZeroes = int(np.ceil(np.log10(Repeats + 1)))
    print("numZeroes:", numZeroes)
    
    # cylinder volume is v = np.pi * r^2 * h, so h = v/np.pi * r^2
    # generate the necessary cylinder height for each specified radius 
    CylZs = [ CylVolume / ( np.pi  * (r**2) )  for r in CylRadii ]
    
    # dump the input values
    [ print("Cylinders: r: " + str(r) + ", z:" + str(z) + ", v_act: " + str(z * np.pi * (r**2)) + ", v_target: " + str(CylVolume) ) for r, z in zip(CylRadii, CylZs) ]
    

    # loop through each packing cylinder and fit N repeats of each polymer. Save each polymer in an Xyz file for the cylinder parameters.
    for r, z in zip(CylRadii, CylZs):
        for i in range(Repeats):
            # set up the dictionary
            polymer = {}
            
            
            
            # pack the dictionaries
            polymer['inpFilename'] = 'RandomPolymer.txt'  
            polymer['outFileroot'] = 'Polymer_' + str(r) + '_' + str(z) + '_' + str(CylVolume) + '_' + str(i).zfill(numZeroes)
            polymer['numMonomers'] = NumberOfParticlesInChain 
            polymer['alpha1'] = -50.0
            polymer['alpha2'] = 50.0
            polymer['beta1'] = 100.0
            polymer['beta2'] = 120.0
            polymer['minDist'] = 1.0
            polymer['bondLength'] = 1.0
            polymer['Z1'] = 0
            polymer['R1'] = r
            polymer['Z2'] = z 
            polymer['R2'] = r
            polymer['envScale'] = max(2.0 * r, z)
            polymer['numEnv'] = min(1000000, int(np.floor(100000*polymer['envScale']**3/CylVolume)))
            
            # pick a random point in a unit cylinder (r=1, h=1) 
            (r0, phi0 ,z0) = coords.pickRandomPointInUnitCylinder()
            
            # scale point to size of current cylinder - make sure starting point is in the first 10% of the cylinder
            polymer['startPoint'] = coords.cyl2XYZ(np.array([r * r0, phi0, 0.1 * z * z0])) 
            print(r0, phi0, z0, polymer['startPoint'])
        
        
            PolymerSpacePack(polymer)
       
    print("done")