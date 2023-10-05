import numpy as np
from Library.VolumePackCylinder import VolumePackCylinderBBG as VPCBBG
import Utilities.fileIO as fIO
import itertools as iter
import matplotlib.pyplot as plt

# pinched from pdblib.py
def RadiusOfGyration(vectors):
    arr = np.array(vectors)
    center_of_mass = np.mean(arr, axis=0)
    distances = np.sqrt(np.sum((arr - center_of_mass)**2, axis=1))
    return np.sqrt(np.sum(distances**2) / len(vectors))

# returns the two points that are furthest away from each other in an input list of np arrays.
def maxDistInList(listofpoints):
    dist_array = np.array([ [p, q, np.linalg.norm(listofpoints[p] - listofpoints[q]) ] for p, q in iter.combinations(range(len(listofpoints)), 2) ]) 
    maxIndices = tuple(dist_array[ np.argmax( dist_array[:,2] )][0:2])   
    return listofpoints[int(maxIndices[0])], listofpoints[int(maxIndices[1])] 

    
if __name__=="__main__":

    # Packs a range of cylindrical regions with the same volume with spheres. The number density of all the regions is the same.  
    
    # Plots RG, maxEuclid, anistropy (rg/max) and Rg/radius for a range of radii for the resulting distributions.    


    # remove all previous xyz file     
    import os
    import glob
    files = glob.glob("*.xyz")
    for f in files:
        os.remove(f)
    
    NumberOfParticles = 3000
    CylVolume = (np.pi * (2**2)) * 50  # typical volume of 2nm radius tube with 50 nm length.
    particleVolume = CylVolume/NumberOfParticles
    minDist = (3.0 * particleVolume/ (4.0 * np.pi))**(1.0/3.0)
    CylRadii = [0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0, 50.0, 75.0, 100]
    Repeats = 3
 
    with open('data.out',"w") as file:
        file.write("CylVolume: " + f'{CylVolume:8.2f}' + ", minDist:" + f'{minDist:8.2f}' + '\n')
    print("CylVolume: ", f'{CylVolume:8.2f}', ", minDist:", f'{minDist:8.2f}')
    
    numZeroes = int(np.ceil(np.log10(Repeats + 1)))
    print("numZeroes:", numZeroes)
    
    # cylinder volume is v = np.pi * r^2 * h, so h = v/np.pi * r^2
    # generate the necessary cylinder height for each specified radius 
    CylZs = [ CylVolume / ( np.pi  * (r**2) )  for r in CylRadii ]
    
    # dump the input values
    [ print("Cylinders: r: " + f'{r:8.1f}' + ", z:" + f'{z:8.2f}' + ", v_act: " + f'{z * np.pi * (r**2):8.2f}' + ", v_target: " + f'{CylVolume:8.2f}' ) for r, z in zip(CylRadii, CylZs) ]
    with open('data.out',"a") as file:
        [ file.write("r: " + f'{r:8.1f}' + ", z:" + f'{z:8.2f}' + ", v_act: " + f'{z * np.pi * (r**2):8.2f}' + ", v_target: " + f'{CylVolume:8.2f}' + '\n') for r, z in zip(CylRadii, CylZs) ]
 
    rs = []
    zs = []
    Rgs = []
    MaxEs = []
    Ans = []
    RgRad = []
    densities = []

    with open('data.out',"a") as file:
        file.write(" r, z, n, v, Rg, MaxEuc, Rg/MaxE, Rg/r, r/z \n")
            
    # loop through each packing cylinder and perform N repeats of each pack. Save each dump in an Xyz file for the cylinder parameters.
    for r, z in zip(CylRadii, CylZs):
        for i in range(Repeats):

            
            # create the NPack object.
            CylinderPackBBG = VPCBBG('VolumePackCylinder.txt')
       
            # generate the building block
            CylPackBB = CylinderPackBBG.generateBuildingBlock(NumberOfParticles, r, z, minDist)
 
            # extract xyz Vals
            xyzVals = CylPackBB.blockXYZVals
    
            outFilename = 'VPC_' + f'{r:.1f}' + '_' + f'{z:.0f}' + '_' + f'{CylVolume:.0f}' + '_' + str(i).zfill(numZeroes) + '.xyz'
    
            # dump to file
            fIO.saveXYZ(CylPackBB.blockXYZVals, "C", outFilename)


            print("Analysing Structure")
            Rg = RadiusOfGyration(xyzVals)
            MaxPair = maxDistInList(xyzVals)
            MaxEuclid = np.linalg.norm(MaxPair[0] - MaxPair[1])
            
            print(f'{r:6.1f}', ' ', f'{z:6.2f}', ' ', str(NumberOfParticles), ' ', f'{CylVolume:6.2f}', ' ', f'{Rg:6.2f}', ' ', f'{MaxEuclid:6.2f}', ' ', f'{Rg/MaxEuclid:9.2f}', ' ', f'{Rg/r:9.2f}', ' ', f'{z/r:9.2f}')
            with open('data.out',"a") as file:
                file.write(f'{r:6.1f}' + ' ' + f'{z:6.2f}' + ' ' + str(NumberOfParticles) + ' ' + f'{CylVolume:6.2f}' + ' ' + f'{Rg:6.2f}' + ' ' + f'{MaxEuclid:6.2f}' + ' ' + f'{Rg/MaxEuclid:9.2f}' + ' ' + f'{Rg/r:9.2f}'+ ' ' + f'{z/r:9.2f}' + '\n')
            
            rs.append(r)
            zs.append(z)
            Rgs.append(Rg)
            MaxEs.append(MaxEuclid)
            Ans.append(Rg/MaxEuclid)
            RgRad.append(Rg/r)
            densities.append(NumberOfParticles/CylVolume)
            

    rs = np.average(np.array(rs).reshape((int(np.floor(len(rs)/3)),3)),axis=1)
    zs = np.average(np.array(zs).reshape((int(np.floor(len(zs)/3)),3)),axis=1)
    Rgs = np.average(np.array(Rgs).reshape((int(np.floor(len(Rgs)/3)),3)),axis=1)
    MaxEs = np.average(np.array(MaxEs).reshape((int(np.floor(len(MaxEs)/3)),3)),axis=1)
    Ans = np.average(np.array(Ans).reshape((int(np.floor(len(Ans)/3)),3)),axis=1)
    RgRad = np.average(np.array(RgRad).reshape((int(np.floor(len(RgRad)/3)),3)),axis=1)
    densities = np.average(np.array(densities).reshape((int(np.floor(len(densities)/3)),3)),axis=1)

    #set the current figure
    fig = plt.figure()
    plt.plot( rs, Rgs, label='Rg vs R', 
              markeredgewidth=1,
              markeredgecolor="#ff0000",
              markerfacecolor="#ff0000",
              marker='+',
              linestyle='',
              color="#ff0000",
              zorder=1)
    plt.legend(loc='upper right')
    plt.xlabel('r (a.u.)')
    plt.ylabel('Rg (a.u.)')
    plt.savefig('RgVsR.png')
    plt.show()

    fig = plt.figure()
    plt.plot( rs, MaxEs, label='MaxEuclid vs R', 
              markeredgewidth=1,
              markeredgecolor="#ff0000",
              markerfacecolor="#ff0000",
              marker='+',
              linestyle='',
              color="#ff0000",
              zorder=1)
    plt.legend(loc='upper right')
    plt.xlabel('r (a.u.)')
    plt.ylabel('MaxEuclid (a.u.)')
    plt.savefig('MaxEucVsR.png')
    plt.show()


    fig = plt.figure()
    plt.plot( rs, Ans, label='Anisotropy vs R', 
              markeredgewidth=1,
              markeredgecolor="#ff0000",
              markerfacecolor="#ff0000",
              marker='+',
              linestyle='',
              color="#ff0000",
              zorder=1)
    plt.legend(loc='lower right')
    plt.xlabel('r (a.u.)')
    plt.ylabel('Anisotropy (a.u.)')
    plt.savefig('AnisotropyVsR.png')
    plt.show()

    fig = plt.figure()
    plt.plot( rs, RgRad, label='Rg/r vs R', 
              markeredgewidth=1,
              markeredgecolor="#ff0000",
              markerfacecolor="#ff0000",
              marker='+',
              linestyle='',
              color="#ff0000",
              zorder=1)
    plt.legend(loc='upper right')
    plt.xlabel('r (a.u.)')
    plt.ylabel('Rg/r (a.u.)')
    plt.savefig('RgoverR_VsR.png')
    plt.show()

    fig = plt.figure()
    plt.plot( rs, densities, label='Density vs R', 
              markeredgewidth=1,
              markeredgecolor="#ff0000",
              markerfacecolor="#ff0000",
              marker='+',
              linestyle='',
              color="#ff0000",
              zorder=1)
    plt.legend(loc='upper right')
    plt.xlabel('r (a.u.)')
    plt.ylabel('Density (a.u.)')
    plt.savefig('DensityVsR.png')
    plt.show()

    fig = plt.figure()
    plt.plot( rs, [z/r for r,z in zip(rs,zs)], label='Z vs R', 
              markeredgewidth=1,
              markeredgecolor="#ff0000",
              markerfacecolor="#ff0000",
              marker='+',
              linestyle='',
              color="#ff0000",
              zorder=1)
    plt.legend(loc='upper right')
    plt.xlabel('r (a.u.)')
    plt.ylabel('Z/r (a.u.)')
    plt.savefig('ZOverRVsR.png')
    plt.show()



    print("done")