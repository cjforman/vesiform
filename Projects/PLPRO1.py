'''
Created on Aug 17, 2021

@author: Chris

Project classes for the Tau/PLP RO1 proposal.

'''
from Projects.GeneralPolymerBrush import GeneralBlockPolymer as GBP
from Library.SurfacePackSphere import SurfacePackSphereBBG as SPSBBG
from Library.SurfacePackPlane import SurfacePackPlaneBBG as SPPBBG
import Utilities.fileIO as fIO
import Utilities.coordSystems as coords
import Utilities.cartesian as carts
from PDBProc.pdbLib import PDB
import numpy as np
import copy as cp

class UnfoldedTau():
    def __init__(self):
        # set up the dictionaries
        backboneDict1 = {}
        brushDict1 = {}
        
        # pack the dictionaries
        backboneDict1['filename'] = 'RandomPolymer.txt'  
        backboneDict1['mode'] = 'Polymer'
        backboneDict1['name'] = 'BlockA'
        backboneDict1['numMonomers'] = 441 
        backboneDict1['monomerNames'] = 441 * ['H'] 
        # backboneDict1['monomerNames'] += (74 - 45 + 1) * ['He'] 
        # backboneDict1['monomerNames'] += (103 - 75 + 1) * ['Li'] 
        # backboneDict1['monomerNames'] += (150 - 104 + 1) * ['Be'] 
        # backboneDict1['monomerNames'] += (198 - 151 + 1) * ['B'] 
        # backboneDict1['monomerNames'] += (243 - 199 + 1) * ['C'] 
        # backboneDict1['monomerNames'] += (274 - 244 + 1) * ['N']
        # backboneDict1['monomerNames'] += (305 - 275 + 1) * ['O'] 
        # backboneDict1['monomerNames'] += (336 - 306 + 1) * ['F'] 
        # backboneDict1['monomerNames'] += (368 - 337 + 1) * ['Ne'] 
        # backboneDict1['monomerNames'] += (441 - 369 + 1) * ['Na']

        print(len(backboneDict1['monomerNames']))
                                                 
        backboneDict1['monomerNames']
        backboneDict1['alpha1'] = -20
        backboneDict1['alpha2'] = 20
        backboneDict1['beta1'] = 165
        backboneDict1['beta2'] = 185
        backboneDict1['minDist'] = 2.1
        backboneDict1['bondLength'] = 2.0
        backboneDict1['Z1'] = 0
        backboneDict1['R1'] = 200
        backboneDict1['Z2'] = 0.5 * 1.5 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
        backboneDict1['R2'] = 200
    
        brushDict1['filename'] = 'RandomPolymer.txt'  
        brushDict1['mode'] = 'PolymerAlternate'
        brushDict1['name'] ='brush1'
        brushDict1['numMonomers'] = 0
        brushDict1['monomerNames'] = []
        brushDict1['indexOfFirstBrush'] = 1
        brushDict1['spacing'] = 1
        brushDict1['phaseRange'] = 15.0 
        brushDict1['alpha1'] = 40
        brushDict1['alpha2'] = 50
        brushDict1['beta1'] = 145
        brushDict1['beta2'] = 155
        brushDict1['minDist'] = 1.0
        brushDict1['bondLength'] = 2.0
        brushDict1['Z1'] = 2.0
        brushDict1['R1'] = 10.0
        brushDict1['Z2'] = brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
        brushDict1['R2'] = 30.0
    
        polymerBrushDict={}
        polymerBrushDict['backbones'] = [backboneDict1]#, backboneDict1, backboneDict2]
        polymerBrushDict['brushes'] = [brushDict1]#, brushDict1, brushDict2]
        polymerBrushDict['connectors'] = []#, connectorDict12, connectorDict12]
    
        xyz, names = GBP(polymerBrushDict)
        fIO.saveXYZList(xyz, names, "TauUnfolded.xyz")

class TauFolded():
    def __init__(self):
        pass

class ProtoFilamentCluster():
    def __init__(self, numPreFibrils, clusterRange, clusterIndex, alphaMin=-20, alphaMax=20, betaMin=165, betaMax=185, bondLength=2.0, minDist=2.1, clusterMode='loose', sphericalParamsDict={}):
        if clusterMode=='spherical':
            sphereBBG = SPSBBG('SurfacePackSphere.txt')
            # generate the sphere base Points centered on origin by default
            SphereBB = sphereBBG.generateBuildingBlock(numPreFibrils, 
                                                       sphericalParamsDict['radius'], 
                                                       sphericalParamsDict['theta1'], 
                                                       sphericalParamsDict['theta2'],
                                                       sphericalParamsDict['phi1'],
                                                       sphericalParamsDict['phi2'],
                                                       sphericalParamsDict['minDist'])
            SpherePoints = SphereBB.blockXYZVals
            
            # computer directors
            directors = [spoint/np.linalg.norm(spoint) for spoint in SpherePoints]
        
        # generate numPreFibril filaments. 
        for fibNum in range(0, numPreFibrils):
            # compute the new filament
            newXYZVals, newNames = self.unfoldedTau(alphaMin=alphaMin, 
                                                    alphaMax=alphaMax,
                                                    betaMin=betaMin, 
                                                    betaMax=betaMax,
                                                    bondLength=bondLength,
                                                    minDist=minDist)
            
            if clusterMode=='loose':
                # pick an index around the index that we want the filament to cluster loosely.
                randomIndex = np.random.choice(np.arange(clusterIndex - 5, clusterIndex + 5))
    
                # add a random vector to the filament so they don't intersect perfectly            
                randomVector = np.random.uniform(low=-clusterRange, high=clusterRange, size=3)
                newXYZVals = newXYZVals - newXYZVals[randomIndex] + randomVector
            
            if clusterMode=='spherical':
                newXYZVals = coords.transformFromBlockFrameToLabFrame( directors[fibNum], 
                                                                       SpherePoints[fibNum],
                                                                       0.0,
                                                                       np.array([0.0, 0.0, 1.0]),
                                                                       newXYZVals[-1],
                                                                       newXYZVals)
            if fibNum==0:
                xyzVals = cp.copy(newXYZVals)
                names = cp.copy(newNames)
            else:
                # add the newVals to the ever growing output vals 
                xyzVals= np.concatenate( (xyzVals, newXYZVals ), 0 )
                names = np.concatenate(  (names, newNames), 0 )
            
        fIO.saveXYZList(xyzVals, names, "TauUnfoldedCluster.xyz")
            

class ProtoFilamentClusterAligned():
    
    def __init__(self, numPreFibrils, clusterIndex, offsetRange, r1=10, bondLength=2.0, minDist=2.1):
      
        backBoneLengths = [29, 12, 7]
        backbones = ['H', 'He', 'Li']
        sequences = [[],[],[]]
        BBAlphas = [(50, 80), (45, 50), (50, 80)]
        BBBetas =[ (140, 170), (165, 175), (140, 170)]
      
        xyzVals = []
      
        # generate numPreFibril filaments. 
        for fibNum in range(0, numPreFibrils):

            PLPName = 'WompRat_' + str(fibNum)
            
            # compute the new filament
            GPLP = GenericPLP(backBoneLengths, backbones, sequences, PLPName, BBAlphas=BBAlphas, BBBetas=BBBetas)
            newXYZVals, oldNames = GPLP.getVals() 
            newNames = 5 * ['H'] 
            newNames += 3 * ['He'] 
            newNames += 3 * ['Li'] 
            newNames += 6 * ['Be'] 
            newNames += 6 * ['B'] 
            newNames += 6 * ['C'] 
            newNames += 3 * ['N']
            newNames += 3 * ['O'] 
            newNames += 3 * ['F'] 
            newNames += 3 * ['Ne'] 
            newNames += 7 * ['Na']  
            
            offsetVector = np.array([fibNum * minDist * offsetRange, 0, 0])
            newXYZVals = newXYZVals - newXYZVals[clusterIndex] + offsetVector 
      
            if fibNum==0:
                xyzVals = cp.copy(newXYZVals)
                names = cp.copy(newNames)
            else:
                # add the newVals to the ever growing output vals 
                xyzVals= np.concatenate( (xyzVals, newXYZVals ), 0 )
                names = np.concatenate(  (names, newNames), 0 )
            
        fIO.saveXYZList(xyzVals, names, "AlignedCluster.xyz")
        
    def unfoldedTau(self, alphaMin=-20, alphaMax=20, betaMin=165, betaMax=185, bondLength=2.0, minDist=2.1, r1=2):
        # set up the dictionaries
        backboneDict1 = {}
        brushDict1 = {}
        
        # pack the dictionaries
        backboneDict1['filename'] = 'RandomPolymer.txt'  
        backboneDict1['mode'] = 'Polymer'
        backboneDict1['name'] = 'BlockA'
        backboneDict1['monomerNames'] = 5 * ['H'] 
        backboneDict1['monomerNames'] += 3 * ['He'] 
        backboneDict1['monomerNames'] += 3 * ['Li'] 
        backboneDict1['monomerNames'] += 6 * ['Be'] 
        backboneDict1['monomerNames'] += 6 * ['B'] 
        backboneDict1['monomerNames'] += 6 * ['C'] 
        backboneDict1['monomerNames'] += 3 * ['N']
        backboneDict1['monomerNames'] += 3 * ['O'] 
        backboneDict1['monomerNames'] += 3 * ['F'] 
        backboneDict1['monomerNames'] += 3 * ['Ne'] 
        backboneDict1['monomerNames'] += 7 * ['Na']
        backboneDict1['numMonomers'] = len(backboneDict1['monomerNames']) 
        backboneDict1['alpha1'] = alphaMin
        backboneDict1['alpha2'] = alphaMax
        backboneDict1['beta1'] = betaMin
        backboneDict1['beta2'] = betaMax
        backboneDict1['minDist'] = minDist
        backboneDict1['bondLength'] = bondLength
        backboneDict1['Z1'] = bondLength
        backboneDict1['R1'] = 2
        backboneDict1['Z2'] = 0.5 * 1.5 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
        backboneDict1['R2'] = 200
    
        brushDict1['filename'] = 'RandomPolymer.txt'  
        brushDict1['mode'] = 'PolymerAlternate'
        brushDict1['name'] ='brush1'
        brushDict1['numMonomers'] = 0
        brushDict1['monomerNames'] = []
        brushDict1['indexOfFirstBrush'] = 1
        brushDict1['spacing'] = 1
        brushDict1['phaseRange'] = 15.0 
        brushDict1['alpha1'] = 40
        brushDict1['alpha2'] = 50
        brushDict1['beta1'] = 145
        brushDict1['beta2'] = 155
        brushDict1['minDist'] = minDist
        brushDict1['bondLength'] = bondLength
        brushDict1['Z1'] = bondLength
        brushDict1['R1'] = 10.0
        brushDict1['Z2'] = brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
        brushDict1['R2'] = 30.0
    
        polymerBrushDict={}
        polymerBrushDict['backbones'] = [backboneDict1]#, backboneDict1, backboneDict2]
        polymerBrushDict['brushes'] = [brushDict1]#, brushDict1, brushDict2]
        polymerBrushDict['connectors'] = []#, connectorDict12, connectorDict12]
    
        return GBP(polymerBrushDict)
        
# this function finds the rms distance between the overlapping atoms of two adjacent helical bits
def computeOffset(x, xyz, names, axis, COM):
    
    # generate an amyloid with the given offset
    outBlock, nameBlock = makeAmyloid(xyz, names, 2, x[0], x[1], axis, COM)
 
    minBr = 0
    maxBr = minBr + 1136 * 2
    minTr = 20448  
    maxTr = minTr + 1136 * 2

    rms = sum( np.linalg.norm( outBlock[minBr:maxBr,:] - outBlock[minTr:maxTr, :], axis=1) )

    print(rms, x[0]/4, x[1]/4)
    # rmsd distance between equivalent points           
    return rms  

# function creates a new block of NumCopies of the XYZ and Names by rotating 
# XYZ input by angle about an axis through P0 and then offsetting by vOffset
def makeAmyloid(XYZ, Names, numCopies, angle, vOffset, axis, P0):
    # seed the output with a copy of the input at the blockNum = 1 level
    outBlock = [ carts.rotPAboutAxis(p - P0, axis, -1 * angle * np.pi/180) + P0 + vOffset * axis for p in XYZ ]
    nameBlock = cp.copy(Names)
    
    # create N-1 copies of the block of atoms, 
    # rotates each point in the block by an amount 'angle' about an axis at P0
    # translates each rotated point by an amount vOffset along the axis
    for blockNum in range(2, numCopies):
        print(blockNum, len(outBlock))
        newBlock = [ carts.rotPAboutAxis(p - P0, axis, -blockNum * angle * np.pi/180) + P0 + blockNum * vOffset * axis for p in XYZ ]
        
        outBlock = np.concatenate( (outBlock, newBlock), 0 )
        nameBlock = np.concatenate( (nameBlock, cp.copy(Names)), 0 )
    
    return outBlock, nameBlock

class TauAmyloid():
    def __init__(self, filename, numCopies, angle, vOffset):
        # get the PDB object
        PDBObj = PDB(filename)

        # extract atoms from the PDB
        atomsXYZ = np.asarray( PDBObj.extractCoords() )

        # extract the names from the PDB
        atomNames = np.asarray([PDBObj.getFieldFromPDB(1)])[0]

        # replace CAs with Cs
        atomNames = np.asarray([an if not an=='CA' else 'C' for an in atomNames ])


        # output the xyz file of the raw PDB file
        fIO.saveXYZList(atomsXYZ, atomNames, filename + '.xyz')

        # compute centre of mass of object     
        COM = carts.getCentreOfMass(atomsXYZ)
        
        #get principal axis
        # axis  = carts.getPrincipalAxis(atomsXYZ)
        #axis  = coords.axisFromHelix(atomsXYZ)
        
        # axis = np.array([0.0, 0.0, 1.0])
        # axis = np.array([0.00590127, -0.21029529, 0.9776201]) #obtained from a large amyloid Mostly by hacking. 
        axis = np.array([0.02159924, -0.18547132,  0.98241227])
        
        # normalize axis
        axis = axis/np.linalg.norm(axis)
        print(axis)

        # res = opt.minimize(computeOffset, 
        #                   [angle, vOffset], 
        #                   args = (atomsXYZ, atomNames, axis, COM), 
        #                   method='nelder-mead', 
        #                    )        

        # FineAngle = res.x[0]
        # FineOffset = res.x[1]
        # print("FineAngle:", FineAngle/4)
        # print("FineOffset:", FineOffset/4)
       
        # extract backboneIndices
        bbIndices = np.asarray( PDBObj.extractBackBoneIndices() )
 
        # get the backbone names 
        bbNames = atomNames[bbIndices]

        bbNames = [bb if not bb=='CA' else 'C' for bb in bbNames ]
        
        # get the backboneXYZ coords
        bbXYZ = atomsXYZ[bbIndices] 

        caIndices = np.asarray( PDBObj.extractCas() )
        caXYZ = atomsXYZ[caIndices]
        caNames = len(caXYZ) * ['C']

        # compute new amyloid
        outBlockBB, nameBlockBB = makeAmyloid(bbXYZ, bbNames, numCopies, angle, vOffset, axis, COM)
        outBlockCA, nameBlockCA = makeAmyloid(caXYZ, caNames, numCopies, angle, vOffset, axis, COM)
        outBlockAll, nameBlockAll = makeAmyloid(atomsXYZ, atomNames, numCopies, angle, vOffset, axis, COM)

        # compute a axis based on the large amyloid       
        # newAxis  = coords.axisFromHelix(outBlockAll)
       
        # print(newAxis)
       
        # save the amyloid        
        fIO.saveXYZList(outBlockAll, nameBlockAll, "AmyloidAll.xyz")
        fIO.saveXYZList(outBlockCA, nameBlockCA, "AmyloidCA.xyz")
        fIO.saveXYZList(outBlockBB, nameBlockBB, "AmyloidBB.xyz")

class GenericPLP():
    def __init__(self, backboneLengths, backbones, sequences, PLPName, doProlineKink=False, BBAlphas=[(40,50)], BBBetas=[(165,185)], BrAlphas=[(40,50)], BrBetas=[(145,155)]   ):

        polymerBrushDict={}
        polymerBrushDict['backbones'] = []
        polymerBrushDict['brushes'] = []
        polymerBrushDict['connectors'] = []
        
        if len(BBAlphas)==1:
            BBAlphas = len(backbones) * BBAlphas
        if len(BBBetas)==1:
            BBBetas = len(backbones) * BBBetas
        if len(BrAlphas)==1:
            BrAlphas = len(backbones) * BrAlphas
        if len(BrBetas)==1:
            BrBetas = len(backbones) * BrBetas
        
        for sequence, backbone, length, alphaBBRange, betaBBRange, alphaBrRange, betaBrRange in zip(sequences, backbones, backboneLengths, BBAlphas, BBBetas, BrAlphas, BrBetas):
            
            # set up the dictionaries
            backboneDict = {}
            brushDict = {}

            # pack the dictionaries
            backboneDict['filename'] = 'RandomPolymer.txt'  
            backboneDict['mode'] = 'Polymer'
            backboneDict['name'] = 'BlockA'
            backboneDict['numMonomers'] = length 
            backboneDict['monomerNames'] = length * [backbone]
            backboneDict['alpha1'] = alphaBBRange[0]
            backboneDict['alpha2'] = alphaBBRange[1]
            backboneDict['beta1'] = betaBBRange[0]
            backboneDict['beta2'] = betaBBRange[1]
            backboneDict['minDist'] = 2.0
            backboneDict['bondLength'] = 2.0
            backboneDict['Z1'] = 2
            backboneDict['R1'] = 20
            backboneDict['Z2'] = 1.5 *backboneDict['numMonomers'] * backboneDict['bondLength'] + backboneDict['Z1']
            backboneDict['R2'] = 70
        
            brushDict['filename'] = 'RandomPolymer.txt'  
            brushDict['mode'] = 'PolymerAlternate'
            brushDict['name'] ='brush'
            brushDict['numMonomers'] = len(sequence)         
            brushDict['monomerNames'] = [ self.convertPeptideToElementsForBlender(peptide) for peptide in sequence ] # if specified in peptide mode then overrides the C-CA-N peptide names
            brushDict['indexOfFirstBrush'] = 0
            brushDict['spacing'] = 1
            brushDict['phaseRange'] = 10.0 
            brushDict['alpha1'] = alphaBrRange[0]
            brushDict['alpha2'] = alphaBrRange[1]
            brushDict['beta1'] = betaBrRange[0]
            brushDict['beta2'] = betaBrRange[1]
            brushDict['minDist'] = 1.0
            brushDict['bondLength'] = 2.0
            brushDict['Z1'] = 2
            brushDict['R1'] = 50
            brushDict['Z2'] = 1.5 * brushDict['numMonomers'] * brushDict['bondLength'] + brushDict['Z1']
            brushDict['R2'] = 50
            brushDict['doProlineKink']=doProlineKink
            brushDict['prolineKinkAngle']=90 * np.pi/180

            polymerBrushDict['backbones'].append(cp.copy(backboneDict))
            polymerBrushDict['brushes'].append(cp.copy(brushDict))

        # set up any connectors that we need (loop through the first N-1 backbones
        for backbone in (polymerBrushDict['backbones'])[0:-1]:
            connectorDict = {}
            connectorDict['displacement'] = backbone['bondLength']
            connectorDict['alpha'] = 0.0
            connectorDict['beta'] = 170.0

            polymerBrushDict['connectors'].append(cp.copy(connectorDict))
        
        self.xyz, self.names = GBP(polymerBrushDict)
        fIO.saveXYZList(self.xyz, self.names, PLPName + ".xyz")

    def getVals(self):
        return self.xyz, self.names

    def convertPeptideToElementsForBlender(self, peptide):
        peptides = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'Fl', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']
        return elements[ peptides.index(peptide) ]

class CappedGrowthPLPProtoFilamentCluster():
    def __init__(self):
        pass

class SequestrationPLP():
    def __init__(self, Randomize=False, save=True, scale=1.0):
        # set up the dictionaries
        backboneDict1 = {}
        backboneDict2 = {}
        connectorDict12 = {}
        brushDict1 = {}
        brushDict2 = {}
        
        # pack the dictionaries
        backboneDict1['filename'] = 'RandomPolymer.txt'  
        backboneDict1['mode'] = 'Polymer'
        backboneDict1['name'] = 'BlockA'
        backboneDict1['numMonomers'] = 15 
        backboneDict1['monomerNames'] = backboneDict1['numMonomers'] * ['Pd']
        backboneDict1['alpha1'] = 40
        backboneDict1['alpha2'] = 50
        backboneDict1['beta1'] = 165
        backboneDict1['beta2'] = 185
        backboneDict1['minDist'] = scale * 2.0
        backboneDict1['bondLength'] = scale * 2.0
        backboneDict1['Z1'] = scale * 2
        backboneDict1['R1'] = scale * 20
        backboneDict1['Z2'] = scale * 1.5 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
        backboneDict1['R2'] = scale * 70
    
        backboneDict2['filename'] = 'RandomPolymer.txt'  
        backboneDict2['mode'] = 'Polymer'
        backboneDict2['name'] = 'BlockB'
        backboneDict2['numMonomers'] = 15
        backboneDict2['monomerNames'] = backboneDict2['numMonomers'] *['Pb']
        backboneDict2['indexOfFirstBrush'] = 0
        backboneDict2['spacing'] = 0        
        backboneDict2['alpha1'] = 40
        backboneDict2['alpha2'] = 50
        backboneDict2['beta1'] = 165
        backboneDict2['beta2'] = 185
        backboneDict2['minDist'] = scale * 2.0
        backboneDict2['bondLength'] = scale * 2.0
        backboneDict2['Z1'] = scale * 2
        backboneDict2['R1'] = scale * 20
        backboneDict2['Z2'] = scale * 1.5 * backboneDict2['numMonomers'] * backboneDict2['bondLength'] + backboneDict2['Z1']
        backboneDict2['R2'] = scale * 70
    
        connectorDict12['displacement'] =  scale * 1.0 * backboneDict1['bondLength']
        connectorDict12['alpha'] = 0.0
        connectorDict12['beta'] = 160.0
    
        brushDict1['filename'] = 'RandomPolymer.txt'  

        brushDict1['name'] ='brush1'
        brushDict1['numMonomers'] = 10
        brushDict1['monomerNames'] = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
        brushDict1['indexOfFirstBrush'] = 1
        brushDict1['spacing'] = 1
        if Randomize:
            brushDict1['phaseRange'] = 360
            brushDict1['mode'] = 'Polymer'
        else:
            brushDict1['phaseRange'] = 15
            brushDict1['mode'] = 'PolymerAlternate'
        brushDict1['alpha1'] = 40
        brushDict1['alpha2'] = 50
        brushDict1['beta1'] = 145
        brushDict1['beta2'] = 155
        brushDict1['minDist'] = scale * 1.0
        brushDict1['bondLength'] = scale * 2.0
        brushDict1['Z1'] = scale * 2.0
        brushDict1['R1'] = scale * 10.0
        brushDict1['Z2'] = scale * brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
        brushDict1['R2'] = scale * 30.0
    
        brushDict2['filename'] = 'randomPolymer.txt'  
        brushDict2['name'] ='brush2'
        brushDict2['numMonomers'] = 5         
        brushDict2['monomerNames'] = ['Na', 'Mg', 'Al', 'Si', 'P']#, 'S', 'Cl', 'Ar', 'K', 'Ca']
        brushDict2['indexOfFirstBrush'] = 0
        brushDict2['spacing'] = 1
        if Randomize:
            brushDict2['phaseRange'] = 360
            brushDict2['mode'] = 'Polymer'
        else:
            brushDict2['phaseRange'] = 15
            brushDict2['mode'] = 'PolymerAlternate'
        
        brushDict2['alpha1'] = 40.0
        brushDict2['alpha2'] = 50.0
        brushDict2['beta1'] = 155.0
        brushDict2['beta2'] = 175.0
        brushDict2['minDist'] = scale * 1.0
        brushDict2['bondLength'] = scale * 2.0
        brushDict2['Z1'] = scale * 2.0
        brushDict2['R1'] = scale * 50.0
        brushDict2['Z2'] = scale * 1.5 * brushDict2['numMonomers'] * brushDict2['bondLength'] + brushDict2['Z1']
        brushDict2['R2'] = scale * 50.0
                  
        polymerBrushDict={}
        polymerBrushDict['backbones'] = [backboneDict1, backboneDict2]
        polymerBrushDict['brushes'] = [brushDict2, brushDict1]
        polymerBrushDict['connectors'] = [connectorDict12]
       
        self.xyz, self.names = GBP(polymerBrushDict)

        if save:
            fIO.saveXYZList(self.xyz, self.names, "SequestrationPLP.xyz")
    
    def getData(self):
        return self.xyz, self.names

class circles():
    
    def __init__(self, N):
        # create the NPack object.
        PlanePackBBG = SPPBBG('SurfacePackPlane.txt')

        xRange1 = -20
        xRange2 = 20
        yRange1 = -20
        yRange2 = 20
        distFromOrig = 0
        minDist = 2

        env = ['outersphere ' + str(xRange2)]

        # generate the building block
        PlanePackBB = PlanePackBBG.generateBuildingBlock(N, xRange1, xRange2, yRange1, yRange2, distFromOrig, minDist, envelopeList = env)
        
        newNames = 5 * ['H'] 
        newNames += 3 * ['He'] 
        newNames += 3 * ['Li'] 
        newNames += 6 * ['Be'] 
        newNames += 6 * ['B'] 
        newNames += 6 * ['C'] 
        newNames += 3 * ['N']
        newNames += 3 * ['O'] 
        newNames += 3 * ['F'] 
        newNames += 3 * ['Ne'] 
        newNames += 7 * ['Na']
    
        print(len(newNames))
        print(len(PlanePackBB.blockXYZVals))
    
        fIO.saveXYZList(PlanePackBB.blockXYZVals, newNames, "circles.xyz")

class SequestrationPLPCluster():
    def __init__(self):
        sphereBBG = SPSBBG('SurfacePackSphere.txt')
        numTauFilaments = 20
        ClusterNum  = 2000
        radius = 50
        theta1 = -90
        theta2 = 90
        phi1 = -135
        phi2 = 135
        minDistPLP = 0.5
        minDistTau = 15
        particleMinDist = 2.1
        bondLength = 2.0
        
        # generate the sphere base Points centered on origin by default
        SphereBB = sphereBBG.generateBuildingBlock(ClusterNum, radius, theta1, theta2, phi1, phi2, minDistPLP)
        SpherePoints = SphereBB.blockXYZVals

        # generate clusterNumber PLP Strand objects
        PLPStrands = [ self.SequestrationPLP(bondLength, particleMinDist) for _ in range(ClusterNum) ]
        directors = [ pos/np.linalg.norm(pos) for pos in SpherePoints] # points away from origin 

        # transform each PLP to the sphere base point aligned with director
        vesicleXYZ = [ coords.transformFromBlockFrameToLabFrame(director, pos, 0.0, np.array([0,0,1.0]), strand[0][0], strand[0]) for director, pos, strand in zip(directors, SpherePoints, PLPStrands) ]
 
        # concatentate output arrays        
        xyzVals = vesicleXYZ[0]
        names =  PLPStrands[0][1]
        for strand, vesicleStrand in zip(PLPStrands[1:], vesicleXYZ[1:]):
            xyzVals = np.concatenate( (xyzVals, vesicleStrand), 0)
            names = np.concatenate( (names, strand[1]), 0)

        # dump to file
        fIO.saveXYZList(xyzVals, names, "PLPVesicle.xyz")
        
        # generate the tau cluster in the sequestered place
        sphereDict = {'radius': radius,
                      'theta1': theta1,
                      'theta2':theta2,
                      'phi1': phi1,
                      'phi2': phi2,
                      'minDist': minDistTau } 
        
        ProtoFilamentCluster(numTauFilaments, 
                             radius/2, 
                             35, 
                             alphaMin=50,
                             alphaMax=70,
                             betaMin=120, 
                             betaMax=185,
                             bondLength=bondLength,
                             minDist=particleMinDist,
                             clusterMode='spherical',
                             sphericalParamsDict=sphereDict)


    def SequestrationPLP(self, bondLength):
        return np.array([[0.0, 0.0, 0.0], [0.0, 0.0, bondLength]]), np.array(['Pd', 'Pb'])

if __name__=='__main__':
    #circles(48)
    ProtoFilamentClusterAligned(10, 30, 3)
    #SequestrationPLPCluster()
    #SequestrationPLP()
    #ProtoFilamentCluster(10, 10, 35) # numFIlaments, cluster range, cluster Index
    # UnfoldedTau()
    # TauAmyloid('../pdb structures/5o3l.pdb', 60, 5 * 1.2005610112987413, 5 * 4.717551560924715)
    #GenericPLP([15], ['Pd'], ['VVVVVV'], 'Class1_PLP1_VQICYK')
    #GenericPLP([15], ['Pd'], ['DDVVVVVVDD'], 'Class1_PLP2_DRVQIVYKRR')
    #GenericPLP([15], ['Pd'], ['VVVVVVRR'], 'Class1_PLP3_VQIINKRR ')
    #GenericPLP([15], ['Pd'], ['VVVVVVVVVVVR'], 'Class1_PLP4_KVQIINKKLDRR ')
    #GenericPLP([15], ['Pd'], ['WWWWW'], 'Class1_PLP5_WMINK')
    #GenericPLP([15], ['Pd'], ['VVPVVV'], 'Class1_PLP6_VQPINK', doProlineKink=True)
    #GenericPLP([15, 15], ['Pd', 'Pb'], ['WWWWW', 'VVVVVV'], 'Class1_PLP7_WMINK_VQIVYK')
    #GenericPLP([15, 15], ['Pd', 'Pb'], ['WWPWWW', 'VVVVVVRR'], 'Class1_PLP8_VQPINK_VQIINKRR', doProlineKink=True)
    GenericPLP([15], ['Pd'], ['YYYYYYYYYYYYYRRR'], 'Class2_PLP9_YQQYQDATADEQGRRR')
    #GenericPLP([26, 4], ['Pd', 'Pb'], ['YYYYYYYYYYYYYRRR', 'AAAAAPRR'], 'Class3_PLP10_YQQYQDATADEQGRRR_ALAYIPRR', doProlineKink=True)
    #GenericPLP([15, 15], ['Pd', 'Pb'], ['AYQYAYQYAYQY', 'TITITIT'], 'mockup', doProlineKink=True)
    print("example done")