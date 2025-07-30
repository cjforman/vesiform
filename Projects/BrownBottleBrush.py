import numpy as np
import copy as cp
import random as rnd
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
import Utilities.cartesian as carts 
import Utilities.coordSystems as coords  
import Utilities.fileIO as fIO

# takes a dictionary which contains polymer definitions - a backbone and the side chains of a bottle brush
def BrushPolymer(polymerDict):
    # decorate the polymer backbone with brushes using the brush information
    return decoratePolymerBackbone(polymerDict['brush'], makeBackboneBB(polymerDict['backbone']))

       
def makeBackboneBB(backbone):
    ''' Takes the backbone definition and returns as back bone polymer building block object defined in its own block frame ''' 
    
    # create a backbone generator object for the polymer
    backboneBBG = RPPBBG(backbone['filename'])
    
    # create the starting point of the polymer
    try:
        startPointBackBone = backbone['pointToStart']
    except KeyError:
        startPointBackBone = np.array([ 0.0, 0.0, backbone['Z1']])

    # generate the frustrum containing the polymer
    envelopeList = ['frustum ' + str(backbone['Z1']) + ' ' + str(backbone['R1']) + ' ' + str(backbone['Z2']) + ' ' + str(backbone['R2'])]

    try:
        pointsToAvoid = backbone['PointsToAvoid']
    except KeyError:
        pointsToAvoid = []

    # assume failure
    rejectPolymer = True

    # loop until we get a viable polymer    
    while rejectPolymer:
        # generate polymerBB of appropriate length. 
        backboneBB = backboneBBG.generateBuildingBlock( backbone['numMonomers'],  
                                                        startPointBackBone,
                                                        backbone['alpha1'],
                                                        backbone['alpha2'],
                                                        backbone['beta1'],
                                                        backbone['beta2'],
                                                        backbone['minDist'],
                                                        backbone['bondLength'],
                                                        envelopeList=envelopeList,
                                                        visualiseEnvelope=(0, 2*backbone['R2'] , backbone['name'] + '_envelope.xyz'),
                                                        pointsToAvoid=pointsToAvoid)
        
        backboneBB.blockAtomNames = backbone['monomerNames']
    
        # check that Z coord of point 0 is below z last point of polymer  otherwise reject polymer and choose another one.  
        if backboneBB.blockXYZVals[0][2] < backboneBB.blockXYZVals[-1][2]:
            rejectPolymer = False
        else:
            print("Polymer Rejected. Doing another.")

    return backboneBB


def decoratePolymerBackbone(brush, backboneInfo):
    # the dictionary defines the nature of the side chains on the main polymer.
    # The mode Keyword can contain any combination of the following control phrases in any order
    # There can be multiple side chains per building backbone entry
    
    # unpack input
    points = backboneInfo.blockXYZVals
    backboneNames = backboneInfo.blockAtomNames
    backboneLength = len(points)


    # set up output - it is important to ensure the correlation between name and xyz points is always maintained. 
    # begin the final output which is just a list of points and point names for the xyz file 
    polymerXYZ = cp.copy(points)
    polymerNames = cp.copy(backboneNames)

    # compute local TNB frame for each point in the back bone
    TNBs = [ coords.constructTNBFrame(p1, p2, p3) for pi, p2, p3 in zip(points[0:-2], points[1:-1], points[2:])


    # for each point in the backbone, define a bunch of brushes.
    for pointNum, point in enumerate(points):
        
        #for each brush at the backbone point, define the right number of brushes. 
        for brushNum in arange(brush['numBrushesPerMonomer']):
     
            # Create a generator for each brush. 
            brushBBG = RPPBBG(brush['filename'])

            # Set the angle of each brush around the director defined relative to the Binormal and normal in the TNB frame for that polymer.
            # The range of angles is defined in the dictionary. Evenly distribute the brushes around the phase range, for each monomer take a phase step, allowing staggered brushes on adjacent phases. 
            brushPhaseAngle = brushNum * brush['phaseRange'] * np.pi/180.0 + pointNum * brush['PhaseStepPerMonomer'] * np.pi/180.0

            # generate a director in direction of phase angles
            targetBrushDirector = np.cos(brushPhaseAngle) * TNB[2] + np.sin(brushPhaseAngle) * TNB[1] 
            targetBrushDirectorNorm = targetBrushDirector/np.linalg.norm(targetBrushDirector)
                    
            # generate the local brush frustum
            envelopeList = ['frustum ' + str(brush['Z1'] - brush['minDist']) + ' ' + str(brush['R1']) + ' ' + str(brush['Z2']) + ' ' + str(brush['R2'])]

            # add a brush, avoiding all points already in the brush        
            brushBB = brushBBG.generateBuildingBlock( brush['numMonomers'],
                                                      np.array([ 0.0, 0.0, brush['Z1']]),
                                                      brush['alpha1'],
                                                      brush['alpha2'],
                                                      brush['beta1'],
                                                      brush['beta2'],
                                                      brush['minDist'],
                                                      brush['bondLength'],
                                                      envelopeList=envelopeList,
                                                      visualiseEnvelope=(0, 100, brush['name'] + '_envelope.xyz'),
                                                      pointsToAvoid = polymerXYZ)
            
            # extract the xyzs and names from the BB
            brushXYZ = brushBB.blockXYZVals
            brushNames = brushBB.blockAtomNames
                        
            # overrides the names created in the building block generator with monomerNames in BrushDictionary
            try:
                if len(brushXYZ)==len(brush['monomerNames']):
                    brushNames = brush['monomerNames']
                else:
                    print("Warning: Number of XYZ vals and number of specified monomer names inconsistent. Not overriding generator.") 
            except KeyError:
                pass

            # work out the actual brush director from the brushXYZ vals
            if brush['numMonomers']>2:
                brushDirector = coords.axisFromHelix(brushXYZ)
            elif brush['numMonomers']==2:
                brushDirector = brushXYZ[1] - brushXYZ[0]
                brushDirector = brushDirector/np.linalg.norm(brushDirector)
            else:
                brushDirector = np.array([0.0, 0.0, 1.0])

            

            # rotate, translate the coords and concatenate the brush positions and brush names to the output arrays
            polymerXYZ = np.concatenate( (polymerXYZ, 
                                          coords.transformFromBlockFrameToLabFrame( targetBrushDirector,
                                                                                    point + brush['bondLength'] * targetBrushDirector,
                                                                                    0.0,
                                                                                    brushDirector,
                                                                                    brushXYZ[0],
                                                                                    brushXYZ)), 0)
            polymerNames = np.concatenate( (polymerNames, brushNames), 0)
        
    return polymerXYZ, polymerNames, polymerBrushesPerBlock

if __name__=="__main__":

    # set up the dictionaries
    backboneDict1 = {}
    brushDict1 = {}
    
    # pack the dictionaries
    backboneDict1['filename'] = 'RandomPolymer.txt'  
    backboneDict1['name'] = 'BackBone'
    backboneDict1['numMonomers'] = 18 
    backboneDict1['monomerNames'] = backboneDict1['numMonomers'] * ['P']
    backboneDict1['alpha1'] = 40
    backboneDict1['alpha2'] = 50
    backboneDict1['beta1'] = 165
    backboneDict1['beta2'] = 185
    backboneDict1['minDist'] = 2.0
    backboneDict1['bondLength'] = 2.0
    backboneDict1['Z1'] = 2
    backboneDict1['R1'] = 20
    backboneDict1['Z2'] = 1.5 *backboneDict1['numMonomers'] * backboneDict1['bondLength'] + backboneDict1['Z1']
    backboneDict1['R2'] = 70

    brushDict1['filename'] = 'RandomPolymer.txt'  
    brushDict1['name'] ='brush1'
    brushDict1['numMonomers'] = 1 
    brushDict1['monomerNames'] = brushDict1['numMonomers'] * ['Fe']
    brushDict1['phaseRange'] = 5.0 
    brushDict1['alpha1'] = 40
    brushDict1['alpha2'] = 50
    brushDict1['beta1'] = 145
    brushDict1['beta2'] = 155
    brushDict1['minDist'] = 1.0
    brushDict1['bondLength'] = 2.5
    brushDict1['Z1'] = 2.0
    brushDict1['R1'] = 10.0
    brushDict1['Z2'] = brushDict1['numMonomers'] * brushDict1['bondLength'] + brushDict1['Z1']
    brushDict1['R2'] = 30.0

              
    polymerBrushDict={}
    polymerBrushDict['backbone'] = backboneDict1
    polymerBrushDict['brush'] = brushDict1
       
    strand = BrushPolymer(polymerBrushDict)
    fIO.saveXYZList(strand[0], strand[1], "Polymer.xyz")
    
    print("example done")