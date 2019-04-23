import numpy as np
import random as rnd
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
from Projects.BlockCopolymer import makeBlockCopolymer as MBCP
import Utilities.coordSystems as coords  
import Utilities.fileIO as fIO

def polymerBrush(brushDict, mode="RandomPolymer"):

    # unpack the brush dictionary
    filenameRandom = brushDict['filenameBlock']
    filenameBrush = brushDict['filenameBrush']
    mode = brushDict['mode']
    
    ABlock = brushDict['ABlock']
    num_A = ABlock['num'] 
    alpha1_A = ABlock['alpha1']
    alpha2_A = ABlock['alpha2'] 
    beta1_A = ABlock['beta1'] 
    beta2_A = ABlock['beta2'] 
    minDist_A = ABlock['minDist'] 
    bondLength_A = ABlock['bondLength']
    Z1_A = ABlock['Z1'] 
    R1_A = ABlock['R1'] 
    Z2_A = ABlock['Z2'] 
    R2_A = ABlock['R2'] 
    
    BBlock = brushDict['BBlock']
    num_B = BBlock['num'] 
    alpha1_B = BBlock['alpha1']
    alpha2_B = BBlock['alpha2'] 
    beta1_B = BBlock['beta1'] 
    beta2_B = BBlock['beta2'] 
    minDist_B = BBlock['minDist'] 
    bondLength_B = BBlock['bondLength']
    Z1_B = BBlock['Z1'] 
    R1_B = BBlock['R1'] 
    Z2_B = BBlock['Z2'] 
    R2_B = BBlock['R2']
    
    brushBlock = brushDict['brushBlock']
    num_brush = brushBlock['num'] 
    alpha1_brush = brushBlock['alpha1']
    alpha2_brush = brushBlock['alpha2'] 
    beta1_brush = brushBlock['beta1'] 
    beta2_brush = brushBlock['beta2'] 
    minDist_brush = brushBlock['minDist'] 
    bondLength_brush = brushBlock['bondLength']
    Z1_brush = brushBlock['Z1'] 
    R1_brush = brushBlock['R1'] 
    Z2_brush = brushBlock['Z2'] 
    R2_brush = brushBlock['R2']
    

    # Mode word must specify one of Polymer or Peptide.  Can also include modifiers Sheet or Random to specify phase of brushes around polymer backbone axis.
    # defaults to RandomPolymer.  Can supply brush polymer parameters via parameter polymer params. These are alpha1, alpha2, beta1, beta2, minDist and bond length. 
    # atom name supplied in filenameBrush.  This is a mess. I know. I don't care.  
    
    # if one of Peptide or Polymer is not in mode then add Polymer to the end of mode. Anything else is dandy. If you specify both peptide and polymer, then Peptide will override.
    if not ("Peptide" in mode or "Polymer" in mode):
        mode = mode + "Polymer"
    
    # generate the block Copolymer backbone
    (polymerXyzVals, polymerNames) = MBCP(num_A, num_B,   Z1_A, R1_A, Z2_A, R2_A, alpha1_A, alpha2_A, beta1_A, beta2_A, minDist_A, bondLength_A, 
                                                        Z1_B, R1_B, Z2_B, R2_B, alpha1_B, alpha2_B, beta1_B, beta2_B, minDist_B, bondLength_B, filenameRandom)
    
    # get axis of block A
    BlockAAxis = coords.axisFromHelix(polymerXyzVals[0:num_A])
    BlockAAxisNorm = BlockAAxis/np.linalg.norm(BlockAAxis) 
    
    # get an orthogonal vector to BlockAAxis which is random. 
    randDirXYZ = BlockAAxis
    randXYZIsAxizXYZ = True
    
    while randXYZIsAxizXYZ: 
        randVecDir = coords.pickRandomPointOnUnitSphere()
        randDirXYZ = coords.sphericalPolar2XYZ(np.array([1.0, randVecDir[0], randVecDir[1]]))
        randXYZIsAxizXYZ = not False in [ np.abs(a - b) < 1e-7 for a, b in zip(randDirXYZ, BlockAAxis) ]
        
    OrthVec = randDirXYZ - np.dot(BlockAAxisNorm, randDirXYZ) * BlockAAxisNorm
    OrthVec1Norm = OrthVec/np.linalg.norm(OrthVec) 
    OrthVec2Norm  = np.cross(OrthVec1Norm, BlockAAxisNorm)

    # now have orthonormal basis for the polymer backbone for the brush part of system 
    
    # create the brush generator. Default is polymer mode. If peptide is included in mode word then over ride to use peptide instead 
    if "Peptide" in mode:
        # create a peptide generator using supplied filename to specify parameters of peptide - filename overides polymerParams
        brushGenerator = PBG(filenameBrush)
        brushObject = brushGenerator.generateBuildingBlock(num_brush) # only need to create this once for peptides as all are the same
    else:
        brushGenerator = RPPBBG(filenameBrush)
    
    # choose the phase angles of each brush
    # initially make the phase angle a little bit random
    brushPhaseAngles = [ rnd.uniform(0, 0.2 * np.pi) for _ in range(0, num_A) ]

    # if Random is in mode then make it very random
    if "Random" in mode:
        brushPhaseAngles = [ rnd.uniform(0, 2 * np.pi) for _ in range(0, num_A) ] 
    
    # if sheet is in mode then make phase zero.
    if "Sheet" in mode:
        brushPhaseAngles = [ 0.0 for _ in range(0, num_A) ]

    # generate directors in direction of phase angles
    brushDirectors = [ np.cos(angle) * OrthVec1Norm + np.sin(angle) * OrthVec2Norm for angle in brushPhaseAngles] 
    brushDirectorsNorm = [ d/np.linalg.norm(d) for d in brushDirectors]

    # for each of the directors (defined by the length of block A) figure out the final xyz vals    
    for point, labDirector in zip(polymerXyzVals[0:num_A], brushDirectorsNorm):
        if "Polymer" in mode and not "Peptide" in mode:
            # if we're doing a polymer then generate a new polymer with the given polymer parameters for each pass.
            polyStart = np.array([ 0.0, 0.0, Z1_brush])
            envelopeList = ['frustum ' + str(Z1_brush - minDist_brush) + ' ' + str(R1_brush) + ' ' + str(Z2_brush) + ' ' + str(R2_brush)]
            
            brushObject = brushGenerator.generateBuildingBlock( num_brush,
                                                                polyStart,
                                                                alpha1_brush,
                                                                alpha2_brush,
                                                                beta1_brush,
                                                                beta2_brush,
                                                                minDist_brush,
                                                                bondLength_brush,
                                                                envelopeList = envelopeList)
        
        newBrushXYZ = coords.transformFromBlockFrameToLabFrame(labDirector, point + minDist_A * labDirector, 0.0, brushObject.blockDirectorHat, brushObject.blockXYZVals[0], brushObject.blockXYZVals)
        polymerXyzVals = np.concatenate((polymerXyzVals, newBrushXYZ), 0) 
        polymerNames = np.concatenate( (polymerNames, brushObject.blockAtomNames), 0 )


    # direction of brush is z-axis.  
    # reference point is the first of the B polymers, which is the numAth entry in the list of xyzVals
    return (np.array([0.0, 0.0, 1.0]), polymerXyzVals[num_A], polymerXyzVals, polymerNames)

if __name__=="__main__":

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
    brushDict['ABlock']['minDist'] = 2.0
    brushDict['ABlock']['bondLength'] = 3.0
    brushDict['ABlock']['Z1'] = 2
    brushDict['ABlock']['R1'] = 20
    brushDict['ABlock']['Z2'] = brushDict['ABlock']['num'] * brushDict['ABlock']['bondLength'] + brushDict['ABlock']['Z1']
    brushDict['ABlock']['R2'] = 20


    brushDict['BBlock']['num'] = 40 
    brushDict['BBlock']['alpha1'] = 40.0
    brushDict['BBlock']['alpha2'] = 50.0
    brushDict['BBlock']['beta1'] = 145.0
    brushDict['BBlock']['beta2'] = 155.0
    brushDict['BBlock']['minDist'] = 2.0
    brushDict['BBlock']['bondLength'] = 3.0
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
      
    (brushDir, brushPoint, xyz, names) = polymerBrush(brushDict)
    fIO.saveXYZList(xyz, names, "brush.xyz")

    print("example done")