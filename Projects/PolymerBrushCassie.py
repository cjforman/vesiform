import numpy as np
import random as rnd
from Library.peptideBackbone import peptideBackboneGenerator as PBG
from Library.randomPolymer import RandomPolymerPackBBG as RPPBBG
from Projects.PolymerBrush import polymerBrush as PB
import Utilities.coordSystems as coords  
import Utilities.fileIO as fIO


if __name__=="__main__":

    brushDict={}
    # pack the brush dictionary
    brushDict['filenameBlock'] = 'RandomPolymer.txt'  
    brushDict['filenameBrush'] = 'RandomPolymer.txt' # use 'betastrand.txt' for mode peptide
    brushDict['mode'] = 'PolymerRandom'
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
      
    (brushDir, brushPoint, xyz, names) = PB(brushDict)
    fIO.saveXYZList(xyz, names, "brush.xyz")

    print("example done")