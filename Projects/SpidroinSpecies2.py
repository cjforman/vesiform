import sys
import numpy as np
from Projects.Spidroin3 import spidroinProteinGenerator as SPG

filename = sys.argv[1]

# individual spidroins
spidroinProteinGenerator = SPG(filename)
species2 = 'SP2'
startPoint = np.array([0.0, -0.0, 0.0])
direction = np.array([-0.0, 0.0, 1.0])
rotation = 0.0
minDist = 1.0

SpidroinSpecies2BB = spidroinProteinGenerator.generateBuildingBlock( species2,
                                                                     minDist,
                                                                     showBlockDirector=False, sheared=True)  
SpidroinSpecies2BB.exportBBK('species2Spidroin.xyz')

print("Spidroin Species 2 Done.")
        
        
    