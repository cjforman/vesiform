import sys
import numpy as np
from Projects.Spidroin3 import spidroinProteinGenerator as SPG

filename = sys.argv[1]
# individual spidroins
spidroinProteinGenerator = SPG(filename)
species1 = 'SP1'
startPoint = np.array([0.0, -0.0, 0.0])
direction = np.array([-0.0, 0.0, 1.0])
rotation = 0.0
minDist = 1.0

SpidroinSpecies1BB = spidroinProteinGenerator.generateBuildingBlock( species1,
                                                                     minDist,
                                                                     showBlockDirector=False)  
SpidroinSpecies1BB.exportBBK('species1Spidroin.xyz')

print("Spidroin Species 1 Done.")
        
        
    