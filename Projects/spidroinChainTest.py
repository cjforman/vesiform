from Projects.chainOfChains import chainOfChainsGenerator as CCGen
import Utilities.fileIO as fIO 

filename = "chainOfChains.txt"
nameList, xyzVals = fIO.loadXYZ("species2SpidroinHPin.xyz")
    
# create the backbone generator object using static file parameters
chainOfChainsGen = CCGen(filename)

# generate backbone realtime parameters
minDist = 1.0
residues = [19, 15] * 131 + [19]
radii = [minDist * 6] * 263
pointsA = [ xyzVal for xyzVal,name in zip(xyzVals, nameList) if name == 'N' ]
pointsB = [ xyzVal for xyzVal,name in zip(xyzVals, nameList) if name == 'C' ]
    
fIO.saveXYZList(pointsA + pointsB, ['Ca'] * len(pointsA) + ['O'] * len(pointsB), "labPoints.xyz")
    
cOfChainsBB = chainOfChainsGen.generateBuildingBlock(residues, pointsA, pointsB, radii, minDist, visualiseEnvelope=(0, 50, 'envelope.xyz'))
cOfChainsBB.exportBBK("species2Chain")
    
print "chainOfChains done"
    
    
        
    