import Utilities.fileIO as fIO
import Utilities.coordSystems as coords

blockNames, xyzVals = fIO.loadXYZ('species2Spidroin.hPin.xyz')
eAtomNames, ePositions, eSizes, eRs, eRotVecs = coords.convertTriplesToEllipsoids(blockNames, xyzVals, 1.0, 1.0)
fIO.saveEllipsoidXYZList(eAtomNames, ePositions, eSizes, eRs, eRotVecs, 'species2Spidroin.e.xyz')

print "example done"