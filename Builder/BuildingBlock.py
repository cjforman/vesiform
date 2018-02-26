'''
Created on 16 Jan 2018

@author: chris
'''
import numpy as np
import Utilities.fileIO as fIO
import Utilities.cartesian as cart
import Utilities.coordSystems as coord

import copy as cp

class BuildingBlock():
    # A building block class is a wrapper for a list of xyz data points each point has an atomic name.
    # 
    # A series of member functions have been written to modify the data in the building block
    # in various ways. 
    #  
    # Such as : 
    #     Join the building blocks to another building block by stapling objects 
    #     together at connections points and returning a single new combined building block
    #     whose atomic co-ordinates, blockRefPoint and blockDirector are defined. 
    #
    #     Return the building block coords relative to it's centre of mass.
    #     Return the building block coords translated so it's reference point is at a target position.
    #     Return the building block coords oriented to it's body director is aligned with an external director. 
    #     Return the building block coords oriented, and positioned in a lab frame with a body axis rotation.
    #       
    #     Atoms can be renamed.
    # 
    #     Connection points can be added or removed.


    def __init__(self, blockXYZVals, blockAtomNames, blockConnectors, blockRefPoint, blockDirector):
        # the xyz Data must be an np.array of 1 x 3 np.arrays
        # connector data must a np.array of 1 x 3 np.arrays
        # refPoint is a 1 x 3 np. array
        # orientation is a 1 x 3 np.array
        # atomNames must be list same length as xyzVals
        # The blockDirector and the blockRefPoint are the director and ref points
        # used in the construction of the xyz points.
        self.blockAtomNames = blockAtomNames
        self.blockXYZVals = blockXYZVals
        self.blockConnectors = blockConnectors
        self.blockRefPoint = blockRefPoint
        self.blockDirectorHat = blockDirector/np.linalg.norm(blockDirector)
        self.blockRotation = 0
      
        if len(blockXYZVals)!=len(blockAtomNames):
            print "Building Block Error: Number of building block atom names must match number of atom vectors."
        
        if False in [len(vec)==3 for vec in blockXYZVals]:
            print "Warning: Building Block Error: All xyz vals must be 3 values long."

        if len(blockConnectors)>1:
            if False in [len(vec)==3 for vec in blockConnectors]:
                print "Warning: Building Block Error: All connectors must be 3 indices long."

    def append(self, BB):
        # simple add function for adding data from another building block to this one. Just copies 
        # the connectors, xyz vals and names but leaves refpoints and directors as
        # per the original building block. Adds an offset to the connector indices from the new block
        # which is the number of atoms in this building block
        self.addConnectors(BB.getAllConnectionIndices(), self.countAtoms())
        self.addAtoms(BB.getAtoms())

    def addBuildingBlock(self, newBB, selfConnector, newBBConnector, displacement, alpha1, beta1, alpha2, beta2, alpha3, newBlockRefPoint, newBlockOrientation):
        # Define what happens when two building blocks are added. 
        # Returns a building block with both sets of xyzvals incorporated and positioned correctly relative 
        # to each other. The atomic coordinates of the object being added are transformed into the
        # reference frame of the host object (the static object) using the angular and displacement 
        # parameters to correctly set the relative orientation.  
        #
        # 
        # A new combined building block is constructed and returned with the existing connection
        # points renumbered according to a single set of indices. The connectors used to join the blocks
        # are deleted, although the points used to join the building blocks are retained for future use.  
        # This aids with debugging and construction.   
        
        # Take the new building block and make a temporary copy of it which is 
        # oriented correctly (mobile building block). This leaves the original mobile building block intact.  
        mBB = self.makeConnection(selfConnector, newBB, newBBConnector, displacement, alpha1, beta1, alpha2, beta2, alpha3)

        # Create a new building block which will become the returned object
        # and starts as a carbon copy of the static object.
        retBB = cp.copy(self)

        # Generate a "staple" building block which is returned along side the new joined up
        # building block. This massively aids with debugging and system design and is simply
        # the six reference atoms used to connect the two building blocks in their final configuration.  
        stapleXYZ = np.append(self.getConnectionAtoms(selfConnector), mBB.getConnectionAtoms(newBBConnector), 0) 
        stapleNames = np.append(self.getConnectionNames(selfConnector), mBB.getConnectionNames(newBBConnector), 0)
        sBB = BuildingBlock(stapleXYZ, stapleNames, np.array([[0, 1, 2], [3, 4, 5]]), stapleXYZ[0], stapleXYZ[5] - stapleXYZ[0])

        # remove the connectors we used in makeConnection from the building blocks 
        # (they are now inside the joined block)
        retBB.removeConnector(selfConnector)
        mBB.removeConnector(newBBConnector)

        # count the number of atoms in retBB
        indexOffset = retBB.countAtoms()
    
        # add any remaining connectors from the mobile connector to the 
        # return building block renumbering them correctly for the new building block.
        retBB.addConnectors(mBB.getAllConnectionIndices(), indexOffset = indexOffset)
    
        # add the new atoms to the return buildingblocks
        retBB.addAtoms(mBB.getAtoms())
        
        # set the new ref point and orientation for the joint building block which are supplied by the user
        retBB.setBlockDirector(newBlockOrientation)
        retBB.setBlockRefPoint(newBlockRefPoint)
        
        # return the staple and the final composite building block
        return retBB, sBB 
    
    def makeConnection(self, selfConnector, newBB, newBBConnector, displacement, alpha1, beta1, alpha2, beta2, alpha3):
        # Returns a copy of the newBB that is correctly positioned relative to self, as defined
        # by the six parameters provided. Six degrees of freedom are necessary to specify the global
        # position and orientation of the mobile block. These are specified in terms of the dihedral and bond angles
        # formed by the juxtaposition of the two connectors used to connect the blocks. Input angles are in Radians. 
        # 
        # Six atoms are used to connect the blocks which are specified using 2 arrays of 3 indices that refer to the atoms in 
        # both of the building blocks. There are three atoms from each block: called 0s, 1s, 2s and 0m, 1m, 2m where m stands for 
        # mobile and s for static.
        # 
        # It is envisaged that 2m and 2s are the two atoms that form the connecting "bond" between the blocks so 2 is 
        # at the surface of a structure and 1, and 0 are increasingly deeper into the structures with m and s specifying which 
        # block they come from.
        # 
        # The algorithm never moves the static block and proceeds by first establishing the correct position for 2m in the 
        # mobile block, then 1m and then 0m. This is performed by a translation and three rotations.
        
        # Alpha1 is the dihedral angle formed by the atoms referenced by connector indices 0s 1s 2s 2m
        # beta1 is the bond angle formed by the atoms 1s, 2s and 2m.
        # displacement is the distance between 2s and 2m which form the connector points.
        # Together displacement, alpha1 and beta1 define where atom 2m goes in the static TNB frame formed 
        # from 0s, 1s and 2s, when this frame is positioned at 2s.
        #  
        # Thus the first translation places the mobile block so that atom 2m is at the location specified.
        #
        # Alpha2, beta2 and alpha3 then define the re-orientation of the entire second block.
        #
        # Beta2 is set first and is the the bond angle that should exist between the 2s 2m and 1m atoms.
        # This angle is set by first measuring the existing bond angle so formed after the translation.  
        # The second block is then rotated by the necessary angle about an axis through 2m which is Normal 
        # to the plane formed by 2s, 2m and 1m. (this is the N vector of the TNBframe define by the points 
        # 2s, 2m, 1m). 
        #
        # Alpha2 is set next and is the dihedral angle formed by the point 1s 2s 2m 1m. 
        # The existing angle is measured and the mobile block is then
        # rotated by the appropriate angle about a vector through 2m which is parallel to 2m - 2s. 
        #
        # point 1m is now positioned correctly.
        #  
        # THe final dihedral that must be set is the one defined by 2s, 2m, 1m and 0m, which is a rotation about the 2m to 1m axis.
        # The dihedral is measured and adjusted accordingly. 
        # Only moves point 0m and that means that the second block is now correctly positioned. 
        
        # first create a copy of the new BB so as not to modify the original building block.
        # mBB = mobile Building Block
        mBB = cp.copy(newBB)
        
        # extract the connection points for the two building blocks (self and mobile).  
        selfConnectXYZ = self.getConnectionAtoms(selfConnector)
        mobileConnectXYZ = mBB.getConnectionAtoms(newBBConnector)

        # construct a TNB frame from 0s, 1s and 2s. 
        TNB1 = coord.constructTNBFrame(selfConnectXYZ[0], selfConnectXYZ[1], selfConnectXYZ[2])
        
        # Compute the displacement of the mobileBlock so that the mobile BB's connection point 2m  
        # is at the right place in the TNB frame of the first connector when it is position with it's origin on the 
        # third point. The vector in this frame is defined by displacement, alpha1 and beta1.
        displacementVec = selfConnectXYZ[2] + displacement * coord.generateTNBVecXYZ(TNB1, beta1, alpha1) - mobileConnectXYZ[2]     
    
        # translate the mobile block so it's in the right place
        mBB.translateBB(displacementVec)
        
        
        ##### Rot 1: Set alpha 2 - dihedral about 2s to 2m axis ######
        # get the new repositioned mobile connector block
        mobileConnectXYZ = mBB.getConnectionAtoms(newBBConnector)
        
        # measure the dihedral between s1, s2, m2 and m1, keep order the same to maintain signs.
        alpha2Meas = coord.Dihedral(selfConnectXYZ[1], selfConnectXYZ[2], mobileConnectXYZ[2], mobileConnectXYZ[1])
        if alpha2Meas==None:
            alpha2Meas = alpha2
        
        # compute the alpha2rotation axis (unti vector from 2s to 2m)
        alpha2RotAxis = mobileConnectXYZ[2] - selfConnectXYZ[2]
        alpha2RotAxisHat = alpha2RotAxis/np.linalg.norm(alpha2RotAxis)
        
        # rotate the second block by an angle (alpha2 - alpha2Meas) about a vector through the point 2m which is parallel to the 2s 2m axis 
        # Only 1m and 0m should move from the static and mobile connector sets.
        mBB.rotateBBArbitraryAxis(mobileConnectXYZ[2], alpha2RotAxisHat, alpha2 - alpha2Meas)
        
        
        ##### Rot 2: Set Beta2- bond angle from point m2 ######
        # get the new repositioned mobile connector block
        mobileConnectXYZ = mBB.getConnectionAtoms(newBBConnector)
        
        # measure the bondangle between s2, m2 and m1
        beta2Meas = coord.bondAngle(selfConnectXYZ[2], mobileConnectXYZ[2], mobileConnectXYZ[1])
        
        # construct a TNB frame from s2, m2 and m1 
        TNB2 = coord.constructTNBFrame(selfConnectXYZ[2], mobileConnectXYZ[2], mobileConnectXYZ[1])
        
        # rotate mobile connector by angle (beta2 - beta2Meas) about a vector through the point 2m which is parallel with TNB2 frame N-axis
        # i.e. normal to plane of s2, m2 and m1  - the way normal is defined versus measurements of angle means to put a minus sign in there.
        mBB.rotateBBArbitraryAxis(mobileConnectXYZ[2], TNB2[1], -(beta2 - beta2Meas))
        
        
        ###### rot3:  set alpha3 - dihedral ab out 2m to 1m axis.
        # extract the new mobile connectors
        mobileConnectXYZ = mBB.getConnectionAtoms(newBBConnector)
        
        # construct rot3Axis from 2m to 1m. 
        rot3Axis = mobileConnectXYZ[1] - mobileConnectXYZ[2]
        rot3AxisHat = rot3Axis/np.linalg.norm(rot3Axis)

        # measure the dihedral betwee s2, m2, m1 and m0.
        alpha3Meas = coord.Dihedral(selfConnectXYZ[2], mobileConnectXYZ[2], mobileConnectXYZ[1], mobileConnectXYZ[0])
        
        # rotate the second block by an angle (alpha3 - alpha3Meas) about a vector through the point 2m which is parallel to the 2m 1m axis 
        # Only 0m should move to it's final place from the static and mobile connector sets.
        mBB.rotateBBArbitraryAxis(mobileConnectXYZ[2], rot3AxisHat, alpha3 - alpha3Meas)
        
        # return the rotated Building Block
        return mBB

    def replaceSingleAtomName(self, index, newName):
        self.blockAtomNames[index] = newName
        
    def replaceNames(self, newName):
        self.blockAtomNames = [newName] * self.countAtoms()         

    def setBlockRefPoint(self, newRefPoint):
        # changes the position of the reference point in the core representation
        self.blockRefPoint = newRefPoint
    
    def setBlockDirector(self, newBlockDirector):
        # sets a vector which defines a direction for the building block. 
        self.blockDirectorHat = newBlockDirector/np.linalg.norm(newBlockDirector)
    
    def getCOM(self):
        return cart.getCentreOfMass(self.blockXYZVals)
    
    def centerBB(self):
        # moves blockXYZValues so they are relative to their centre of mass. 
        self.blockXYZVals = self.blockXYZVals - self.getCOM()
        self.blockRefPoint = self.blockRefPoint - self.getCOM()
    
    def translateRefPointToTarget(self, targetPoint):
        # sets the block xyz vals translated such that the current lab 
        # ref point would be paced a the new target point
        self.blockXYZVals = self.blockXYZVals - self.blockRefPoint + targetPoint
        self.blockRefPoint = targetPoint  
        
    def orientToDirector(self, Director):
        # sets the block XYZ vals oriented to a given director.
        DirectorHat = Director/np.linalg.norm(Director)
        axisToRotate = np.cross(self.blockDirectorHat, DirectorHat)
        angleToRotate = np.arccos(np.dot(DirectorHat, self.blockDirectorHat))
        # rotates block xyzVals to be parallel with director - 
        # also rotates the orientor and refPoints.
        self.rotateBBArbitraryAxis(self.blockRefPoint, axisToRotate, angleToRotate)
        
    def translateBB(self, transVec):
        # returns the building block moved by specified vector - must also translate the reference point
        self.blockXYZVals  = self.blockXYZVals + transVec
        self.blockRefPoint = self.blockRefPoint + transVec
        
    def centerAtom(self, atomNum):
        # centers the building block with a particular atom at the origin.
        # also modifies the ref_point
        self.blockXYZVals  = self.blockXYZVals - self.xyzVals[atomNum]
        self.blockRefPoint = self.blockRefPoint - self.xyzVals[atomNum]
        
    def placeAtom(self, atomNum, newVec):
        # places the building block positions with a particular atom placed at a specific vector
        # also modifies the ref point 
        self.blockXYZVals = self.blockXYZVals - self.blockXYZVals[atomNum] + newVec
        self.blockRefPoint = self.blockRefPoint - self.blockXYZVals[atomNum] + newVec
        
    def rotateBBAtomAxis(self, atom1, atom2, angle):
        # sets the blockXYZVals to be rotated about an axis defined by two of the atoms
        # also modifies the refPoint and orienter
        axis = cart.vectorBetweenTwoPoints(self.blockXYZVals[atom1], self.blockXYZVals[atom2], direction='OneToTwo')
        self.blockXYZVals = [ cart.rotPAboutAxisAtPoint(p, self.blockXYZVals[atom1], axis, angle) for p in self.blockXYZVals]
        self.blockRefPoint = cart.rotPAboutAxisAtPoint(self.blockRefPoint, self.blockXYZVals[atom1], axis, angle)
        self.blockDirectorHat = cart.rotPAboutAxis(self.blockDirectorHat, self.blockXYZVals[atom1], axis, angle)
    
    def rotateBBArbitraryAxis(self, p0, axis, angle):
        # sets building block coords rotated about an axis through p0
        # also modifies the refPoint and orientation 
        self.blockXYZVals = [ cart.rotPAboutAxisAtPoint(p, p0, axis, angle) for p in self.blockXYZVals]
        self.blockRefPoint = cart.rotPAboutAxisAtPoint(self.blockRefPoint, p0, axis, angle)
        self.blockDirectorHat = cart.rotPAboutAxis(self.blockDirectorHat, axis, angle)
    
    def cloneBuildingBlock(self):
        return BuildingBlock(self.blockXYZVals, self.blockAtomNames, self.blockConnectors, self.blockRefPoint, self.blockDirectorHat) 
        
    def addAtoms(self, newAtoms):
        # Adds atoms and atomnames to core representation of the building 
        # defined in the tuple of the xyz vals and the atom names
        if len(self.blockXYZVals)==0:
            self.blockXYZVals = newAtoms[0]
            self.blockAtomNames = newAtoms[1]            
        else:
            self.blockXYZVals = np.concatenate((self.blockXYZVals, newAtoms[0]), axis=0)
            self.blockAtomNames = self.blockAtomNames + newAtoms[1]
    
    def getAtoms(self):
        return (self.blockXYZVals, self.blockAtomNames)
    
    def getAtomsXYZ(self):
        return self.blockXYZVals
    
    def getConnectionAtoms(self, connectionIndex):
        return np.array([self.blockXYZVals[atomIndex] for atomIndex in self.blockConnectors[connectionIndex]])
    
    def getConnectionNames(self, connectionIndex):
        return np.array([self.blockAtomNames[atomIndex] for atomIndex in self.blockConnectors[connectionIndex]])
        
    def getConnectionIndices(self, connectionIndex):
        return self.blockConnectors[connectionIndex]

    def getAllConnectionIndices(self):
        return self.blockConnectors
         
    def replaceConnector(self, connectionIndex, newConnectorArray):
        self.blockConnectors[connectionIndex] = newConnectorArray
        
    def addConnectors(self, connectors, indexOffset=0):
        # checks to see if we are adding connectors for first time (connectors list has no entry in it)
        if len(self.blockConnectors)==1: 
            if len(connectors[0]==0):
                self.blockConnectors = connectors + indexOffset
        else:
            # adds a list of connectors to the building block with an index offset to be applied to all the indices in the new connectors
            self.blockConnectors = np.append(self.blockConnectors, [ [connectorIndex + indexOffset  for connectorIndex in connector] for connector in connectors], 0)
        
    def removeConnector(self, connectorIndex):
        self.blockConnectors = np.delete(self.blockConnectors, connectorIndex, 0)
        
    def countAtoms(self):
        return len(self.blockXYZVals)

    def transformBBToLabFrame(self, labDirector, labRefPoint, labRotation):
        self.blockXYZVals = coord.transformFromBlockFrameToLabFrame(labDirector, labRefPoint, labRotation, self.blockDirectorHat, self.blockRefPoint, self.blockXYZVals)
        self.blockDirectorHat = coord.transformFromBlockFrameToLabFrame(labDirector, labRefPoint, labRotation, self.blockDirectorHat, self.blockRefPoint, [self.blockRefPoint + self.blockDirectorHat])[0] - labRefPoint
        self.blockRefPoint = labRefPoint

    def importBBK(self, fileroot):
        # replaces the core representation of the building block with the data stored in a BBK file.
        xyzVals, self.blockAtomNames = fIO.loadXYZ(fileroot + '.xyz')
        self.blockXYZVals = np.array([xyzVals])
        inputVals = fIO.readTextFile(fileroot+'.bbk')
        self.blockRefPoint = np.array([float(refVal) for refVal in inputVals[0].split()[2:]])
        director = np.array([float(refVal) for refVal in inputVals[1].split()[2:]])
        self.blockDirectorHat = director/np.linalg.norm(director)
        self.blockConnectors = np.array([int(val) for val in inputVals[2].split()[2:]])
        self.blockConnectors = np.reshape(self.blockConnectors, (self.blockConnectors/3, 3))
    
    def exportBBK(self, fileroot):
        fIO.saveXYZList(self.blockXYZVals, self.blockAtomNames, fileroot + '.xyz')
        line1 = "floatlist blockRefPoint " + str(self.blockRefPoint[0]) + " " + str(self.blockRefPoint[1]) + " " + str(self.blockRefPoint[2]) + "\n"
        line2 = "floatlist blockOrientation " + str(self.blockDirectorHat[0]) + " " + str(self.blockDirectorHat[1]) + " " + str(self.blockDirectorHat[2]) + "\n"
        lines = "intlist blockConnectors" 
        for connection in self.blockConnectors:
            lines += " " + str(connection[0]) + " " + str(connection[1]) + " " + str(connection[2])
        lines +='\n'
        fIO.writeTextFile([line1, line2, lines], fileroot + '.bbk')

if __name__=="__main__":

    from Library.peptideBackbone import peptideBackboneGenerator as PBG
        
    # set up an alpha helix generator object using the alphaHelix parameters
    helixGenerator = PBG("alphaHelix.txt")
     
    # specify the desired location and orientations of two specific helices
    numResidues = 16
    refPoint1 = np.array([0.0, 0.0, 0.0])
    director1 = np.array([0.0, 0.0, 1.0])
    refPoint2 = np.array([20.0, 0.0, 0.0])
    director2 = np.array([1.0, 0.0, 1.0])
    
    # generate helix1 using the generator function 
    helix1 = helixGenerator.generateBuildingBlock(numResidues, showBlockDirector=False)
    helix1.transformBBToLabFrame(director1, refPoint1, 0)
    
    # generate helix 2 using the generator function
    helix2 = helixGenerator.generateBuildingBlock(numResidues, showBlockDirector=False)
    helix2.transformBBToLabFrame(director2, refPoint2, 0)
     
    # save the helices to file
    helix1.exportBBK('alphaHelix1')
    helix2.exportBBK('alphaHelix2')
     
    # set up the angles and displacement between the two building blocks
    alpha1Inp = -45.0
    beta1Inp = 179.0  
    alpha2Inp = -0.0
    beta2Inp = 135.0  
    alpha3Inp = -0.0
    displacementInp = 15.0
     
    sConnectionIndex = 1
    mConnectionIndex = 0
     
    # connect the building blocks with the specified information
    BB3, SB = helix1.addBuildingBlock(helix2, sConnectionIndex, mConnectionIndex, 
                                      displacementInp,
                                      alpha1Inp * np.pi / 180,
                                      beta1Inp * np.pi / 180,
                                      alpha2Inp * np.pi / 180,
                                      beta2Inp * np.pi / 180,
                                      alpha3Inp * np.pi / 180,                               
                                      refPoint1, director1)
    
    SB.exportBBK("helixStaple12")
    BB3.exportBBK("jointHelix")
     
    # check the bond angles and dihedrals of the connecting bonds are the same 
    # in both the returned staple and the returned building block, and that they match the input values.
    XYZBB = BB3.getAtoms()[0]
    XYZSB = SB.getAtoms()[0]
    
    # figure out the indices in BB3 of the original connectors in the two helix building blocks
    sConnectorIndices = helix1.getConnectionIndices(sConnectionIndex)
    mConnectorIndices = [ ci + helix1.countAtoms() for ci in helix2.getConnectionIndices(mConnectionIndex)] 
     
    print "dead static connector: ", sConnectorIndices
    print "dead mobile connector: ", mConnectorIndices
     
    sXYZBB = [ XYZBB[index] for index in sConnectorIndices ]
    mXYZBB = [ XYZBB[index] for index in mConnectorIndices ]
    sXYZSB = [ XYZSB[0], XYZSB[1], XYZSB[2] ]
    mXYZSB = [ XYZSB[3], XYZSB[4], XYZSB[5] ] 
    
    # measure the angles for the new BBuilding block - should be the same as the input angles
    alpha1BB = coord.Dihedral(sXYZBB[0], sXYZBB[1], sXYZBB[2], mXYZBB[2]) * 180 / np.pi 
    beta1BB = coord.bondAngle(sXYZBB[1], sXYZBB[2], mXYZBB[2]) * 180 / np.pi
    alpha2BB = coord.Dihedral(sXYZBB[1], sXYZBB[2], mXYZBB[2], mXYZBB[1]) * 180 / np.pi
    beta2BB = coord.bondAngle(sXYZBB[2], mXYZBB[2], mXYZBB[1]) * 180 / np.pi
    alpha3BB = coord.Dihedral(sXYZBB[2], mXYZBB[2], mXYZBB[1], mXYZBB[0]) * 180 / np.pi
    displacementBB = np.linalg.norm(mXYZBB[2]- sXYZBB[2])
    
    # measure the angles for the staple Building block - should be the same as the input angles
    alpha1SB = coord.Dihedral(sXYZSB[0], sXYZSB[1], sXYZSB[2], mXYZSB[2]) * 180 / np.pi
    beta1SB = coord.bondAngle(sXYZSB[1], sXYZSB[2], mXYZSB[2]) * 180 / np.pi
    alpha2SB = coord.Dihedral(sXYZSB[1], sXYZSB[2], mXYZSB[2], mXYZSB[1]) * 180 / np.pi
    beta2SB = coord.bondAngle(sXYZSB[2], mXYZSB[2], mXYZSB[1]) * 180 / np.pi
    alpha3SB = coord.Dihedral(sXYZSB[2], mXYZSB[2], mXYZSB[1], mXYZSB[0]) * 180 / np.pi
    displacementSB = np.linalg.norm(mXYZSB[2]- sXYZSB[2])
    
    print "input : ", alpha1Inp, beta1Inp, alpha2Inp, beta2Inp, alpha3Inp, displacementInp
    print "new BB: ", alpha1BB, beta1BB, alpha2BB, beta2BB, alpha3BB, displacementBB
    print "new SB: ", alpha1SB, beta1SB, alpha2SB, beta2SB, alpha3SB, displacementSB
    
    
    print "done"