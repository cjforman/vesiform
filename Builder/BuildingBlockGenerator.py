import sys
from Utilities.keyProc import keyProc
import Utilities.cartesian as cart
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO

import numpy as np
import random as rnd
from Builder.BuildingBlock import BuildingBlock 

class BuildingBlockGenerator(keyProc):
    # This is the class definition for the base building block generator.
    #
    # Any class which inherits this definition will be capable of generating 
    # building blocks used by the vesiform builder system. In essence it is an interface
    # between the vesiform builder system and any arbitrary functional method for 
    # generating xyz values, such as Npack objects or the PDB reader.
    #
    # Static parameters of the building block type are loaded from file via 
    # the keyproc system when the generator is constructed. It is envisaged that
    # for each building block type used in a project these file driven parameters 
    # will never change across a project.
    #
    # Dynamic features such as length, number of points etc 
    # are passed as parameters in to the generator function when it is called.
    #   
    # The class provides functions for generating a list of xyz points, 
    # an associated list of names as well defining the connectors of the building block
    # and the orientation and principle reference points of a building block.
    # 
    # Finally there is a function for generating and returning a building block once all this
    # information has been computed.
     
    def __init__(self, paramFilename):
        # initialise the parameter dictionary
        keyProc.__init__(self, paramFilename)
    
    def initialiseParameters(self):

        keyProc.initialiseParameters(self)

        self.verbose = self.getParam('verbose')
        self.dumpInterimFiles = self.getParam('dumpInterimFiles')

        if self.noLoadErrors == False:            
            print("Critical Parameters are undefined for BuildingBlockGenerator")
            sys.exit()        

    def generateBuildingBlock(self, numPoints, minDist, showBlockDirector=False, visualiseEnvelope=(0, 20, 'envelope.xyz'), envelopeList=['None'], pointsToAvoid=[], defaultBlockRefPoint=None, warp=False, **kwds):
        # Returns a building block object with numPoints xyz values and names.  
        # A reference point and principal direction of the building block are defined in the block Ref Frame.
        
        if defaultBlockRefPoint is not None:
            self.blockRefPoint = defaultBlockRefPoint
            
        # store num points
        self.numPoints = numPoints
        
        # store min Dist
        self.minDist = minDist

        # envelopeList
        self.parseEnvelopeList(envelopeList)

        # visualise envelope if requested
        if visualiseEnvelope[0]>0:
            print("Visualising envelope")
            self.visualiseEnvelopeFunc(visualiseEnvelope[0], visualiseEnvelope[1], visualiseEnvelope[1], visualiseEnvelope[1], visualiseEnvelope[2])

        # from the points to avoid list remove those points that would fail the envelope test
        self.pointsToAvoid=[]
        indicesOfPointsToRemoveFromPointsToAvoid = self.checkAllPointsInBounds(pointsToAvoid)
        if len(indicesOfPointsToRemoveFromPointsToAvoid)>0:
            self.pointsToAvoid = [ pointsToAvoid[index] for index in range(0, len(pointsToAvoid)) if not index in indicesOfPointsToRemoveFromPointsToAvoid]
                

        # save the pointsToAvoid if we're dumping files
        if self.dumpInterimFiles==1:
            fIO.saveXYZ(self.pointsToAvoid, 'Li', "pointsToAvoid.xyz")
        
        # generate the xyz positions
        self.buildingBlockXYZ = self.generateBuildingBlockXYZ()
        
        # warp the xyz positions if required
        if warp:
            fIO.saveXYZ(self.buildingBlockXYZ, 'N', "preWarp.xyz")
            self.buildBlockXYZ = self.warpBuildingBlockXYZ(**kwds)
            fIO.saveXYZ(self.buildingBlockXYZ, 'H', "postWarp.xyz")
        
        # generate the BuildingBlock reference point
        self.blockRefPoint = self.generateBuildingBlockRefPoint()
        
        # generate the BuildingBlock director unit vector
        self.blockDirectorHat = self.generateBuildingBlockDirector()
        
        # generate the BuildingBlock names
        self.blockNames = self.generateBuildingBlockNames()

        if showBlockDirector:
            # This is is a little hack for visualising the building Block Director in 
            # blender for debugging purposes. Insert this sequence into the middle of the chain (at point 15)
            # because the end of the chain is used to connect coils etc.
            insertIndex = int(np.floor(len(self.buildingBlockXYZ)/2))
            [self.buildingBlockXYZ.insert(insertIndex, self.blockRefPoint + i * self.blockDirectorHat)  for i in range(0, 10)]
            [self.blockNames.insert(insertIndex, 'CA') for i in range(0, 10)]

        # generate the BuildingBlock connector indices
        self.blockConnectors = self.generateBuildingBlockConnectors()

        # call the function to construct and return a building block from the 
        # information in the generator object. 
        return self.getBuildingBlock()
    
    # a function to warp a set of coordinates according to a transformation defined 
    # in the input param file
    def warpBuildingBlockXYZ(self, **kwds):
        warpType = kwds['warpType']
        
        if warpType=='Cylindrical':
            try:
                self.buildingBlockXYZ = self.CylinderWarp(self.buildingBlockXYZ,
                                                          kwds['CylinderRadius'],
                                                          kwds['CylinderHeight'],
                                                          kwds['wrapNumber'],
                                                          axis=np.array([ float(component) for component in kwds['initialAxis']]),
                                                          basePoint=np.array([ float(component) for component in kwds['basePoint']]),
                                                          zeroPoint=np.array([ float(component) for component in kwds['zeroPoint']]))
            except KeyError as e:
                print("Missing parameter:", e, " need cylinderRadius, cylinderHeight, initialiAxis, basePoint, zeroPoint.")
                sys.exit()
    
    # Computes the coordinates of each input point in cylindrical coords (rho, phi, z) relative to a cylinder 
    # defined by an axis, positioned at basePosition and a zero point (to define a zero phi). 
    # The Z parameter along this cylindrical system is taken to be the arclength position of a space curve along which we will 
    # warp the entire set of points. THe relative rho and phi positions define the azimuth and distance of the point from the new space curve in a plane
    # perpendicular to the space curve at each point. 
    def CylinderWarp(self, posList, radius, height, wrapNumber, axis=np.array([0.0 ,0.0, 1.0]), basePoint=np.array([0.0, 0.0, 0.0]), zeroPoint=np.array([1.0, 0.0, 0.0])): 
       
        # compute the relative coords of each given point in the global cylindrical polar system defined by the input params. 
        posListCyl = [ coords.XYZ2Cyl(pos, axis=axis, basePoint=basePoint, zeroPoint=zeroPoint) for pos in posList ]
        
        # fIO.saveXYZ(posListCyl, 'C', 'cylPolCoords.xyz')
            
        # now have rho, phi and z of each atom. z acts as arc length parameter for a space curve, and rho and phi the relative position in a plane
        # perpendicular to space curve.         
        return self.computeNewPos(posListCyl, radius, height, wrapNumber)  
        
    def computeNewPos(self, cylPosList, R, H, N):
        # compute a scaling constant alpha, which comes from parameterizing the helix via arc length.
        # This is necessary so that the spacing of the arclength parameters along the helix is the same as the 
        # spacing along the axis of the unwrapped chain
        # t = s/( (2 N pi R)^2 + h^2)^0.5  t goes from 0 to 1, which maps a helix around a cylinder
        # of height H, radius R with wrapping number N.
        alpha = np.power( np.power( 2.0 * np.pi * float(N) * R, 2 ) + np.power( H, 2 ), 0.5 ) 

        # debugging code to replace input cylPos with a more easily understood system.
        # cylPosList = [  np.array([0.1, 2.0 * np.pi * phi, z] ) for phi, z in zip(np.linspace(0, 1, len(cylPosList) ), np.linspace(0, 10, len(cylPosList) ) ) ]
        # fIO.saveXYZ(cylPosList, 'Pt', 'cylPosListTest.xyz')
        
        # more debugging code to help figure out difference between t and arc parameters
        # tTerms = [ np.array([ R * np.cos( 2.0 * np.pi * float(N) * t), 
        #                      R * np.sin( 2.0 * np.pi * float(N) * t), 
        #                      H * t]) for t in np.linspace(0.0, 1.0, 50)]
        # fIO.saveXYZ(tTerms, 'Pt', 'tTerms.xyz')

        # Compute the base point along the arc for each position. This is the XYZ coords of 
        # the origin of each TNB frame along the arc.
        basePointList = [ np.array([ R * np.cos( 2.0 * np.pi * float(N) * pos[2] / alpha ), 
                                     R * np.sin( 2.0 * np.pi * float(N) * pos[2] / alpha ), 
                                     H * pos[2] / alpha]) for pos in cylPosList ]

        # fIO.saveXYZ(basePointList, 'Pt', 'basePointList.xyz')

        # compute the tangent vector at the current arclength parameter
        TList = [ np.array([ - 2.0 * np.pi * float(N) * R * np.sin( 2.0 * np.pi * float(N) * pos[2] / alpha) / alpha, 
                               2.0 * np.pi * float(N) * R * np.cos( 2.0 * np.pi * float(N) * pos[2] / alpha) / alpha,
                               H / alpha ]) for pos in cylPosList ]
        
        #normalize the tangent vector
        THatList = [ NewT/np.linalg.norm(NewT) for NewT in TList ]
        
        # fIO.saveXYZ(THatList, 'Ca', 'THatList.xyz')        
        
        # compute the Normal vector (equivalent to the zero point vector)
        NList = [ np.array([ np.cos(2.0 * np.pi * float(N) * pos[2]/alpha), 
                             np.sin(2.0 * np.pi * float(N) * pos[2]/alpha),
                             0]) for pos in cylPosList ]
        # normalize the normal vector
        NHatList = [ NewN/np.linalg.norm(NewN) for NewN in NList ]
        
        # fIO.saveXYZ(NHatList, 'Li', 'NHatList.xyz')
        
        # compute the Binormal vector  
        BHatList = [ np.cross(newT, newN) for newT, newN in zip(THatList, NHatList) ]
    
        # fIO.saveXYZ(BHatList, 'Be', 'BHatList.xyz')
    
        # In the TNB frame, apply the relevant relative coords rho and phi (pos[0] and pos[1]. the T coord is always zero. 
        return [ basePoint + pos[0] * np.cos(pos[1]) * NHat + pos[0] * np.sin(pos[1]) * BHat for pos, NHat, BHat, basePoint in zip(cylPosList, NHatList, BHatList, basePointList)]
    
    
    
    def parseEnvelopeList(self, envelopeList):
        self.envelopeSummary = {}
        for envelopeString in envelopeList:
            envelope = envelopeString.split()
            self.envelopeSummary[envelope[0].lower()]=tuple([float(param) for param in envelope[1:]])
    
    def generateBuildingBlockXYZ(self):
        # TODO overload this function for generating the xyz positions
        return [np.array([i * .5, rnd.uniform(-3,3), rnd.uniform(-1,1)]) for i in range(0, self.numPoints)] 

    def generateBuildingBlockNames(self):
        # TODO overload this function for generating the names
        return ['C'] * self.numPoints 

    def generateBuildingBlockDirector(self):
        # TODO overload this for generating a principal direction for the backbone
        return np.array([1,0,0]) 

    def generateBuildingBlockConnectors(self):
        # TODO overload this function for defining the connectors of the backbone
        return np.array( [ [2 , 1, 0],
                           [len(self.buildingBlockXYZ) - 3, 
                            len(self.buildingBlockXYZ) - 2,
                            len(self.buildingBlockXYZ) - 1]] )

    def generateBuildingBlockRefPoint(self):
        # TODO overload this function defining a principal reference point in the backbone 
        return cart.getCentreOfMass(self.buildingBlockXYZ)

    def getBuildingBlock(self):
        # access function that returns a building block based on the backbone
        return BuildingBlock(self.buildingBlockXYZ, 
                             self.blockNames, 
                             self.blockConnectors, 
                             self.blockRefPoint,
                             self.blockDirectorHat)

    def checkEnvelope(self, xyzVals, ignorePTA=False, pointsToAvoid=[]):
        # returns a list of indices of the points which violate the envelope 
        # and pointsToAvoid list.
        indicesOutside = []
        for n, pos in enumerate(xyzVals):
            if not self.checkPointInBounds(pos, ignorePTA=ignorePTA, pointsToAvoid=pointsToAvoid):
                indicesOutside.append(n)
        return indicesOutside

    def checkPointInBounds(self, pos, ignorePTA=False, pointsToAvoid=[]):
        # checks the external conditions that we wish to enforce on the points
        # These are either pointsToAvoid or envelopes (regions of space that are out of bounds).
        # Can supply an auxiliary set of pointsToAvoid which overides the pointsToAvoid that are
        # the member variable. Can shoose to ignore pointsToAvoid, whether overwritten or not. 
        # pos is in xyz coords
        
        # assume the pos is inBounds to begin with.
        inBounds = True

        if ignorePTA==False:
            
            # if points to avoid are not supplied as a parameter (i.e. have the default value of an empty list)
            # then set the local pointsToAvoid Variable to point to self.pointsToAvoid
            if not pointsToAvoid:
                pointsToAvoid=self.pointsToAvoid
                
            # check pos against whatever points to avoid we have determined above 
            for pta in pointsToAvoid:
                if np.linalg.norm(pta - pos) < self.minDist:
                    inBounds = False
                    if self.verbose==1:
                        print("pointsToAvoid Violation")
                    break
        
        # Now perform envelop processing.
        
        # We check the envelope summary dictionary created at run time against every known type of envelope.
        # If the checked envelope type is in the supplied envelopeSummary dictionary then 
        # we enforce that envelope. The value of the dictionary is a list of the parameters 
        # associated with that envelope type. Keywords are always all lower case. 

        # Always check to see if we are still inBounds before performing what can be a lengthy test.
        if inBounds:
            # in this envelope any point outside the sphere centered at blockRefPoint is out of bounds
            try:
                envelopeParams = self.envelopeSummary['outersphere']
                # sphere centered at refPoint with input radius
                inBounds = True
                if np.linalg.norm( pos - self.blockRefPoint) > envelopeParams[0]:
                    if self.verbose==1:
                        print("outerSphere Violation")
                    inBounds=False
            except KeyError:
                pass

        if inBounds:        
            # in this envelope any point inside the inner sphere centered at blockRefPoint is out of bounds
            try:
                envelopeParams = self.envelopeSummary['innersphere']
                # sphere is centered at blockRefPoint.
                inBounds = True
                if np.linalg.norm(pos - self.blockRefPoint) < envelopeParams[0]:
                    if self.verbose==1:
                        print("innerSphere Violation")
                    inBounds=False
            except KeyError:
                pass

        if inBounds:        
            # ensure that the wedge is inside the spherical polar region
            # defined by theta1, theta2, phi1 and phi2
            # -pi/2 < theta < pi/2
            # -pi < phi < pi  
            try:
                envelopeParams = self.envelopeSummary['sphericalWedge']
                
                theta1 = envelopeParams[0]
                theta2 = envelopeParams[1]
                phi1 = envelopeParams[2]
                phi2 = envelopeParams[3]
                inside = envelopeParams[4] # determines where point should be inside or outside defined region
                thetaMin = min(theta1, theta2)
                thetaMax = max(theta1, theta2)
                phiMin = min(phi1, phi2)
                phiMax = max(phi1, phi2)                
                xyzSpherical = coords.XYZ2SphericalPolar(pos)
                
                inBounds = True
                if inside >= 0:
                    if xyzSpherical[1] < thetaMin or xyzSpherical[1] > thetaMax:
                        if self.verbose==1:
                            print("Spherical Wedge Violation")
                        inBounds=False
                    if xyzSpherical[2] < phiMin or xyzSpherical[2] > phiMax:
                        if self.verbose==1:
                            print("Spherical Wedge Violation")
                        inBounds=False

                if inside < 0:
                    if xyzSpherical[1] > thetaMin and xyzSpherical[1] < thetaMax:
                        if self.verbose==1:
                            print("Spherical Wedge Violation")
                        inBounds=False
                    if xyzSpherical[2] > phiMin and xyzSpherical[2] < phiMax:
                        if self.verbose==1:
                            print("Spherical Wedge Violation")
                        inBounds=False
                    
            except KeyError:
                pass
        
        if inBounds:
            # Envelope is frustum with radius1 at z1 and radius2 at z2.
            try:
                envelopeParams = self.envelopeSummary['frustum']
                
                Z1 = envelopeParams[0]
                Radius1 = envelopeParams[1]
                Z2 = envelopeParams[2]
                Radius2 = envelopeParams[3]
                inBounds = coords.checkPointInFrustum(pos, Z1, Radius1, Z2, Radius2, self.verbose)

            except KeyError:
                pass
        
        if inBounds:
            # makes sure that the pos is in the +ve Z half space
            try:
                envelopeParams = self.envelopeSummary['halfspace']
                if pos[2]>envelopeParams[0]:
                    inBounds=False
                    if self.verbose==1:
                        print("halfspace Violation")                    
            except KeyError:
                pass

        if inBounds:
            # in this envelope any point outside the sphere is out of bounds
            try:
                envelopeParams = self.envelopeSummary['betasphere']
                # the params define a sphere of radius envelopeParams[1]
                # which is centered on z-axis a distance envelopeParams[0] from the origin.
                # Any pos inside that sphere is accepted.
                inBounds = True
                if np.linalg.norm(pos - np.array([0.0, 0.0, envelopeParams[0]])) > envelopeParams[1]:
                    inBounds=False
            except KeyError:
                pass

        if inBounds:
            # in this envelope any point outside the spherecylinder is out of bounds
            try:
                envelopeParams = self.envelopeSummary['spherocylinder']
                # the params define two points and a radius
                # the two points define the line segment swept by a sphere of radius minDist that form a spherocylinder of radius minDist.
                p1 = np.array([envelopeParams[0], envelopeParams[1], envelopeParams[2]])
                q1 = np.array([envelopeParams[3], envelopeParams[4], envelopeParams[5]])
                minDist = envelopeParams[6]
                inBounds = True
                # abuse the well writen closestApproachTwoLineSegmentsSquared function which can deal with degenerate points.
                if cart.closestApproachTwoLineSegmentsSquared(p1, q1, pos, pos) > minDist**2:
                    inBounds=False
            except KeyError:
                pass

        
        if inBounds:
            # in this envelope any point outside the sphere is out of bounds
            try:
                envelopeParams = self.envelopeSummary['endsphere']
                # the params define a sphere of radius envelopeParams[1]
                # which is centered on z-axis a distance envelopeParams[0] from the origin.
                # Any pos inside that sphere is rejected.
                inBounds = True
                if np.linalg.norm(pos - np.array([0.0, 0.0, envelopeParams[0]])) < envelopeParams[1]:
                    inBounds=False
            except KeyError:
                pass

        if inBounds:
            # in this envelope any point outside the a cuboid range is out of bounds
            try:
                envelopeParams = self.envelopeSummary['cuboid']
                # the six params define xmin, xmax, ymin, ymax, zmin and zmax
                inBounds = True
                if pos[0] < envelopeParams[0]: # xmin
                    inBounds=False
                if pos[0] > envelopeParams[1]: # xmax
                    inBounds=False                    
                if pos[1] < envelopeParams[2]: # ymin
                    inBounds=False
                if pos[1] > envelopeParams[3]: # ymax
                    inBounds=False                   
                if pos[2] < envelopeParams[4]: # zmin
                    inBounds=False
                if pos[2] > envelopeParams[5]: # zmax
                    inBounds=False                   
            
            except KeyError:
                pass

        return inBounds
    
    def visualiseEnvelopeFunc(self, N, X, Y, Z, filename):
        
        posList = [ np.array([rnd.uniform(-X, X), rnd.uniform(-Y, Y), rnd.uniform(-Z, Z) ]) for _ in range(N) ]

        posOut = []
        n = 0
        for pos in posList:
            if n % 10000==0:
                print("Visualising Envelope:", str(n), "out of", str(len(posList)))
            if self.checkPointInBounds(pos, ignorePTA=True):
                posOut.append(pos)
            n += 1
        fIO.saveXYZ(posOut, 'Ne', filename)
        
        return
    
    def checkAllPointsInBounds(self, xyzVals):
        ''' Creates a list of indices of all the points in XYZ that failed the checkPointInBounds tests.'''
        # set up output arrays        
        pointGood = True
        badPoints = []
        
        for n, pointToTest in enumerate(xyzVals):
            # perform spatial test
            pointGood = self.checkPointInBounds(pointToTest)
    
            # if the point is bad then add the index to the array of bad points        
            if not pointGood:
                badPoints.append(n)

        return badPoints
    
            
if __name__ == "__main__":
    
    # specify filename
    filename = "buildingBlockGeneratorExample.txt"
    
    # create the backbone generator object using static file parameters
    BBG = BuildingBlockGenerator(filename)

    # generate building block realtime parameters
    numPos = 30
    startPos = np.array([0.0, 0.0, 0.0])
    director = np.array([0.0, 0.0, 1.0])
    rotation = 45 * np.pi/180
    minDist = 1.0
    
    testBuildBlock = BBG.generateBuildingBlock(numPos, minDist)
    testBuildBlock.transformBBToLabFrame(director, startPos, rotation) 
    testBuildBlock .exportBBK(fIO.fileRootFromInfile(filename, 'txt'))
    
    print("building block done")