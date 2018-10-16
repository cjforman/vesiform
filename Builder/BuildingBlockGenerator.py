import sys
from Utilities.keyProc import keyProc
import Utilities.cartesian as cart
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO

import numpy as np
import random as rnd
from BuildingBlock import BuildingBlock 

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
            print "Critical Parameters are undefined for BuildingBlockGenerator"
            sys.exit()        

    def generateBuildingBlock(self, numPoints, minDist, showBlockDirector=False, visualiseEnvelope=(0, 20, 'envelope.xyz'), envelopeList=['None'], pointsToAvoid=[], defaultBlockRefPoint=None):
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
            print "Visualising envelope"
            self.visualiseEnvelope(visualiseEnvelope[0], visualiseEnvelope[1], visualiseEnvelope[1], visualiseEnvelope[1], visualiseEnvelope[2])

        # from the points to avoid list only remember those that would pass the specified envelope test
        self.pointsToAvoid =[]
        self.pointsToAvoid = [ pointsToAvoid[index] for index in range(len(pointsToAvoid)) if not index in self.checkAllPointsInBounds(pointsToAvoid) ]

        # save the pointsToAvoid if we're dumping files
        if self.dumpInterimFiles==1:
            fIO.saveXYZ(self.pointsToAvoid, 'Li', "pointsToAvoid.xyz")
        
        # generate the xyz positions
        self.buildingBlockXYZ = self.generateBuildingBlockXYZ()

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

    def checkEnvelope(self, xyzVals):
        # returns a list of indices of the points which violate the envelope 
        # and pointsToAvoid list.
        indicesOutside = []
        for n, pos in enumerate(xyzVals):
            if not self.checkPointInBounds(pos):
                indicesOutside.append(n)
        return indicesOutside

    def checkPointInBounds(self, pos, ignorePTA=False):
        # checks the external conditions that we wish to enforce on the points
        # These are either points to avoid or envelopes (regions of space that are out of bounds). 
        # pos is in xyz coords
        
        # assume the pos is inBounds to begin with.
        inBounds = True

        if ignorePTA==False:
            # check pos against points to avoid 
            for pta in self.pointsToAvoid:
                if np.linalg.norm(pta - pos) < self.minDist:
                    inBounds = False
                    if self.verbose==1:
                        print "pointsToAvoid Violation"
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
                        print "outerSphere Violation"
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
                        print "innerSphere Violation"
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
                            print "Spherical Wedge Violation"
                        inBounds=False
                    if xyzSpherical[2] < phiMin or xyzSpherical[2] > phiMax:
                        if self.verbose==1:
                            print "Spherical Wedge Violation"
                        inBounds=False

                if inside < 0:
                    if xyzSpherical[1] > thetaMin and xyzSpherical[1] < thetaMax:
                        if self.verbose==1:
                            print "Spherical Wedge Violation"
                        inBounds=False
                    if xyzSpherical[2] > phiMin and xyzSpherical[2] < phiMax:
                        if self.verbose==1:
                            print "Spherical Wedge Violation"
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
                
                ZMin = min(Z1, Z2)
                ZMax = max(Z1, Z2)
                
                # principal axis of envelope is aligned with Z axis. 
                # Compute radial distance from z axis and z height
                zComponent = pos[2] 
                r = np.sqrt(pos[0]**2 + pos[1]**2)
        
                if np.abs(zComponent - 0.0) < 1e-10:
                    zComponent = 0.0
        
                inBounds = True
                if  zComponent > ZMax or zComponent < ZMin:
                    inBounds = False
                    if self.verbose==1:
                        print "frustum Z Violation"

                
                if r > coords.FrustumRadiusAtGivenZ(zComponent, Z1, Z2, Radius1, Radius2):
                    inBounds = False
                    if self.verbose==1:
                        print "frustum R Violation" 
            except KeyError:
                pass
        
        if inBounds:
            # makes sure that the pos is in the +ve Z half space
            try:
                envelopeParams = self.envelopeSummary['halfspace']
                if pos[2]>envelopeParams[0]:
                    inBounds=False
                    if self.verbose==1:
                        print "halfspace Violation"                    
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

        return inBounds
    
    def visualiseEnvelope(self, N, X, Y, Z, filename):
        
        posList = [ np.array([rnd.uniform(-X, X), rnd.uniform(-Y, Y), rnd.uniform(-Z, Z) ]) for _ in range(N) ]

        posOut = []
        n = 0
        for pos in posList:
            if n % 1000==0:
                print "Visualising Envelope:", str(n), "out of", str(len(posList))
            if self.checkPointInBounds(pos, ignorePTA=True):
                posOut.append(pos)
            n += 1
        fIO.saveXYZ(posOut, 'Ne', filename)
        
        return
    
    def checkAllPointsInBounds(self, xyzVals):
        ''' Creates a list of all the points in XYZ that failed the spatial tests.'''
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
    
    print "building block done"