#!/home/wilmer/anaconda3/bin/python

__author__ = 'wilmer'
try:
    import mkl
    have_mkl = True
    print("Running with MKL Acceleration")
except ImportError:
    have_mkl = False
    print("Running with normal backends")

from abc import ABCMeta, abstractmethod
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pylab
import itertools
import math
from scipy.sparse import csr_matrix, lil_matrix

## Class with static information about the case
class caseinfo:
    ## Center of the body (which also will contain a tumour) Coord. X
    isoX = 0.0
    ## Center of the body (which also will contain a tumour) Coord. X
    isoY = 0.0
    ## Default radius of the body. Will be overriden
    R = 20.0
    ## Number of beamlets in the fan
    N = 64
    ## Interbeamlet distance in the fan. Which amounts to 6mms.
    interleaf = 0.6
    ## Source to axis distance calibration in cms
    SAD = 80
    ## Original fan with N positions (64)
    genFan2D = None
    ## Constructor
    def __init__(self, x = 0.0, y = 0.0, radio = 20.0):
        self.isoX = x
        self.isoY = y
        self.R = radio
        self.genFan2D = np.matrix([[self.interleaf * (i - self.N/2 + 1/2), self.SAD] for i in range(0, self.N)]).transpose()

## Class that uses a data pair and implements some geographical operations. Depth of the voxel given beam.
class voxelbeamletpair:
    plotsomething = 0
    ## Constructor Function
    def __init__(self, v):
        ## Voxel center x coordinate
        self.x = v.x
        ## Voxel Center y coordinate
        self.y = v.y
        ## Variable containing the depth of the voxel in the direction from the beamlet or how much the beam travels.
        self.depth = None
    ## This function calculates the distance from my geographical location to the center of a beamlet (not used)
    def distToBeamC(self, xBeamC, yBeamC):
        d = np.sqrt(np.sum((self.x - xBeamC )**2 + (self.y - yBeamC )**2))
        return(d)
    ## Function to calculate distance from this point to isocenter (not used)
    def distToIsoC(self, thiscase):
        d = np.sqrt(np.sum((self.x - thiscase.isoX)**2 + (self.y - thiscase.isoY)**2))
        return(d)
    ## This function finds whether a point lies INSIDE the line SEGMENT between the beamlet and the voxel or not.
    def isinterior(self, xinterp, xBeamC):
        ininterior = False
        if (min(xBeamC, self.x) <= xinterp and xinterp <= max(xBeamC, self.x)):
            ininterior = True
        return(ininterior)
    ## Find the depth of this voxel inside the body. The depth will be the only factor used in order to calculate
    # accumulated dose.
    def depthBeamC(self, xBeamC, yBeamC, R):
        ## To understand the methodology look at http://mathworld.wolfram.com/Circle-LineIntersection.html
        # First I initialize some variables
        x2 = xBeamC
        y2 = yBeamC
        x1 = self.x
        y1 = self.y
        dx = x2 - x1
        dy = y2 - y1
        dr = np.sqrt(dx**2 + dy**2)
        D = x1 * y2 - x2 * y1
        # There are two point of intersection
        xinterp = (D * dy + np.sign(dy) * dx * np.sqrt(R**2 * dr**2 - D**2))/(dr**2)
        xinterm = (D * dy - np.sign(dy) * dx * np.sqrt(R**2 * dr**2 - D**2))/(dr**2)
        yinterp = (-D * dx + np.abs(dy) * np.sqrt(R ** 2 * dr ** 2 - D ** 2)) / (dr ** 2)
        yinterm = (-D * dx - np.abs(dy) * np.sqrt(R ** 2 * dr ** 2 - D ** 2)) / (dr ** 2)
        # Check which one of the intersection points lies in the segment. Only 1 coordinate necessary.
        if (self.isinterior(xinterp, xBeamC)):
            intX = xinterp
            intY = yinterp
        else:
            intX = xinterm
            intY = yinterm
        # Check that indeed you did everything right
        assert(min(xBeamC, self.x) <= intX and intX <= max(xBeamC, self.x))
        ## Use the point of intersection of line and circle to calculate the depth of the voxel
        self.depth = np.sqrt((intX - self.x)**2 + (intY - self.y)**2)
        if 2810 == self.plotsomething:
            print('plotting something')
            circlemain = plt.Circle((0, 0), 20, color='blue', fill=False)
            fig2 = plt.gcf()
            fig2.gca().add_artist(circlemain)
            pylab.xlim([min(xBeamC, self.x) - 2, max(xBeamC, self.x) + 2])
            pylab.ylim([min(yBeamC, self.y) - 2, max(yBeamC, self.y) + 2])
            fig2.suptitle('Ray Plot. Beamlet source located at ('+ str(xBeamC) + ', ' + str(yBeamC) + ')')
            plt.plot([self.x, xBeamC], [self.y, yBeamC], 'ro', mew = 10)
            # Plot the line segment from the voxel center to the beamlet
            # Calculate directionvector
            vx = xBeamC - self.x
            vy = yBeamC - self.y
            nm = np.sqrt(vx**2 + vy**2)
            vx, vy = self.depth * vx/nm, self.depth * vy/nm
            plt.plot((self.x, self.x + vx), (self.y, self.y + vy), 'r', lw = 1)
            plt.plot((self.x + vx, xBeamC), (self.y + vy, yBeamC), 'g--', lw = 1)
            fig2.savefig('voxelandray.png')
        voxelbeamletpair.plotsomething = voxelbeamletpair.plotsomething + 1

## Abstract class that implements a volume of interest with common location and radius. Parent of OAR and TARGET
class VOI:
    ## Static counter of Volumes of Interest
    numVOIs = 0
    ## Definition necessary for an abstract class
    __metaclass__ = ABCMeta
    ## Constructor function
    def __init__(self, thiscase, x = 0.0, y = 0.0, r = 0.0):
        ## Object with general information about the case
        self.tc = thiscase
        ## X location of center
        self.xcenter = x
        ## Y location of ceter
        self.ycenter = y
        ## Radius of the Volume of Interest. All of them are circumferences
        self.radius = r
        ## Boolean that determines whether this is a target or not (in case of False, it is an OAR)
        self.isTarget = None
        ## Is this region contained inside the body?
        self.isinside = self.isContained()
        ## Unique ID for each particular Volume of Interest
        self.VOIID = VOI.numVOIs
        VOI.numVOIs = VOI.numVOIs + 1
    ## This method finds whether the attribute is viable, given its center and radius and given the center and radius
    ## of the original body that contains it
    def isContained(self):
        isv = True
        ## Find radius from center of VOI to center of structure
        distcenter = np.sqrt((self.xcenter - self.tc.isoX)**2 + (self.ycenter - self.tc.isoY)**2)
        if distcenter + self.radius > self.tc.R:
            isv = False
        return (isv)
    ## This method takes a location in space and returns whether this location exists in this VOI or not
    def isInThisVOI(self, x, y):
        isinit = False
        distcenter = np.sqrt((self.xcenter - x) ** 2 + (self.ycenter - y) ** 2)
        if distcenter <= self.radius:
            isinit = True
        return(isinit)
    # Abstract method to be implemented by classes deriving from here
    @abstractmethod
    def printVOI(self):
        pass

class OAR(VOI):
    numOARS = 0
    ## Constructor function that also calls the constructor of VOI. Notice that a VOI object is instantiated first, and
    # then an OAR object is instantiated later
    def __init__(self, thiscase, x = 0.0, y = 0.0, r = 0.0):
        ## Boolean. Is this a target structure?
        self.isTarget = True
        super(OAR, self).__init__(thiscase, x, y, r)
        ## Assign an ID to each of the different OARs
        self.OARID = OAR.numOARS
        OAR.numOARS = OAR.numOARS + 1
    ## Print the characteristics of this structure on screen
    def printVOI(self):
        print('OAR with center (', self.xcenter, ', ', self.ycenter, '); and radius ', self.radius)

class TARGET(VOI):
    numTARGETS = 0
    ## Constructor function that also calls the constructor of VOI. Notice that a VOI object is instantiated first, and
    # then a TARGET object is instantiated later
    def __init__(self, thiscase,  x = 0.0, y = 0.0, r = 0.0):
        ## Boolean. Is this a target structure?
        self.isTarget = False
        super(TARGET, self).__init__(thiscase, x, y, r)
        ## Assign an ID to each of the different targets
        self.TARGETID = TARGET.numTARGETS
        TARGET.numTARGETS = TARGET.numTARGETS + 1
    ## Print the characteristics of this structure on screen
    def printVOI(self):
        print('Target with center (', self.xcenter, ', ', self.ycenter, '); and radius ', self.radius)

## The next class defines a control point; in particular, the location of all beamlets, notice that angle or zero is
# located in the part above, and the gantry moves counterclockwise
class ControlPoint:
    ## Constructor Function
    def __init__(self, ctrlAngle, thiscase):
        ## Global data structure of the case
        self.tc = thiscase
        ## This control point angle in degrees
        self.angleDegs = ctrlAngle
        ## This control point angle in radians
        self.angleRads = (2 * np.pi * self.angleDegs)/360
        rotMat = np.matrix([[np.cos(self.angleRads), -np.sin(self.angleRads)], [np.sin(self.angleRads), np.cos(self.angleRads)]])
        ## This fan beam. Notice that it gets rotated the number of radians necessary
        self.thisFan = rotMat * thiscase.genFan2D
        ## Find the unit vector that points towards the isocenter
        self.UnitVector = (np.sin(self.angleRads), -np.cos(self.angleRads))
        ## Find the unit vector that is perpendicular to the original unit vector
        self.normaltoUnit = (-self.UnitVector[1], self.UnitVector[0])
    ## Find normal distances to each of the beamlet array centers. This function returns an array of distances from which
    # I will choose those that are small enough to correspond to beams that affect the radiation to the supplied voxel.
    def findNDist(self, x, y):
        distances = []
        for i in range(0, self.tc.N):
            beamlet = (self.thisFan[0,i], self.thisFan[1,i])
            #print('i, beamlet, x, y: ', i, beamlet, x, y)
            #print('unit vector', self.UnitVector)
            vecpos = x - beamlet[0], y - beamlet[1]
            #print('distance:', math.fabs(vecpos[1] * self.normaltoUnit[1] + vecpos[0] * self.normaltoUnit[0]))
            distances.append(math.fabs(vecpos[1] * self.normaltoUnit[1] + vecpos[0] * self.normaltoUnit[0]))
        return(distances)

## This class defines not only the x,y position of a voxel, but also assigns to it a unique ID and maps a structure to it.
class voxel:
    ## Static variable that serves as a counter of how many voxels are actively being used at any time.
    numVOXELS = 0
    def __init__(self, vc, OARS, TARGETS):
        ## Indicates a unique ID for each of the voxels
        self.voxelID = voxel.numVOXELS
        ## x location of voxel center
        self.x = vc[0]
        ## y location of voxel center
        self.y = vc[1]
        ## Does this voxel belong to ANY VOI?
        self.belongsToVOI = False
        ## ID of the VOI to which this voxel belongs to. There is a hierarchy that depends on the order of the VOIS with
        # targets taking precedence over OARs.
        self.inStructureID = None
        ## Run this code for all OARs and TARGETs, preference to targets
        for voi in OARS + TARGETS:
            if voi.isInThisVOI(self.x, self.y):
                self.belongsToVOI = True
                self.inStructureID = voi.VOIID

## This function calculates the total dose given a certain depth. It is short but it was created independently so that
# it is easily modified later
def calcDose(depth):
    return(np.exp(-0.04 * depth))

## This is the main function takes a list of numeric values as arguments and produces the list of D matrices
## The box has isocenter on position thiscase.x
# Inputs:
# anglelist = [list of numeric]. These are the control points
# numhozv = Number of horizontal voxel divisions
# numverv = Number of vertical voxel divisions
# xgeoloc = X Geographic location of upper-right side of box
# ygeoloc = Y Geographic location of upper-right side of box
# OARList = List of OAR centers and radiuses
# TARGETList = List of Target centers and radiuses
# Outputs:
# listofD = List of D matrix objects
# upper right corner is (X,Y), lower left corner is (-X, -Y)
def createDosetoPoints(thiscase, anglelist, numhozv, numverv, xgeoloc, ygeoloc, OARS, TARGETS):
    ## Generate 360 control point
    cps = [ControlPoint(i, thiscase) for i in anglelist]
    ## Create voxels
    voxelhoz = np.arange(-xgeoloc, xgeoloc, 2 * xgeoloc/numhozv) + xgeoloc/numhozv
    voxelvec = np.arange(-ygeoloc, ygeoloc, 2 * ygeoloc/numverv) + ygeoloc/numverv
    ## Create cartesian product to find voxel centers
    voxelcenters = itertools.product(voxelhoz, voxelvec)
    ## Limit the list only to those voxels that are included in the body and assign a organ to them
    allvoxels = [voxel(voxelcenter, OARS, TARGETS) for voxelcenter in voxelcenters]
    voxel.numVOXELS = 0 # Fix this value at zero because I will have to recount in the next line
    ## Filter only those voxels that belong in any VOI. This is the order that will be preserved
    voxels = [voxel((vxinvoi.x, vxinvoi.y), OARS, TARGETS) for vxinvoi in allvoxels if vxinvoi.belongsToVOI]
    allvoxels = None # Free some memory
    Dlist = []
    ## Create a matrix for each of the control points.
    for cp in cps:
        D = lil_matrix((len(voxels), caseinfo.N), dtype = np.float)
        #print(cp.thisFan, cp.thisFan[0,25])
        for v in voxels:
            bsst = findvoxbeamlets(v.x, v.y, cp)
            if not bsst:
                continue
            else:
                bs = bsst[0]
                for blet in bs:
                    vbpair = voxelbeamletpair(v)
                    vbpair.depthBeamC(cp.thisFan[0, bs[0]], cp.thisFan[1, bs[0]], thiscase.R)
                    D[v.voxelID, bs[0]] = calcDose(vbpair.depth)
        Dlist.append(D)
    return(Dlist, voxels)

## This function is here to find the beamlets associated with a particular voxel at a certain control point
## You should be able to change it in case the team requires a different function
def findvoxbeamlets(x, y, cp):
    dists = cp.findNDist(x, y)
    return([(i, dists[i]) for i in range(0, len(dists)) if dists[i] < 0.6/2])

## Implementation part that should be separated later
def plotstructure(thiscase, OARlist, TARGETlist, xgeo, ygeo, voxels):
    ## This function plots the case to make sure that everything is understood
    numOARS = len(OARlist)
    numTARGETS = len(TARGETlist)
    # Plot the outside circle
    circlemain = plt.Circle((thiscase.isoX, thiscase.isoY), thiscase.R, color = 'blue', fill = False)
    fig = plt.gcf()
    fig.gca().add_artist(circlemain)
    pylab.xlim([-xgeo, xgeo])
    pylab.ylim([-ygeo, ygeo])
    for i in range(0, numOARS):
        circle = plt.Circle((OARlist[i].xcenter, OARlist[i].ycenter), OARlist[i].radius, color = 'g', fill = False)
        fig.gca().add_artist(circle)
    for i in range(0, numTARGETS):
        circle = plt.Circle((TARGETlist[i].xcenter, TARGETlist[i].ycenter), TARGETlist[i].radius, color = 'r', fill = False)
        fig.gca().add_artist(circle)
    fig.suptitle('Case Plot')
    ## Plot points in the voxel structure
    for v in voxels:
        plt.plot(v.x, v.y, 'ro')
    fig.savefig('plotcase.png')
