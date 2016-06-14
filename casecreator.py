
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

## Class with static information about the case
class caseinfo:
    isoX = 0.0
    isoY = 0.0
    R = 1.0
    ## Number of beamlets in the fan
    N = 64
    interleaf = 0.6
    ## Source to axis distance calibration in cms
    SAD = 80
    ## Original fan
    genFan2D = None
    def __init__(self, x = 0.0, y = 0.0, radio = 1.0):
        self.isoX = x
        self.isoY = y
        self.R = radio
        self.genFan2D = np.transpose(np.matrix([[self.interleaf * (i - self.N/2 + 1/2), self.SAD] for i in range(0, self.N)]))

## Class that uses a data pair and implements some geographical operations. Depth of the voxel given beam.
class geoloc:
    def __init__(self, x = 0.0, y = 0.0):
        self.x = x
        self.y = y
        self.depth = None
    ## This function calculates the distance from my geographical location to the center of a beamlet
    def distToBeamC(self, xBeamC, yBeamC , zBeamC ):
        d = np.sqrt(np.sum((self.x - xBeamC )**2 + (self.y - yBeamC )**2))
        return(d)
    ## Function to calculate distance from this point to isocenter
    def distToIsoC(self):
        d = np.sqrt(np.sum((self.x - thiscase.isoX)**2 + (self.y - thiscase.isoY)**2))
        return(d)

    ## This function finds whether a point lies INSIDE the line SEGMENT between the beamlet and the voxel or not.
    def isinterior(self, xinterp, xBeamC):
        ininterior = False
        if (np.min(xBeamC, self.x) <= xinterp and xinterp <= np.max(xBeamC, self.x)):
            ininterior = True
        return(ininterior)

    ## Find the depth of this voxel inside the body
    def depthBeamC(self, xBeamC, yBeamC, zBeamC, dBeamlettoIso):
        ## To understand the methodology look at http://mathworld.wolfram.com/Circle-LineIntersection.html
        x2 = xBeamC
        y2 = yBeamC
        x1 = self.x
        y1 = self.y
        dx = x2 - x1
        dy = y2 - y1
        dr = np.sqrt(dx**2 + dy**2)
        D = x1 * y2 - x2 * y1
        ## There are two point of intersection
        xinterp = (D * dy + np.sign(dy) * dx * np.sqrt(thiscase.R**2 * dr**2 - D**2))/(dr**2)
        xinterm = (D * dy - np.sign(dy) * dx * np.sqrt(thiscase.R**2 * dr**2 - D**2))/(dr**2)
        yinterp = (-D * dx + np.abs(dy) * np.sqrt(thiscase.R ** 2 * dr ** 2 - D ** 2)) / (dr ** 2)
        yinterm = (-D * dx - np.abs(dy) * np.sqrt(thiscase.R ** 2 * dr ** 2 - D ** 2)) / (dr ** 2)
        ## Check which one of the intersection points lies in the segment
        if (self.isinterior(xinterp, xBeamC)):
            intX = xinterp
            intY = yinterp
        else:
            intX = xinterm
            intY = yinterm
        ## Check that indeed you did everything right
        assert(np.min(xBeamC, self.x) <= intX and intX <= np.max(xBeamC, self.x))
        ## Use the point of intersection to calculate the depth of the voxel
        self.depth = np.sqrt((intX - self.x)**2 + (intY - self.y)**2)

## Abstract class that implements a volume of interest with common location and radius. Parent of OAR and TARGET
class VOI:
    numVOIs = 0
    __metaclass__ = ABCMeta
    def __init__(self, x = 0.0, y = 0.0, r = 0.0):
        ## X location of center
        self.xcenter = x
        ## Y location of ceter
        self.ycenter = y
        ## Radius of the Volume of Interest. All of them are circumferences
        self.radius = r
        self.isTarget = None
        self.isinside = self.isContained()
        self.VOIID = VOI.numVOIs
        VOI.numVOIs = VOI.numVOIs + 1
    ## This method finds whether the attribute is viable, given its center and radius and given the center and radius
    ## of the original body that contains it
    def isContained(self):
        isv = True
        ## Find radius from center of VOI to center of structure
        distcenter = np.sqrt((self.xcenter - thiscase.isoX)**2 + (self.ycenter - thiscase.isoY)**2)
        if distcenter + self.radius > thiscase.R:
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
    def __init__(self, x = 0.0, y = 0.0, r = 0.0):
        ## Boolean. Is this a target structure?
        self.isTarget = True
        super(OAR, self).__init__(x, y, r)
        ## Assign an ID to each of the different OARs
        self.OARID = OAR.numOARS
        OAR.numOARS = OAR.numOARS + 1

    def printVOI(self):
        print('OAR with center (', self.xcenter, ', ', self.ycenter, '); and radius ', self.radius)

class TARGET(VOI):
    numTARGETS = 0
    def __init__(self, x = 0.0, y = 0.0, r = 0.0):
        ## Boolean. Is this a target structure?
        self.isTarget = False
        super(TARGET, self).__init__(x, y, r)
        ## Assign an ID to each of the different targets
        self.TARGETID = TARGET.numTARGETS
        TARGET.numTARGETS = TARGET.numTARGETS + 1
    def printVOI(self):
        print('Target with center (', self.xcenter, ', ', self.ycenter, '); and radius ', self.radius)

## The next class defines a control point; in particular, the location of all beamlets
class ControlPoint:
    def __init__(self, ctrlAngle):
        angleDegs = ctrlAngle
        angleRads = (2 * np.pi * angleDegs)/360
        rotMat = np.matrix([[np.cos(angleRads), -np.sin(angleRads)], [np.sin(angleRads), np.cos(angleRads)]])
        thisFan = rotMat * thiscase.genFan2D

class voxel:
    def __init__(self, vc, OARS, TARGETS):
        self.x = vc[0]
        self.y = vc[1]
        self.belongsToVOI = False
        self.inStructureID = None
        ## Run this code for all OARs and TARGETs, preference to targets
        for voi in OARS + TARGETS:
            if voi.isInThisVOI(self.x, self.y):
                self.belongsToVOI = True
                self.inStructureID = voi.VOIID
## This function takes a list of numeric values as arguments and produces the list of D matrices
## The box has isocenter on position thiscase.x
# Inputs:
# anglelist = [list of numeric]. These are the control points
# numhozv = Number of horizontal voxel divisions
# numverv = Number of vertical voxel divisions
# xgeoloc = X Geographic location of upper-right side of box
# ygeoloc = Y Geographic location of upper-right side of box
# radius  = Radius of the body
# OARList = List of OAR centers and radiuses
# TARGETList = List of Target centers and radiuses
# Outputs:
# listofD = List of D matrix objects
# upper right corner is (X,Y), lower left corner is (-X, -Y)
def createDosetoPoints(anglelist, numhozv, numverv, xgeoloc, ygeoloc, radius, OARS, TARGETS):
    ## Generate 360 control points
    cps = [ControlPoint(i) for i in anglelist]
    ## Create voxels
    voxelhoz = np.arange(-xgeoloc, xgeoloc, 2 * xgeoloc/numhozv) + xgeoloc/numhozv
    voxelvec = np.arange(-ygeoloc, ygeoloc, 2 * ygeoloc/numverv) + ygeoloc/numverv
    ## Create cartesian product to find voxel centers
    voxelcenters = itertools.product(voxelhoz, voxelvec)
    ## Limit the list only to those voxels that are included in the body and assign a organ to them
    allvoxels = [voxel(voxelcenter, OARS, TARGETS) for voxelcenter in voxelcenters]
    ## Filter only those voxels that belong in any VOI
    voxels = [vxinvoi for vxinvoi in allvoxels if vxinvoi.belongsToVOI]
    [print(i.inStructureID) for i in voxels]
    allvoxels = None # Free some memory
    return(1)

## Implementation part that should be separated later
def plotstructure(OARlist, TARGETlist, xgeo, ygeo):
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
    fig.savefig('plotcase.png')

bodyradius = 7.0
thiscase = caseinfo(0.0, 0.0, bodyradius)
## Create 3 OARs around the body
OARlist = []
OARlist.append(OAR(5.0, 0.0, 1.0))
OARlist.append(OAR(3.0, 2.0, 2.5))
OARlist.append(OAR(-3.0, 2.0, 1.3))
TARGETlist = []
TARGETlist.append(TARGET(0.0, 0.0, 2.0))
TARGETlist.append(TARGET(-4.0, -4.0, 1.2))

## Initialize an instance. This is not necessary but it is good practice.
xgeoloc = bodyradius
ygeoloc = bodyradius
plotstructure(OARlist, TARGETlist, xgeoloc, ygeoloc)
numhozv = 20
numverv = 20
xgeoloc = bodyradius
ygeoloc = bodyradius
radius = bodyradius
anglelist = [i * 360 / 51 for i in range(0, 51)]

createDosetoPoints(anglelist, numhozv, numverv, xgeoloc, ygeoloc, radius, OARlist, TARGETlist)