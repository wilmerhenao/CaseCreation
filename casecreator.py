
#!/home/wilmer/anaconda3/bin/python

__author__ = 'wilmer'
try:
    import mkl
    have_mkl = True
    print("Running with MKL Acceleration")
except ImportError:
    have_mkl = False
    print("Running with normal backends")

from abc import ABCMeta
import numpy as np
import scipy as sp

## Class with static information about the case
class caseinfo:
    isoX = 0.0
    isoY = 0.0
    R = 1.0
    def __init__(self, x = 0.0, y = 0.0, radio):
        self.isoX = x
        self.isoY = y
        self.R = radio
## Initialize an instance. This is not necessary but it is good practice.
thiscase = caseinfo(0.0, 0.0, 1.0)

## Class that uses a data pair and implements some geographical operations. Depth of the voxel given beam.
class geoloc:
    def __init__(self, x = 0.0, y = 0.0):
        self.x = x
        self.y = y
        self.depth = None
    ## This function calculates the distance from my geographical location to the center of a beamlet
    def distToBeamC(self, xBeamC, yBeamC , zBeamC ):
        d = np.sqrt(np.sum((self.x - xBeamC )^2 + (self.y - yBeamC )^2))
        return(d)
    ## Function to calculate distance from this point to isocenter
    def distToIsoC(self):
        d = np.sqrt(np.sum((self.x - thiscase.isoX)^2 + (self.y - thiscase.isoY)^2))
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
        dr = np.sqrt(dx^2 + dy^2)
        D = x1 * y2 - x2 * y1
        ## There are two point of intersection
        xinterp = (D * dy + np.sign(dy) * dx * np.sqrt(caseinfo.R^2 * dr^2 - D^2))/(dr^2)
        xinterm = (D * dy - np.sign(dy) * dx * np.sqrt(caseinfo.R^2 * dr^2 - D^2))/(dr^2)
        yinterp = (-D * dx + np.abs(dy) * np.sqrt(caseinfo.R ^ 2 * dr ^ 2 - D ^ 2)) / (dr ^ 2)
        yinterm = (-D * dx - np.abs(dy) * np.sqrt(caseinfo.R ^ 2 * dr ^ 2 - D ^ 2)) / (dr ^ 2)
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
        self.depth = np.sqrt((intX - self.x)^2 + (intY - self.y)^2)


## Abstract class that implements a volume of interest with common location and radius. Parent of OAR and TARGET
class VOI:
    __metaclass__ = ABCMeta
    def __init__(self, x = 0.0, y = 0.0, r = 0.0):
        ## X location of center
        self.xcenter = x
        ## Y location of ceter
        self.ycenter = y
        ## Radius of the Volume of Interest. All of them are circumferences
        self.radius = r
        self.isTarget = None
    ## This method finds whether the attribute is viable, given its center and radius and given the center and radius
    ## of the original body that contains it
    def isContained(self, rbody):
        isv = True
        ## Find radius from center of VOI to center of structure
        distcenter = np.sqrt(self.xcenter * self.xcenter + self.ycenter * self.ycenter)
        if distcenter + self.radius > rbody:
            isv = False
        return (isv)
    # Abstract method to be implemented by classes deriving from here
    @abstractmethod
    def printVOI(self):
        pass

class OAR(VOI):
    def __init__(self, x = 0.0, y = 0.0, r = 0.0):
        ## Boolean. Is this a target structure?
        self.isTarget = True
        super(OAR, self).__init__(x, y, r)

    def printVOI(self):
        print('OAR with center (', self.xcenter, ', ', self.ycenter, '); and radius ', self.radius)

class TARGET(VOI):
    def __init__(self, x = 0.0, y = 0.0, r = 0.0):
        ## Boolean. Is this a target structure?
        self.isTarget = False
        super(TARGET, self).__init__(x, y, r)

    def printVOI(self):
        print('Target with center (', self.xcenter, ', ', self.ycenter, '); and radius ', self.radius)

class voxel:

## This function takes a list of numeric values as arguments and produces the list of D matrices
## The box has isocenter on position 0,0.
# Inputs:
# anglelist = [list of numeric]
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
def createDosetoPoints(anglelist, numhozv, numverv, xgeoloc, ygeoloc, radius):

    return(listofD)

