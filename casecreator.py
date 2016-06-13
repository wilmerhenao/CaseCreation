
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
    isoZ = 0.0
    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        isoX = x
        isoY = y
        isoZ = z
## Initialize an instance. This is not necessary but it is good practice.
thiscase = caseinfo(0.0, 0.0, 0.0)

## Class that uses a data triplet and implements some geographical operations
class geoloc:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z
    ## This function calculates the distance from my geographical location to the center of a beamlet
    def distToBeamC(self, xBeamC, yBeamC , zBeamC ):
        d = np.sqrt(np.sum((self.x - xBeamC )^2 + (self.y - yBeamC )^2 + (self.z - zBeamC )^2))
        return(d)
    ## Function to calculate distance from this point to isocenter
    def distToIsoC(self, xisoC, yisoC , zisoC ):
        d = np.sqrt(np.sum((self.x - xisoC )^2 + (self.y - yisoC )^2 + (self.z - zisoC )^2))
        return(d)

    def depthBeamC(self, xBeamC, yBeamC, zBeamC, dBeamlettoIso):
        ## Use the cosine law to calculate the angle formed by the line from isocenter to this point and the line from
        # beamlet center to this point.
        c = dBeamlettoIso
        a = distToIsoC()

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

