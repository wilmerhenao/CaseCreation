
#!/home/wilmer/anaconda3/bin/python

__author__ = 'wilmer'
try:
    import mkl
    have_mkl = True
    print("Running with MKL Acceleration")
except ImportError:
    have_mkl = False
    print("Running with normal backends")

import numpy as np
import scipy as sp

class VOI:

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

