# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.abspath('..'))
import casecreator
import pickle
from scipy import sparse
import numpy as np
from scipy import array

## How big is the phantom going to be? remember that the aperture has 64 beamlets with a space of 0.6cm in between
bodyradius = 20.0
thiscase = casecreator.caseinfo(0.0, 0.0, bodyradius)
## Create OARs and TARGETs around the body

TARGETlist = []
TARGETlist.append(casecreator.TARGET(thiscase, 0.0, 0.0, 2.2))
TARGETlist.append(casecreator.TARGET(thiscase, -4.0, -4.0, 1.2))
TARGETlist.append(casecreator.TARGET(thiscase, 5.0, 5.0, 2.7))
OARlist = []
OARlist.append(casecreator.OAR(thiscase, 12.5, 0.0, 7.0))
OARlist.append(casecreator.OAR(thiscase, 8.0, 12.0, 3.5))
OARlist.append(casecreator.OAR(thiscase, -3.0, 2.0, 1.3))
OARlist.append(casecreator.OAR(thiscase, 1.1, -11.4, 4.0))
OARlist.append(casecreator.OAR(thiscase, -11.1, 11.4, 2.0))

xgeoloc = bodyradius
ygeoloc = bodyradius

numhozv = 30
numverv = 30
radius = bodyradius
anglelist = [i * 360 / 51 for i in range(0, 51)]

[myDs, voxels] = casecreator.createDosetoPoints(thiscase, anglelist, numhozv, numverv, xgeoloc, ygeoloc, OARlist, TARGETlist)
casecreator.plotstructure(thiscase, OARlist, TARGETlist, xgeoloc, ygeoloc, voxels)

# Print them in pickle format
f = open('voxels.pckl', 'wb')
pickle.dump(voxels, f)
f.close()
f = open('Ds.pckl', 'wb')
pickle.dump(myDs, f)
f.close()

# Print them in txt format
masklist = [voxel.inStructureID for voxel in voxels]

blist = []
jlist = []
dlist = []

for D in myDs:
    [b, j, d] = sparse.find(D)
    blist.extend(b)
    jlist.extend(j)
    dlist.extend(d)

casecreator.savevector('C:/Users/S170452/PycharmProjects/Tomotherapy-Without-Pulse/data/myBixels_out.bin', blist, np.int32)
casecreator.savevector('C:/Users/S170452/PycharmProjects/Tomotherapy-Without-Pulse/data/myVoxels_out.bin', jlist, np.int32)
casecreator.savevector('C:/Users/S170452/PycharmProjects/Tomotherapy-Without-Pulse/data/myDijs_out.bin', dlist, np.float32)
casecreator.savevector('C:/Users/S170452/PycharmProjects/Tomotherapy-Without-Pulse/data/myoptmask.bin', masklist, np.int32)
