# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.abspath('..'))
import casecreator

## How big is the phantom going to be? remember that the aperture has 64 beamlets with a space of 0.6cm in between
bodyradius = 20.0
thiscase = casecreator.caseinfo(0.0, 0.0, bodyradius)
## Create OARs and TARGETs around the body
OARlist = []
OARlist.append(casecreator.OAR(thiscase, 12.5, 0.0, 7.0))
OARlist.append(casecreator.OAR(thiscase, 8.0, 12.0, 3.5))
OARlist.append(casecreator.OAR(thiscase, -3.0, 2.0, 1.3))
OARlist.append(casecreator.OAR(thiscase, 1.1, -11.4, 4.0))
OARlist.append(casecreator.OAR(thiscase, -11.1, 11.4, 2.0))
TARGETlist = []
TARGETlist.append(casecreator.TARGET(thiscase, 0.0, 0.0, 2.2))
TARGETlist.append(casecreator.TARGET(thiscase, -4.0, -4.0, 1.2))
TARGETlist.append(casecreator.TARGET(thiscase, 5.0, 5.0, 2.7))

xgeoloc = bodyradius
ygeoloc = bodyradius
numhozv = 10
numverv = 10
radius = bodyradius
anglelist = [i * 360 / 51 for i in range(0, 51)]

[myDs, voxels] = casecreator.createDosetoPoints(thiscase, anglelist, numhozv, numverv, xgeoloc, ygeoloc, radius, OARlist, TARGETlist)
casecreator.plotstructure(thiscase, OARlist, TARGETlist, xgeoloc, ygeoloc, voxels)
print(myDs)
