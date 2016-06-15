# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.abspath('..'))

from casecreator import *

bodyradius = 20.0
thiscase = caseinfo(0.0, 0.0, bodyradius)
## Create 3 OARs around the body
OARlist = []
OARlist.append(OAR(thiscase, 12.5, 0.0, 7.0))
OARlist.append(OAR(thiscase, 8.0, 12.0, 3.5))
OARlist.append(OAR(thiscase, -3.0, 2.0, 1.3))
OARlist.append(OAR(thiscase, 1.1, -11.4, 4.0))
OARlist.append(OAR(thiscase, -11.1, 11.4, 2.0))
TARGETlist = []
TARGETlist.append(TARGET(thiscase, 0.0, 0.0, 2.2))
TARGETlist.append(TARGET(thiscase, -4.0, -4.0, 1.2))
TARGETlist.append(TARGET(thiscase, 5.0, 5.0, 2.7))

## Initialize an instance. This is not necessary but it is good practice.
xgeoloc = bodyradius
ygeoloc = bodyradius
numhozv = 20
numverv = 20
xgeoloc = bodyradius
ygeoloc = bodyradius
radius = bodyradius
anglelist = [i * 360 / 51 for i in range(0, 51)]

[myDs, voxels] = createDosetoPoints(thiscase, anglelist, numhozv, numverv, xgeoloc, ygeoloc, radius, OARlist, TARGETlist)
plotstructure(thiscase, OARlist, TARGETlist, xgeoloc, ygeoloc, voxels)
print(myDs)
