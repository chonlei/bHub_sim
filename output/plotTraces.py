import re
import numpy as np
import matplotlib.pyplot as plt
import sys

# setup data
try:
    fileName = sys.argv[1]
except Exception:
    fileName = "Ca_WT_5phubs_9s_sGJ_100by501.dat"
shape = re.findall('.*_(\w+)by(\w+)\.dat',fileName)[0]
shapeX = (int(shape[0]),int(shape[1]))
varName = r"[Ca]$_i$ [mM]"
title = "short simulation homogeneous GJ @ mouse 40-3"
saveName = fileName[:-3]+'png'

# load data
X = np.memmap(fileName, dtype='float64', mode='r', shape=shapeX)
time = np.arange(shapeX[1])*10.0  # account dt

# hubList (if have)
hubList = [62, 77, 27, 4, 85, 57, 66, 37, 48, 45]#[64, 74, 34, 18, 11, 61, 98, 44, 94, 47]#1025, 570, 1294, 81, 169, 659, 890, 1622, 1486, 1250, 247, 59, 595, 1526, 546, 1008, 1629, 1748, 923, 872, 742, 635, 920, 977, 1333, 867, 1438, 434, 524, 1053, 1235, 420, 718, 1042, 835, 399, 775, 1275, 1573, 148, 407, 365, 1240, 515, 523, 452, 391, 1743, 1049, 753, 1463, 388, 1502, 448, 166, 1718, 687, 1090, 1536, 254, 1219, 1490, 683, 0, 584, 1701, 1747, 11, 1424, 1350, 377, 1062, 1031, 1026, 1266, 1652, 367, 1181, 669, 550, 643, 1137, 1349, 1411, 141, 1247, 95, 1082]

# main plot
aveX = np.zeros(X[0].shape)
for i in xrange(len(X)):
    if i not in hubList:
        plt.plot(time, X[i], 'k', alpha=0.3)
        aveX += X[i]
plt.plot(time, X[i], 'k', alpha=0.3, label='non-hub')
for i in hubList:
    plt.plot(time, X[i], 'r')
    aveX += X[i]
plt.plot(time, X[i], 'r', label='hub')
aveX = aveX/float(shapeX[0])
plt.plot(time, aveX, 'g', linewidth=3.0, label='average all')

# plotting setup
plt.xlabel("t [ms]")
plt.ylabel(varName)
plt.title(title)
plt.savefig(saveName)

