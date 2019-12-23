import shapefile
import numpy as np
import matplotlib.pyplot as plt
import glob
import fnmatch
import os
import gdal
from scipy.optimize import curve_fit
from scipy.special import gamma
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import pandas as pd
plt.ion()

folders = sorted(glob.glob('*mmhr'))
nls = []

# Count the number of landslides in each experiment
for folder in folders:
    shapefiles = sorted(glob.glob(folder+'/*.shp'))
    n = 0
    for _shapefile in shapefiles:
        sf = shapefile.Reader(_shapefile)
        shapes = sf.shapes()
        n += len(shapes)
    nls.append(n)

out = np.array(np.hstack((folders, nls))).transpose()

