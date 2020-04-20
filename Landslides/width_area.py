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

# Run this in the root directory above the six landslide experiment folders.

_save = True

def recursive_glob(rootdir='.', pattern='*'):
	"""Search recursively for files matching a specified pattern.
	
	Adapted from http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
	"""
	matches = []
	for root, dirnames, filenames in os.walk(rootdir):
	  for filename in fnmatch.filter(filenames, pattern):
		  matches.append(os.path.join(root, filename))
	return matches

runtimes = pd.read_csv('runtime.txt', sep='\t')

filenames = recursive_glob(rootdir='.', pattern='Landslides.txt')

frames = []
for filename in filenames:
    frames.append(pd.read_csv(filename))

for frame in frames:
    print np.max(frame['runtime [s]'])/3600.

data = pd.concat(frames)


##############
# WIDTH-AREA #
##############

def plFit(x, a, b):
    return a*x**b

def linFit(x, a):
    return a*x

_valmask = np.isfinite(data['width [m]']) & np.isfinite(data['area [m2]'])
popt, pcov = curve_fit(plFit, data['width [m]'][_valmask],
                              data['area [m2]'][_valmask])
#popt, pcov = curve_fit(linFit, data['width [m]'][_valmask],
#                               data['area [m2]'][_valmask])

plt.figure(figsize=(6.4,4.8*.85)) # new default is 6.4, 4.8
plt.loglog(data['width [m]'], data['area [m2]'], marker='o', color = '0', alpha=.2, linestyle='None')
_x = np.logspace( np.log10(plt.xlim()[0]), np.log10(plt.xlim()[1]), 100 )
plt.loglog(_x, plFit(_x, *popt), '-', color='.75', linewidth=2)
#plt.loglog(_x, linFit(_x, *popt), '-', color='.75', linewidth=2)
plt.xlabel('Landslide width $w$ [m]', fontsize=16)
plt.ylabel('Landslide area $A$ [m$^2$]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.text(6E-2, 2E-4, '$A$ = '+'%.2f' %popt[0]+'$w^{'+'%.2f' %popt[1]+'}$', fontsize=16)
plt.tight_layout()

if _save:
    plt.savefig('wA_log.svg')
    plt.savefig('wA_log.pdf')
    plt.savefig('wA_log.png')
else:
    plt.show()
    
