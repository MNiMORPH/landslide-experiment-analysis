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

_save = False

rho = 2650 * .65 + .1 * 1000
g = 9.8
mu = 0.6

h = np.arange(2E-2, 40E-2, 1E-2)
w = np.arange(1E-2, 40E-2, 1E-2)
H, W = np.meshgrid(h, w)
alpha = 76 * np.pi/180.

data = pd.read_csv('LandslidesBeaulieu_AllDataCals.csv')
H = data['depth [m]']
W = data['width [m]']

#H = .4
#W = .1

# Cohesion = shear strength
C_slide_trapezoid = rho * g * (W - H/(2*np.tan(alpha))) * np.sin(alpha) \
                           * (np.sin(alpha) - mu*np.cos(alpha))
C_slide_triangle = rho * g * W * H * (H - mu*W) / (H**2 + W**2)
C_slide = C_slide_trapezoid.copy()
triangular_failures = H > (W*np.tan(alpha))
print( np.sum(triangular_failures) / float(len(triangular_failures)) )
C_slide[triangular_failures] = C_slide_triangle


# Tensile strength
#TS_topple = rho * g * W**2 / (W**2+H**2)**.5 * (W**2+H**2)**.5
TS_topple = rho * g * W**2 / H


# Relationship between tensile strength and shear cohesion
#relative_tensile_strength = C_slide / .6

# Which mechanism will be seen?
# Topple if TS_topple > C_slide/.6
_x = np.linspace(5E0, 2E4, 10)
_y = _x/.6 # tensile vs. shear failure: cohesion = y-int; tensile str = x-int

# Fraction to fail at each criterion
topple_preferred = TS_topple > C_slide/.6
f_topple_predicted = np.sum(TS_topple > C_slide/.6) / float(len(data))
print(f_topple_predicted)

# Triangular failures among non-topples
triangular_failure_notopple = triangular_failures[topple_preferred == False]
print( np.sum(triangular_failure_notopple) / float(len(triangular_failure_notopple)) )

# Ratio between topple and sliding criteria
print( np.mean(TS_topple / C_slide) )
# Just 1.2x stress at topple

# Cohesion estimate at failure for predicted topples vs. predicted slides
print("MEDIAN:")
print("Topple 'cohesion':", np.nanmedian(C_slide[topple_preferred]))
print("Slide cohesion:", np.nanmedian(C_slide[topple_preferred == False]))
print("MEAN:")
print("Topple 'cohesion':", np.nanmean(C_slide[topple_preferred]))
print("Slide cohesion:", np.nanmean(C_slide[topple_preferred == False]))

# Which ones occur within the predicted range of cohesions?
fig = plt.figure(figsize=(6, 5.4))
plt.loglog(C_slide[triangular_failures == False], TS_topple[triangular_failures == False], 'ks', alpha=.1, label='Trapezoidal slide')
plt.loglog(C_slide[triangular_failures], TS_topple[triangular_failures], 'k^', alpha=.1, label='Triangular slide')
plt.xlim(5E0,2E4)
plt.ylim(5E0,2E4)
plt.loglog(_x, _y, 'k-', linewidth=2)
plt.xlabel('Driving shear stress minus fritional resistance\n(i.e., stress contributing against cohesion) [Pa]', fontsize=16)
plt.ylabel('Stress contributing towards topple [Pa]', fontsize=16)
plt.legend()
plt.tight_layout()

if _save:
    plt.savefig('MassWasting_PhaseDiagram.svg')

# Where do the predicted topples lie?
n, bins = np.histogram(C_slide, bins=np.arange(0, 3201, 100))
mag = (bins[:-1] + bins[1:])/2.
freq = n/float(np.sum(n))/np.diff(bins)

nT, binsT = np.histogram(C_slide[topple_preferred], bins=np.arange(0, 3201, 100))
magT = (bins[:-1] + bins[1:])/2.
freqT = nT/float(np.sum(n))/np.diff(bins)

nS, binsS = np.histogram(C_slide[topple_preferred == False], bins=np.arange(0, 3201, 100))
magS = (bins[:-1] + bins[1:])/2.
freqS = nS/float(np.sum(n))/np.diff(bins)

fig = plt.figure() #figsize=(6, 4.8))
plt.plot( bins[:-1], freqT / freq, 'k-', linewidth=2 )
plt.ylabel('Fraction of failures as topples', fontsize=16)
plt.xlabel('Equivalent landslide cohesion at failure [Pa]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(0, 2300)
plt.ylim(-.01,1.01)
plt.tight_layout()

if _save:
    plt.savefig('PredictedToppleFrequency.svg')

"""
plt.figure(figsize=(6.4,4.8*.85)) # new default is 6.4, 4.8
plt.axvspan(250, 750, fc='None', hatch='//', label='Literature-based cohesion\nof unsaturated sand')
plt.bar(bins[:-1], freq, width=np.diff(bins), ec="k", fc='.8', align="edge")
plt.bar(bins[:-1], freqS, width=np.diff(bins), ec="k", fc='.3', align="edge")
ax = plt.gca()
plt.ylabel('Probability density [Pa$^{-1}$, m$^{-1}$]', fontsize=16)
plt.xlabel('Landslide cohesion at failure [Pa]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(0, ax.get_xlim()[-1])
plt.legend()
plt.twinx()
plt.ylim(np.array(ax.get_ylim()) * float(np.sum(n)) * np.mean(np.diff(bins)))
plt.ylabel('Number of landslides', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()
"""


