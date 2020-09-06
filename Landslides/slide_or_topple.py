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
C_slide_triangle = (rho * g * W**2 * (np.sin(alpha) - mu*np.cos(alpha))) \
                            / (2 * (H**2 * W**2)**.5)
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

fig = plt.figure(figsize=(6, 5.4))
plt.loglog(C_slide[triangular_failures == False], TS_topple[triangular_failures == False], 'ks', alpha=.1, label='Trapezoidal slide')
plt.loglog(C_slide[triangular_failures], TS_topple[triangular_failures], 'k^', alpha=.1, label='Triangular slide')
plt.xlim(5E0,2E4)
plt.ylim(5E0,2E4)
plt.loglog(_x, _y, 'k-', linewidth=2)
plt.xlabel('Cohesion for slide failure [Pa]', fontsize=16)
plt.ylabel('Tensile strength for topple failure [Pa]', fontsize=16)
plt.legend()
plt.tight_layout()

if _save:
    plt.savefig('MassWasting_PhaseDiagram.svg')

