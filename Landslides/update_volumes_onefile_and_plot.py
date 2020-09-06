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

# Within experiment data folder
runtimes = pd.read_csv('runtime.txt', sep='\t')
filenames = sorted(glob.glob(pathname='*/Landslides.txt', recursive=True))

# Concatenate the frames and
# identify the experiment to which the frame belongs
frames = []
for filename in filenames:
    base_level_fall_rate = int( filename.split('_')[-1].split('mmhr')[0] )
    frame = pd.read_csv(filename)
    frame.insert(1, 'base-level fall rate [mm/hr]', base_level_fall_rate)
    frame['base-level fall rate [mm/hr]'] = base_level_fall_rate
    frames.append(frame)
    print( np.max(frame['runtime [s]'])/3600. )

data = pd.concat(frames)

data.to_csv('LandslidesBeaulieu_AllDataCals.csv', index=False)

# Update volume
scarp_angle = 76. * np.pi/180. # degrees to radians

# Trapezoid
Vtrap = data['depth [m]'] * data['area [m2]'] - data['length [m]'] * data['depth [m]']**2 / (2 * np.tan(72*scarp_angle))
# Triangle -- steeper angle required?
Vtri = data['depth [m]'] * data['area [m2]'] / 2.
# Larger of the two
V = np.nanmax([Vtrap, Vtri], axis=0)
# Update column
data['volume [m3]'] = V

###############
# AREA-VOLUME #
###############

def plFit(x, a, b):
    return a*x**b

popt, pcov = curve_fit(plFit, data['area [m2]'], data['volume [m3]'])

plt.figure(figsize=(6.4,4.8*.85)) # new default is 6.4, 4.8
plt.loglog(data['area [m2]'], data['volume [m3]'], marker='o', color = '0', alpha=.2, linestyle='None', markeredgewidth=0)
_x = np.logspace( np.log10(plt.xlim()[0]), np.log10(plt.xlim()[1]), 100 )
plt.loglog(_x, plFit(_x, *popt), '-', color='.75', linewidth=2)
plt.ylabel('Landslide volume $V$ [m$^3$]', fontsize=16)
plt.xlabel('Landslide area $A$ [m$^2$]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.text(1E-2, 4E-6, '$V$ = '+'%.2f' %popt[0]+'$A^{'+'%.2f' %popt[1]+'}$', fontsize=16)
plt.tight_layout()

if _save:
    plt.savefig('AV_log_traptri.svg')

#####################
# COHESIVE STRENGTH #
#####################

rho = 2650 * .65 + .1 * 1000
g = 9.8
coeff_friction = 0.6
#CohesionBlockFall_simple = rho * g * data['width [m]']

Area_XS_trap = data['depth [m]'] * (data['width [m]'] - data['depth [m]'] / (2 * np.tan(scarp_angle)) )
Area_XS_tri = data['depth [m]'] * data['width [m]'] / 2.
Area_XS = np.nanmax([Area_XS_trap, Area_XS_tri], axis=0)
scarp_length = data['depth [m]'] / np.sin(scarp_angle)
#width_function = Area_XS / data['depth [m]']
Cohesion = rho*g*Area_XS * ( np.sin(scarp_angle) - coeff_friction * np.cos(scarp_angle)) / scarp_length
n, bins = np.histogram(Cohesion, bins=np.arange(0, 3201, 100))
mag = (bins[:-1] + bins[1:])/2.
freq = n/float(np.sum(n))/np.diff(bins)
#plt.bar(bins[:-1], freq, width=np.diff(bins), ec='.8', fc='.8', align="edge", alpha=.2)

# Just for width
n_width, bins_width = np.histogram(data['width [m]'][np.isfinite(data['width [m]'])], bins=len(bins)-1)
mag_width = (bins[:-1] + bins[1:])/2.
freq_width = n/float(np.sum(n))/np.diff(bins)

# Shear stress vs. probability density
plt.figure(figsize=(6.4,4.8*.85)) # new default is 6.4, 4.8
#plt.plot(mag, freq, 'ko-')
plt.axvspan(250, 750, fc='None', hatch='//', label='Literature-based cohesion\nof unsaturated sand')
plt.bar(bins[:-1], freq, width=np.diff(bins), ec="k", fc='.8', align="edge")
ax = plt.gca()
plt.ylabel('Probability density [Pa$^{-1}$, m$^{-1}$]', fontsize=16)
plt.xlabel('Landslide cohesion at failure [Pa]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(0, 2500) # ax.get_xlim()[-1])
plt.legend()
# Equivalent landslide width
#plt.twiny()
#plt.xlim(np.array(ax.get_xlim()) / rho / g * 1000)
#plt.xlabel('Trapezoidal landslide width [mm]', fontsize=16)
#plt.bar(bins_width[:-1]*1000., freq_width, width=np.diff(bins_width*1000.), ec="k", fc='b', alpha=.5, align="edge")
#ax = plt.gca()
#plt.tick_params(axis='both', which='major', labelsize=12)
# Equivalent number of landslides
plt.twinx()
plt.ylim(np.array(ax.get_ylim()) * float(np.sum(n)) * np.mean(np.diff(bins)))
plt.ylabel('Number of landslides', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
# Done!
plt.tight_layout()

if _save:
    plt.savefig('StressFrequency.svg')

# Now, area--frequency with colors based on being within the expected values
fig = plt.figure()
ax = plt.subplot(1,1,1)
ax.loglog( data['area [m2]'], Cohesion, 'ko', alpha=.2 )
ax.axhspan(ymin=250, ymax=650, color='.8')


_inrange = (Cohesion >= 200) & (Cohesion <= 750)
_outrange = (1-_inrange).astype(bool)
cdla = 400E-6

n, bins = np.histogram(data['area [m2]'][_inrange], bins=np.logspace(np.log10(0.00012), np.log10(0.66), 26))
mag_inrange = (bins[:-1] + bins[1:])/2.
freq_inrange = n/float(len(data['area [m2]']))/np.diff(bins)

n, bins = np.histogram(data['area [m2]'][_outrange], bins=np.logspace(np.log10(0.00012), np.log10(0.66), 26))
mag_outrange = (bins[:-1] + bins[1:])/2.
freq_outrange = n/float(len(data['area [m2]']))/np.diff(bins)

# Linear plot

plt.figure()
plt.axvspan(0, 400E-6, fc='None', hatch='\\\\\\', label='Failure below limit of consistent detection')
plt.fill_between(mag_inrange, 0, freq_inrange, color='.3', label='Failure stress within expected range of cohesion')
plt.fill_between(mag_outrange, freq_inrange, freq_inrange+freq_outrange, color='.6', label='Failure stress outside expected range of cohesion')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1E-4, 1E-1)
plt.ylabel('Probability density $f$ [m$^{-2}$]', fontsize=16)
plt.xlabel('Landslide area $A$ [m$^2$]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend(loc='lower right')
#plt.twinx()
#plt.ylim(np.array(ax.get_ylim()) * float(len(data['area [m2]'])) * np.mean(np.diff(bins)))
#plt.ylabel('Number of landslides', fontsize=16)
#plt.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()
if _save:
    plt.savefig('Af_cohesion.svg')

# LOOK UP SHEFFIELD THESIS 2015
# Sticking with Lu2009 and using internal friction angle + tensile strength to obtain cohesion


#plt.figure()
#plt.

if not _save:
    plt.show()
else:
    plt.close()
    
