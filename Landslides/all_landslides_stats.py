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


##################
# AREA-FREQUENCY #
##################

def invGammaFit(x, rho, a, s):
    return 1./(a * gamma(rho)) * (a/(x - s)) ** (rho + 1.) * np.exp(-a/(x-s))
def lognormFit(x, s, mu):
    return (1./(s*x*(2*np.pi)**.5)) * np.exp( - (np.log(x)-mu)**2 / (2.*s**2) )
def normFit(x, mu, s2):
    return 1./(2*np.pi*s2)**.5 * np.exp(- (x-mu)**2 / (2*s2) )

# Linear bin spacing
#n, bins = np.histogram(data['area'], bins=np.linspace(cdla, 0.66, 512))
#mag = (bins[:-1] + bins[1:])/2.
#freq = n/float(np.sum(n))#/np.diff(bins)

cdla = 400E-6

n, bins = np.histogram(data['area [m2]'], bins=np.logspace(np.log10(0.00012), np.log10(0.66), 26))
mag = (bins[:-1] + bins[1:])/2.
freq = n/float(np.sum(n))/np.diff(bins)

popt, pcov = curve_fit(invGammaFit, mag[mag>cdla], freq[mag>cdla], p0 = (1.4, 4E-3, -1E-4))
popt1, pcov1 = curve_fit(lognormFit, mag[mag>cdla], freq[mag>cdla])

# Linear plot
plt.figure()
plt.plot(mag[mag <= cdla], freq[mag <= cdla], marker='o', color='.5', linestyle = 'None')
plt.plot(mag[mag >= cdla], freq[mag >= cdla], marker='o', color = '0', linestyle = 'None')
plt.xlim(-0.002, 0.04)
_x = np.logspace( np.log10( np.max((plt.xlim()[0], 1E-6))), np.log10(plt.xlim()[1]), 100 )
plt.plot(_x, invGammaFit(_x, *popt), 'k-')
#plt.semilogx(_x, lognormFit(_x, *popt1))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.ylabel('$f$ [m$^{-2}$]', fontsize=24)
plt.xlabel('$A$ [m$^{2}$]', fontsize=24)
plt.tight_layout()
if _save:
    plt.savefig('Af_linear.svg')

# Log plot
plt.figure()
plt.loglog(mag[mag <= cdla], freq[mag <= cdla], marker='o', color='.5', linestyle = 'None')
plt.loglog(mag[mag >= cdla], freq[mag >= cdla], marker='o', color = '0', linestyle = 'None')
plt.xlim(1E-4, 1E-1)
_x = np.logspace( np.log10(plt.xlim()[0]), np.log10(plt.xlim()[1]), 100 )
plt.loglog(_x, invGammaFit(_x, *popt), 'k-')
plt.ylabel('Probability density $f$ [m$^{-2}$]', fontsize=16)
plt.xlabel('Landslide area $A$ [m$^2$]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()
if _save:
    plt.savefig('Af_log.svg')


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
    plt.savefig('AV_log.svg')


#############
# LAG TIMES #
#############

wts = np.array(data['wait time [s]'])
wts_blc = wts[data['runtime [s]'] < (4*3600)] # base-level constant
wts = wts[np.isfinite(wts)]
wts_blc = wts_blc[np.isfinite(wts_blc)]
# Looks like they are the same! Great.

#plt.figure()
#plt.hist(wts, bins=128)
#plt.hist(wts_blc, bins=64)

# FRACTION OF TOTAL OR PROBABILITY DENSITY?

wts_n, wts_bins = np.histogram(wts, bins=np.arange(0,2400, 24))
wts_bin_centers = (wts_bins[:-1] + wts_bins[1:])/2.
wts_p = wts_n/float(len(wts))#/np.diff(wts_bins) # divide by number of items and bin widths

# Let's fit something here
def plaw(x, a, b):
    return a*x**b
    
popt1, pcov1 = curve_fit(plaw, wts_bin_centers[wts_bin_centers>300], wts_p[wts_bin_centers>300])
print popt1[1]

popt2, pcov2 = curve_fit(plaw, wts_bin_centers[wts_bin_centers<200], wts_p[wts_bin_centers<200])
print popt2[1]


plt.figure(figsize=(6.4,4.8*.85))
plt.semilogy(wts_bin_centers, wts_p, 'o', color='.4')
plt.xlabel('Lag time between landslides [s]', fontsize=16)
plt.ylabel('Fraction of total landslides', fontsize=16)
#_x = np.linspace(100, 1E4, 100)
#plt.plot(_x, plaw(_x, popt1[0], popt1[1]), 'k-', linewidth=2)
#_x = np.linspace(10, 400, 100)
#plt.plot(_x, plaw(_x, popt2[0], popt2[1]), 'k-', linewidth=2)
#plt.xlim(10,4000)
#plt.ylim(5E-4, 5E-1)
plt.tight_layout()

# How does this compare to random integers?
# ~ 60 hrs of total experiment time.... different for different ones.
# 1360 total landslides
# put it all together
"""
N = 6
wts_rls_all = []
for i in range(6):
    rls = sorted(np.random.randint(low=0, high=np.sum(runtimes['t_s'])/N, size=(len(wts)/N,)))
    wts_rls = np.diff(rls)
    wts_rls_all += list(wts_rls)
wts_n, wts_bins = np.histogram(wts_rls_all, bins=np.arange(0,2400, 24))
wts_bin_centers = (wts_bins[:-1] + wts_bins[1:])/2.
wts_p = wts_n/float(len(wts_rls_all)) # divide by bin widths and number of items
plt.semilogy(wts_bin_centers, wts_p, 'o', color='green')
"""

"""
wts_n, wts_bins = np.histogram(wts, bins=64)
wts_bin_centers = (wts_bins[:-1] + wts_bins[1:])/2.
wts_p = wts_n/float(len(wts)) # divide by bin widths and number of items
plt.semilogx(wts_bin_centers, wts_p, 'o', color='.7')
"""

poisson_r = len(wts) / float(np.sum(runtimes['t_s'])) # event density [# events / interval]
# I want to find k events to maximize Poisson
t_lag = np.arange(0, 1201, 1)
pPoisson = np.exp(-poisson_r * t_lag) * poisson_r * np.mean(np.diff(wts_bins))
plt.semilogy(t_lag, pPoisson, 'k-', linewidth=3)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(5E-4, 2E-1)
plt.xlim(0, 2250)

# number of landslides
old_ylim = np.array(plt.ylim())
plt.twinx()
ax = plt.gca()
ax.set_yscale('log')
plt.ylim(old_ylim * len(wts))
plt.ylabel('Number of landslides', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()

if _save:
    plt.savefig('LagTimes.svg')


# PROBABILITY DENSITY INSTEAD
wts_p = wts_n/float(len(wts))/np.diff(wts_bins) # divide by number of items and bin widths

plt.figure(figsize=(6.4,4.8*.85))
plt.semilogy(wts_bin_centers, wts_p, 'o', color='.4')
plt.xlabel('Lag time between landslides [s]', fontsize=16)
plt.ylabel('Probability density [s$^{-1}$]', fontsize=16)
plt.tight_layout()

pPoisson = np.exp(-poisson_r * t_lag) * poisson_r
plt.semilogy(t_lag, pPoisson, 'k-', linewidth=3)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(0, 2250)
plt.ylim(2E-5, 1E-2)

# number of landslides
old_ylim = np.array(plt.ylim())
plt.twinx()
ax = plt.gca()
ax.set_yscale('log')
plt.ylim(old_ylim * len(wts) * np.mean(np.diff(wts_bins)))
plt.ylabel('Number of landslides', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()

if _save:
    plt.savefig('LagTimes_probDensity.svg')




#####################
# COHESIVE STRENGTH #
#####################

rho = 2650 * .65 + .1 * 1000
g = 9.8
CohesionBlockFall_simple = rho * g * data['width [m]']
n, bins = np.histogram(CohesionBlockFall_simple, bins=np.arange(0, 3201, 100))
mag = (bins[:-1] + bins[1:])/2.
freq = n/float(np.sum(n))/np.diff(bins)
#plt.bar(bins[:-1], freq, width=np.diff(bins), ec='.8', fc='.8', align="edge", alpha=.2)

# Shear stress vs. probability density
plt.figure(figsize=(6.4,4.8*.85)) # new default is 6.4, 4.8
#plt.plot(mag, freq, 'ko-')
plt.axvspan(200, 1500, fc='None', hatch='//', label='Literature-based cohesion')
"""
for vwc in np.arange(0.01, 0.26, 0.01):
    rho = 2650 * .65 + vwc * 1000
    g = 9.8
    CohesionBlockFall_simple = rho * g * data['width']
    n, bins = np.histogram(CohesionBlockFall_simple, bins=np.arange(0, 3201, 100))
    mag = (bins[:-1] + bins[1:])/2.
    freq = n/float(np.sum(n))/np.diff(bins)
    plt.bar(bins[:-1], freq, width=np.diff(bins), ec='.8', fc='.8', align="edge", alpha=.2)
"""
plt.bar(bins[:-1], freq, width=np.diff(bins), ec="k", fc='.8', align="edge")
ax = plt.gca()
plt.ylabel('Probability density [Pa$^{-1}$, m$^{-1}$]', fontsize=16)
plt.xlabel('Shear stress at failure [Pa]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(0, ax.get_xlim()[-1])
# Equivalent landslide width
plt.twiny()
plt.xlim(np.array(ax.get_xlim()) / rho / g * 1000)
plt.xlabel('Landslide width [mm]', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
# Equivalent number of landslides
plt.twinx()
plt.ylim(np.array(ax.get_ylim()) * float(np.sum(n)) * np.mean(np.diff(bins)))
plt.ylabel('Number of landslides', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=12)
# Done!
plt.tight_layout()

if _save:
    plt.savefig('StressWidth.svg')

# LOOK UP SHEFFIELD THESIS 2015
# Sticking with Lu2009 and using internal friction angle + tensile strength to obtain cohesion


#plt.figure()
#plt.

if not _save:
    plt.show()
else:
    plt.close()
    
