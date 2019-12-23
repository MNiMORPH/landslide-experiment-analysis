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
import ellipses as el
import pandas as pd
plt.ion()

def recursive_glob(rootdir='.', pattern='*'):
	"""Search recursively for files matching a specified pattern.
	
	Adapted from http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
	"""
	matches = []
	for root, dirnames, filenames in os.walk(rootdir):
	  for filename in fnmatch.filter(filenames, pattern):
		  matches.append(os.path.join(root, filename))
	return matches

def PolyArea(x,y):
    """
    Shoelace formula
    """
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def EllipseFit_LS(xy):
    lsqe = el.LSqEllipse()
    lsqe.fit([xy[:,0], xy[:,1]])
    if np.iscomplex(lsqe.width):
        print "WARN: COMPLEX"
        lsqe.width = np.nan
    if np.iscomplex(lsqe.height):
        print "WARN: COMPLEX"
        lsqe.height = np.nan
    if np.iscomplex(lsqe.phi):
        print "WARN: COMPLEX"
        lsqe.phi = np.nan
    if np.iscomplex(lsqe.center).any():
        print "WARN: COMPLEX"
        lsqe.center = np.array([np.nan, np.nan])
    # "width" and "height" are really half-width and half-height.
    smaja = np.max((lsqe.width, lsqe.height))
    smina = np.min((lsqe.width, lsqe.height))
    angle = lsqe.phi # CCW angle b/t x-axis and "width"
    if smaja == lsqe.height:
        if angle <= 0:
            angle += np.pi/2.
        else:
            angle -= np.pi/2.
    angle *= 180/np.pi
    # Assuming N-S is direction straight to river
    ls_width = 2 * ( lsqe.width**2 * np.cos(np.pi/2. - lsqe.phi)**2 +
                     lsqe.height**2 * np.sin(np.pi/2. - lsqe.phi)**2 )**.5
    ls_length = 2 * ( lsqe.width**2 * np.sin(np.pi/2. - lsqe.phi)**2 +
                      lsqe.height**2 * np.cos(np.pi/2. - lsqe.phi)**2 )**.5
        
    """
    if np.abs(angle) < np.pi/4.:
        width = 2 * lsqe.height # semi-minor axis
    else:
        width = 2 * lsqe.width # semi-major axis
    """
    return lsqe.center[0], lsqe.center[1], smaja, smina, angle, ls_width, \
           ls_length


####################
# INPUT FILE LISTS #
####################

shapefiles = sorted(recursive_glob(pattern='*.shp'))
DEMs = sorted(recursive_glob(pattern='DEM_fullextent_*.tif'))


#try:
#    data = pd.read_csv('LandslideSizes.txt')
#
#except:

##########################
# LANDSLIDE AREA & WIDTH #
##########################

areas = []
widths = []
lengths = []
semimajor_axes = []
semiminor_axes = []
angles_semimajor_axis_to_horizontal = []
x_center = []
y_center = []
lstime = []
for filename in shapefiles:
    sf = shapefile.Reader(filename)
    print filename
    shapes = sf.shapes()
    for shape in shapes:
        points = np.array(shape.points)
        _x_center, _y_center, _smaja, _smina, _angle, _ls_width, \
         _ls_length = EllipseFit_LS(points)
        x_center.append(_x_center)
        y_center.append(_y_center)
        widths.append( _ls_width ) # [m]
        lengths.append( _ls_length ) # [m]
        semimajor_axes.append( _smaja )
        semiminor_axes.append( _smina )
        angles_semimajor_axis_to_horizontal.append( _angle )
        areas.append( PolyArea(points[:,0], points[:,1]) ) # [m**2]
        # Timing of landslides -- experiment-agnostic
        timestr = filename.split('/')[-1].split('.')[0]
        lstime.append(int(timestr))
areas = np.array(areas)
widths = np.array(widths)
lengths = np.array(lengths)
lstime = np.array(lstime)
waiting = np.diff(lstime)
waiting = np.hstack(([np.nan], waiting)) # ignoring first, no earlier landslide

###############################
# DEM-BASED LANDSLIDE VOLUMES #
###############################

#Scan times
dem_seconds = []
for dem in DEMs:
    seconds = dem.split('.')[-2].split('_')[-1]
    dem_seconds.append(int(seconds))
dem_seconds = np.array(dem_seconds)

#shapefile times    
landslide_seconds = []
for sfile in shapefiles:
    shp_seconds = sfile.split('/')[-1].split('.')[0]
    landslide_seconds.append(int(shp_seconds))
landslide_seconds = np.array(landslide_seconds)

#choose the correct DEM 
previous_dem = []
for i in range(len(landslide_seconds)):
    dems_before_landslide = dem_seconds < landslide_seconds[i]
    previous_dem_time = np.max(dem_seconds[dems_before_landslide])
    previous_dem.append('DEMs/DEM_fullextent_'+str(previous_dem_time).zfill(7)+'.tif')

# GeoTIFF
# how do I use the previous_dem with the shapefiles to get volume
vol = []
time = []
for i in range(len(previous_dem)):
    previous = previous_dem[i]
    filename = shapefiles[i]   
    print previous
    ds = gdal.Open(previous)
    DEM = ds.ReadAsArray() #..... (from previous time)
    outarray = np.zeros(DEM.shape)  
    nY, nX = np.array(DEM.shape)
    Y = np.arange(0, nY, 1)[::-1]/1000.
    X = np.arange(0, nX, 1)/1000.
    sf = shapefile.Reader(filename)
    print filename
    shapes = sf.shapes()
    filename_time = filename.split('/')[-1].split('.')[0]
    time.append(filename_time)
    for shape in shapes:
        bbox = np.ceil(np.array(shape.bbox)*1E3)/1E3
        poly = Polygon(shape.points)
        x = np.round(np.arange(bbox[0], bbox[2], 0.001), 3)
        y = np.round(np.arange(bbox[1], bbox[3], 0.001), 3)
        #X, Y = np.meshgrid(x, y)
        #for xi in range(len(x)):
        #    for yi in range(len(y)):
        for xi in x:
            for yi in y:
                if poly.contains(Point(xi, yi)):
                    # Export the height above the minimum cell
                    # in the y-direiction at that x-location
                    # Should be lowest cell in valley
                    outarray[Y == yi, X == xi] = DEM[Y == yi, X == xi] - \
                                                 np.nanmin(DEM[:, X == xi])
        volume = np.nansum(outarray)/1E6 # mm cells to m, check DEM height.
        vol.append(volume)
volumes = np.array(vol)
depths = volumes/areas



##########
# OUTPUT #
##########

outdata = pd.DataFrame()
outdata['x_center [m]'] = x_center
outdata['y_center [m]'] = y_center
outdata['width [m]'] = widths
outdata['length [m]'] = lengths
outdata['depth [m]'] = depths
outdata['area [m2]'] = areas
outdata['semi-major axis length [m]'] = semimajor_axes
outdata['semi-minor axis length [m]'] = semiminor_axes
outdata['angle: semi-major axis to horizontal [deg]'] = \
    angles_semimajor_axis_to_horizontal # not necessarily CCW...
outdata['volume [m3]'] = volumes
outdata['runtime [s]'] = lstime
outdata['wait time [s]'] = waiting

outdata.to_csv('Landslides.txt', index_label='index')




"""
#####################
# COHESIVE STRENGTH #
#####################

rho = 2650 * .65 + .01 * 1000
g = 9.8

def shearStrength_blockFall(rho=rho, g=g, volumes=volumes, lengths=lengths, 
                            widths=widths, depths=depths):
    F_g = rho*g*volumes
    A_res = lengths*depths + widths*depths
    sigma_cohesion = F_g - A_res
    return sigma_cohesion

def shearStrength_blockFall_simple(rho=rho, g=g, volumes=volumes, 
                                   lengths=lengths, widths=widths,
                                   depths=depths):
    sigma_g = rho*g*widths
    sigma_cohesion = sigma_g
    return sigma_cohesion

def shearStrength_Topple_simple(rho=rho, g=g, widths=widths, depths=depths):
    return 2 * rho*g*widths**2 / (depths * (widths**2 + depths**2)**.5)

CohesionBlockFall = shearStrength_blockFall()
CohesionBlockFall_simple = shearStrength_blockFall_simple()
CohesionTopple_simple = shearStrength_Topple_simple()
"""









"""
def plf(x, a, b):
    return a*x**b

popt, pcov = curve_fit(plf, areas, vol)
print "a =", popt[0]
print "b =", popt[1]
_x = np.linspace(plt.xlim()[0], plt.xlim()[1], 10)
plt.plot(_x, plf(_x, popt[0], popt[1]), '-', color='0.5', linewidth=3)
plt.tight_layout()

#outdata = np.vstack((areas, vol)).transpose()
#np.savetxt('AreaVolume/AreaVolume.txt', outdata)
"""









"""
plt.figure()   
plt.loglog(areas, vol, 'ko')
plt.title('400mm/hr Experiment', fontsize=36, fontweight='bold')
plt.ylabel('Volume (m$^3$)', fontsize=30, fontweight='bold')
plt.xlabel('Area (m$^2$)', fontsize=30, fontweight='bold')
plt.tick_params(axis='both', which='major', labelsize=22)
"""

"""
#Scan times
dem_seconds = []
for dem in dems:
    seconds = dem.split('fullextent_')[-1].split('.')[0]
    dem_seconds.append(int(seconds))
dem_seconds = np.array(dem_seconds)

#shapefile times    
shapefile_seconds = []
for shapefile in shapefiles:
    shp_seconds = shapefile.split('.')[0]
    shapefile_seconds.append(int(shp_seconds))
shapefile_seconds = np.array(shapefile_seconds)

#choose the correct DEM 
previous_dem = []
for i in range(dem_seconds):
    dems_before_landslide = dem_seconds < shapefile_seconds[i]
    previous_dem_time = np.max(dem_seconds[dems_before_landslide])
    previous_dem.append(previous_dem_time)

#find volume from shapefile and preceding DEM
for volume in shapefile_seconds:
    point = Point(range(0, 3936), range(0,2401))
    polygon = Polygon(shapefile_seconds)

#Scan times
dem_seconds = []
for dem in dems:
    seconds = dem.split('.')[0].split('_')[-1]
    dem_seconds.append(int(seconds))
    #print seconds
    #print dem_seconds
dem_seconds = np.array(dem_seconds)

#shapefile times    
shapefile_seconds = []
for shapefile in shapefiles:
    shp_seconds = shapefile.split('/')[-1].split('.')[0]
    shapefile_seconds.append(int(shp_seconds))
    #print shp_seconds
    #print shapefile_seconds
shapefile_seconds = np.array(shapefile_seconds)

#choose the correct DEM 
previous_dem = []
for i in range(len(shapefile_seconds)):
    dems_before_landslide = dem_seconds < shapefile_seconds[i]
    previous_dem_time = np.max(dem_seconds[dems_before_landslide])
    previous_dem.append(previous_dem_time)
"""
