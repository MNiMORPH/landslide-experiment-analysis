# ADW, 07 June 2016

import glob
import datetime as dt
import numpy as np
import os 
import shutil
from matplotlib import pyplot as plt
import scanDateTime

expName = '151117_T_DW_IW_01'

# Run this inside the directory with the images
cwd = 'E:\\Fluvial 2015\\' + expName + '\\' + expName + '_Photos'
os.chdir(cwd)

filenames = sorted(glob.glob('*.jpg'))
f0 = filenames[0]

timestamps = []

# For the "strptime" function for working with time stamps,
# see http://www.tutorialspoint.com/python/time_strptime.htm
experiment_start = dt.datetime.strptime(filenames[0], "Img%Y-%m-%d %H.%M.%S.jpg")
seconds_since_start = []

for fi in filenames:
  timestamp = dt.datetime.strptime(fi, "Img%Y-%m-%d %H.%M.%S.jpg")
  timestamps.append(timestamp)
  # Writing this as an integer instead of timedelta with total_seconds
  seconds_since_start.append( (timestamp - experiment_start).total_seconds() )
seconds_since_start = np.array(seconds_since_start)

timestampsara = np.array(timestamps)    #array of timestamps

# Get start and end times of this experiment run from scanDateTime.py
scan_start_times, scan_end_times = scanDateTime.get_start_end_times(expName)
        
time_elapsed_during_pause = []   #create list of time differentials
#add to list ("scan end times") - ("scan start times") in seconds
for i in range(len(scan_start_times)):
    time_elapsed_during_pause.append(\
        (scan_end_times[i] - scan_start_times[i]).total_seconds() )

#This gives us the differences between the end and the start times.
#scan_end_times - scan_start_times

plt.plot(seconds_since_start, 'k.');

seconds_since_start_modified = seconds_since_start.copy()
#What we want, though, is for the time between the end and the start, 
#instead of being this difference value, to be twenty seconds for all of them.
for i in range(1, len(scan_start_times)+1):
   start_time = scan_start_times[-i]
   after_start = timestampsara > start_time
   # This is going to change the array of the time between the scans
   seconds_since_start_modified[after_start] = \
               seconds_since_start_modified[after_start] - \
               time_elapsed_during_pause[-i] + 20

# Create a new list of filenames with time steps in seconds
filenames_new = []
for timestep in seconds_since_start_modified:
    filenames_new.append( 'ImgSec_' + '%07d' %timestep )

"""
# Copy the files to a processed images folder
try:
    os.mkdir('..\\Photos_timestamped')
except:
    pass

for i in range(len(filenames)):
    print filenames[i], '-->', filenames_new[i]
    shutil.copyfile(filenames[i], '..\\Photos_timestamped\\' + filenames_new[i])
"""

# Plot the results

plt.plot(seconds_since_start_modified, 'g.');
plt.xlabel('Image number (sequential from start)', fontsize=16)
plt.ylabel('Real time since start of experiment [seconds]', fontsize=16)

plt.show()


# ... after this, it will be a matter of looking for legitimate breaks in time
# that are related to pauses in the experiment (e.g., laser scans)
# and replacing them with a number of seconds elapsed that comes
# immediately after those time-stamps that precede the laser scans.

# Also, using this method of working through the time stamps allows you to
# delete pictures before you sort them -- because it just reads the 
# time stamps. Just so long as you don't create a big break in time without
# remembering that it was you who made it!

