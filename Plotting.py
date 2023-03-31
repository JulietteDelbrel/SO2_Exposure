# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 10:21:24 2023

@author: q08996jd
"""

##############################################################################
# This is just to draw the new figures and not calculate the exposure
##############################################################################

##############################################################################
# Import libraries
##############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import os
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs

##############################################################################
# Flight_name
##############################################################################

flight_name0 = 'NT155_21Sept'
flight_name1 = 'NT156_21Sept'
flight_name2 = 'NT222_21Sept'
flight_name3 = 'NT223_21Sept'
flight_name4 = 'NT224_21Sept'
flight_name5 = 'NT227_21Sept'
flight_name6 = 'NT317_21Sept'
flight_name7 = 'NT318_21Sept'
flight_name8 = 'NT516_21Sept'
flight_name9 = 'NT517_21Sept'
flight_name10 = 'NT520_21Sept'
flight_name11 = 'NT591_21Sept'
# flight_name12 = 

##############################################################################
# Open csv files
##############################################################################

# open pixel data csv
prm_dir1 = ('D:/La_Palma_Flights/Sept_21/')
orbit_number = '20418'
pixel_results = '2021_09_21_' + orbit_number + '_pixel_results.csv'

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name0 + '.csv')
prm_dir = pd.read_csv(prm_dir1 + pixel_results)

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name0[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
              urcrnrlon = extent[1], urcrnrlat = extent[3],
              projection = 'cyl', resolution = 'h', fix_aspect = False, 
              ax = axs[1, 0], suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name0 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')

##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name1 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name1[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name1 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')


##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name2 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name2[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name2 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')


##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name3 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name3[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name3 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')


##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name4 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name4[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name4 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')



##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name5 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name5[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name5 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')



##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name6 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name6[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name6 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')

##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name7 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name7[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name7 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')

##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name8 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name8[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name8 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')

##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name9 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name9[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name9 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')


##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name10 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name10[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name10 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')

##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name11 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name11[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name11 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')

##############################################################################
##############################################################################
##############################################################################
# hello
##############################################################################
##############################################################################
##############################################################################

# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name12 + '.csv')

##############################################################################
# Rearrange data in the flight csv
##############################################################################

# extract time and date all in one column
nonfloatt = traj['UTC']
# split time into year, month, etc.
Year = nonfloatt.str[:4].astype(int)
Month = nonfloatt.str[5:7].astype(int)
Day = nonfloatt.str[8:10].astype(int)
Hour = nonfloatt.str[11:13].astype(int)
Minute = nonfloatt.str[14:16].astype(int)
Second = nonfloatt.str[17:19].astype(int)
# extract latitude (y), longitude (x), altitude (z)
nonfloaty = traj['Position'].str.split(',', expand=True)[0]
nonfloatx = traj['Position'].str.split(',', expand=True)[1]
nonfloatz = traj.Altitude*0.0003048*1000 # feet to km to m

x = [float(i) for i in nonfloatx]
y = [float(i) for i in nonfloaty]
z = [float(i) for i in nonfloatz]

# average out the time to closest 1 seconds
def round_time(dt=None, round_to=1):
    if dt == None:
        dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    rounding = (seconds+round_to/2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

# create dataframe with long, lat and altitude
df = pd.DataFrame({'Longitude': x, 'Latitude': y, 'Altitude' : z})
time = []
# extract the latitude, longitude, altitude every x time_step (in seconds). So doesn't plot the points from csv
time_step = 25 # change this number (the lower the closer the points, the longer the code runs - but the more accurate the flight traj)
for t in range(len(traj['UTC'])):
    round1 = (round_time(datetime.datetime(Year[t],Month[t],Day[t],Hour[t],Minute[t],Second[t]),round_to=time_step))
    injecttime = pd.to_datetime(round1,yearfirst = True)
    time.append(injecttime)
df.insert(0, 'RoundTime', time)
s = df.resample('25s', on='RoundTime').mean().interpolate() # change the '800s' to '10s' for eg. !!!Must be the same as time_step!!!
print(s) # prints the interpolated lng, lat, alt and time

lng = s.Longitude
lat = s.Latitude
alt = s.Altitude

##############################################################################
# Get pixels of interest only and arrange into polygons
##############################################################################

# only extract pixels that are within the range of the flight trajectory.
# Most plumes are absolutely huge so get rid of pixels that don't come close to flight
prm_dir = prm_dir[prm_dir.LON4 > min(lng)-0.1]
prm_dir = prm_dir[prm_dir.LON2 < max(lng)+0.1]
prm_dir = prm_dir[prm_dir.LAT3 < max(lat)+0.05]
prm_dir = prm_dir[prm_dir.LAT1 > min(lat)-0.05]
print(prm_dir)
prm_dir = prm_dir.reset_index(drop=True)

# extract the middle point of pixel
LON = prm_dir['LON']
LAT = prm_dir['LAT']

# Open necessary columns
Altitude = prm_dir['Injection_Altitude'] * 1000 #in m. altitude of pixel
ALT1 = prm_dir['Injection_Altitude_Low'] * 1000 #in m. lower error altitude of pixel
ALT2 = prm_dir['Injection_Altitude_High'] * 1000 #in m. higher error altitude of pixel
VCD = prm_dir['SO2_Interp_VCD']#in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
# Long and Lat of corners of pixels
LON1 = prm_dir['LON1']
LON2 = prm_dir['LON2']
LON3 = prm_dir['LON3']
LON4 = prm_dir['LON4']

LAT1 = prm_dir['LAT1']
LAT2 = prm_dir['LAT2']
LAT3 = prm_dir['LAT3']
LAT4 = prm_dir['LAT4']

# Create points
point1 = [LON1, LAT1, Altitude, VCD]
point2 = [LON2, LAT2, Altitude, VCD]
point3 = [LON3, LAT3, Altitude, VCD]
point4 = [LON4, LAT4, Altitude, VCD]

# Create polygons
x_values = [point1[0], point2[0], point3[0], point4[0], point1[0]] # longitude
y_values = [point1[1], point2[1], point3[1], point4[1], point1[1]] # latitude
z_values = [point1[2], point2[2], point3[2], point4[2], point1[2]] # altitude
c_values = [point1[3], point2[3], point3[3], point4[3], point1[3]] # amount of SO2

##############################################################################
# Fill in the pixels with VCD data
##############################################################################

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row")
fig.suptitle('Flight' +' ' + flight_name12[0:5] + ' ' + 'Orbit' + ' ' + orbit_number)
parameters = {'axes.labelsize': 25}
axs[0, 1].set_visible(False)

verts = zip(zip(LON1, LAT1), zip(LON2, LAT2), zip(LON3, LAT3), zip(LON4, LAT4))

cmap = matplotlib.cm.magma_r
bounds = np.linspace(0, 10, 11)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

##############################################################################
# Plot - Pixels
##############################################################################

# Make the collection and add it to the plot
coll = PolyCollection(verts, array = VCD, cmap = mpl.cm.magma_r, norm = norm, edgecolors = 'none')
axs[1, 0].add_collection(coll)
axs[1, 0].autoscale_view()

##############################################################################
# Plot - Basemap
##############################################################################

# only draw basemap where the pixels are
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
extent = [min(LON1) - 0.1, max(LON2) + 0.1, min(LAT1) - 0.1, max(LAT2) + 0.1]

# Create a basemap instance that draws the Earth layer
bm = Basemap(llcrnrlon = extent[0], llcrnrlat = extent[2],
             urcrnrlon = extent[1], urcrnrlat = extent[3],
             projection = 'cyl', resolution = 'h', fix_aspect = False, ax = axs[1, 0],
             suppress_ticks = False)

axs[1, 0].add_collection(bm.drawcoastlines(linewidth = 0.25))
axs[1, 0].add_collection(bm.drawcountries(linewidth = 0.35))

##############################################################################
# Plot - Aircraft (longitude vs. latitude)
##############################################################################

axs[1, 0].plot(x_values, y_values, 'k', alpha = 0.1)
axs[1, 0].scatter(lng, lat, c = alt, cmap = 'Blues', marker = '.')
axs[1, 0].set(xlabel = 'Longitude')
axs[1, 0].set(ylabel = 'Latitude')

##############################################################################
# Plot - Aircrafts (longitude vs. altitude)
##############################################################################

axs[0, 0].scatter(lng, alt, c = alt, cmap = 'Blues', marker = '.')
axs[0, 0].scatter(LON, Altitude, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[0, 0].set(ylabel = 'Altitude (m)')

##############################################################################
# Plot - Aircrafts (altitude vs. latitude)
##############################################################################

s = axs[1, 1].scatter(alt, lat, c = alt, cmap = 'Blues', marker = '.')
t = axs[1, 1].scatter(Altitude, LAT, c = VCD, cmap = 'magma_r', marker = '.', vmin = 0, vmax = 10)
axs[1, 1].set(xlabel = 'Altitude (m)')

##############################################################################
# Plot - colorbars
##############################################################################

fig.subplots_adjust(right = 1)
cbar_ax = fig.add_axes([0.65, 0.55, 0.03, 0.35])
cb = fig.colorbar(s, cax = cbar_ax)
cb.set_label('Altitude [m]')

# Add a colorbar for the PolyCollection
cbar_ax = fig.add_axes([0.85, 0.55, 0.03, 0.35])
cb = fig.colorbar(coll, cax = cbar_ax)
cb.set_label('VCD [DU]')

##############################################################################
# Save figure
##############################################################################

fig.savefig(prm_dir1 + flight_name12 + '_' + orbit_number + '_map' + '.jpeg', format='jpeg', dpi=1200, bbox_inches = 'tight')
