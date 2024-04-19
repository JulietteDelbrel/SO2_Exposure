# -*- coding: utf-8 -*-
"""
Created on Tue May 17 11:40:46 2022

@author: juliette
"""

# import libraries

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haversine import haversine, Unit
from deffunction import doIntersect, line_intersection
from point_in_polygon import is_inside_polygon
import datetime
import math
from numpy import log as ln
import csv
from sympy import Point
from scipy.integrate import simpson

# add gauss1d for points of intersection
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

prm_dir1 = ('D:/mydirectory/') # open directory

##############################################################################
# OPEN CSV FILE FOR TRAJ AND INTERPOLATE
##############################################################################

flight_name = 'NT424_2960a0ac' # example

orbit_number = '20589' # example
pixel_results = '2021_10_03_' + orbit_number + '_pixel_results.csv' # example
# open flight csv
traj = pd.read_csv(prm_dir1 + flight_name + '.csv')
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
# Open csv file FOR PIXELS
##############################################################################

# open files with pixel data (plume data)
prm_dir = pd.read_csv(prm_dir1 + pixel_results)

##############################################################################
# Optimising the code
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
VCD = prm_dir['SO2_Interp_VCD'] #in DU (converted to mol/m2 lower down). how much SO2 is in the pixel
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
# Create csv file
##############################################################################

column_name = ["Seconds", "Cumulative_exposure"] #The name of the columns
data = [] # the data

# creates a csv file
with open(prm_dir1 + flight_name + '_' + orbit_number + '.csv','w', newline='') as f:
    writer = csv.writer(f) # this is the writer object
    writer.writerow(column_name) # this will list out the names of the columns which are always the first entrries
    writer.writerow(data)

##############################################################################
# DRIVER CODE
##############################################################################
# molar mass of SO2 multiplied by the number of molecules in 1 DU
molar_mass = 64.066*1000000 # µg/mol
DU = 0.00044615 # mol/m2

# multiple to use for full width half maximum
multiple = (2*math.sqrt(2*ln(2)))

# Driver program to test above functions:
tot_exp = 0
tot_time = 0

for f in range(len(lng)):
    m = f + 1
    if m > len(lng)-1:
        break

    print(f, '-', m)
    p2 = Point(lng[f], lat[f])
    p2_0 = (lng[f], lat[f])
    q2 = Point(lng[m], lat[m])
    q2_0 = (lng[m], lat[m])
   # 7. Append data to csv file
    dict = {"Seconds": f*time_step, "Cumulative_exposure":tot_exp}
    with open(prm_dir1 + flight_name + '_' + orbit_number + '.csv', 'a', newline='') as csv_file:
        dict_object = csv.DictWriter(csv_file, fieldnames=column_name)
        dict_object.writerow(dict)  

    for g in range(len(LON)):
        success = []

        n = 0
        p1 = Point(LON4[g], LAT4[g])
        p1_0 = (LON4[g], LAT4[g])
        q1 = Point(LON1[g], LAT1[g])
        q1_0 = (LON1[g], LAT1[g])

        if doIntersect(p1, q1, p2, q2):
            n = n + 1
            int1 = (line_intersection((p1_0, q1_0), (p2_0, q2_0)))
            success.append(int1)

        p1 = Point(LON1[g], LAT1[g])
        p1_1 = (LON1[g], LAT1[g])
        q1 = Point(LON2[g], LAT2[g])
        q1_1 = (LON2[g], LAT2[g])

        if doIntersect(p1, q1, p2, q2):
            n = n + 1
            int2 = (line_intersection((p1_1, q1_1), (p2_0, q2_0)))
            success.append(int2)

        p1 = Point(LON2[g], LAT2[g])
        p1_2 = (LON2[g], LAT2[g])
        q1 = Point(LON3[g], LAT3[g])
        q1_2 = (LON3[g], LAT3[g])
       
        if doIntersect(p1, q1, p2, q2):
            n = n + 1
            int3 = (line_intersection((p1_2, q1_2), (p2_0, q2_0)))
            success.append(int3)

        p1 = Point(LON3[g], LAT3[g])
        p1_3 = (LON3[g], LAT3[g])
        q1 = Point(LON4[g], LAT4[g])
        q1_3 = (LON4[g], LAT4[g])
       
        if doIntersect(p1, q1, p2, q2):
            n = n + 1
            int4 = (line_intersection((p1_3, q1_3), (p2_0, q2_0)))
            success.append(int4)
       
        sigma_low = (abs(Altitude[g]-ALT1[g])*2)/(multiple)
        x_values1 = np.linspace(0, Altitude[g], 100)
        sigma_hi = (abs(ALT2[g]-Altitude[g])*2)/(multiple)
        x_values2 = np.linspace(Altitude[g], 10000, 100)
        gaussian1 = gaussian(x_values1, Altitude[g], sigma_low)
        gaussian2 = gaussian(x_values2, Altitude[g], sigma_hi)
        area = simpson(gaussian1, x_values1) + simpson(gaussian2, x_values2)
        scale = VCD[g]/area
        
        if n == 2:
        # the 25 second segment crosses a pixel twice (enters and exits)
            # this will measure the concentration in µg/m3
            Concentration = []
            if alt[f] < Altitude[g]:
                y = np.interp(alt[f], x_values1, gaussian(x_values1, Altitude[g], sigma_low)*scale) # DU/m 
                # line above: without the scale it's DU, so you adjust it and spreads the DU per every metre
                Concentration = (y * molar_mass * DU)
                print('n = 2; Concentration:', Concentration, ' µg/m3')
                # print('VCD and y:', VCD[g], y)
            elif alt[f] > Altitude[g]:
                y = np.interp(alt[f], x_values2, gaussian(x_values2, Altitude[g], sigma_hi)*scale)
                Concentration = (y * molar_mass * DU)
                print('n = 2; Concentration:', Concentration, ' µg/m3')
                # print('VCD and y:', VCD[g], y)
            # to calculate the exposure, we need to do concentration/time of exposure
            int_value1 = success[0][1], success[0][0]
            int_value2 = success[1][1], success[1][0]
            distance1 = haversine(int_value1, int_value2, unit = Unit.KILOMETERS)
        # 1. Measure distance between coordinates
            first_point = lat[f], lng[f]
            second_point = lat[m], lng[m]
            distance = haversine(first_point, second_point, unit = Unit.KILOMETERS)
        # 2. Calculate time between the coordinates
                # time = ti[m] - ti[f]
                #alternative?
            time = time_step
        # 3. Calculate speed between the coordinates
            speed = distance/time
        # 5. Calculate time between the two intersecting points
            time1 = distance1/speed
        # 6. Calculate exposure by multiplying speed with concentration
            tot_exp = ((Concentration*time1)) + tot_exp
            print('Cumulative exposure :' , tot_exp)
              
        if n == 1:
        # if only 1 part of 25 second segment is in a pixel
            inside = []
            Concentration = []
            if alt[f] < Altitude[g]:
                y = np.interp(alt[f], x_values1, gaussian(x_values1, Altitude[g], sigma_low)*scale)
                Concentration = (y * molar_mass * DU)
                print('n = 1; Concentration:', Concentration, ' µg/m3')
                # print('VCD and y:', VCD, y)
            elif alt[f] > Altitude[g]:
                y = np.interp(alt[f], x_values2, gaussian(x_values2, Altitude[g], sigma_hi)*scale)
                Concentration = (y * molar_mass * DU)              
                print('n = 1; Concentration:', Concentration, ' µg/m3')
                # print('VCD and y:', VCD, y)
     
        # 1. Measure distance between coordinates
            first_point = lat[f], lng[f]
            second_point = lat[m], lng[m]
            distance = haversine(first_point, second_point)
        # 2. Calculate time between the coordinates
            # time = ti[m] - ti[f]
            #alternative?
            time = time_step
        # 3. Calculate speed between the coordinates
            speed = distance/time
        # 4. Determine if first_point is inside polygon
            polygon1 = [ (point1[0][g], point1[1][g]), (point2[0][g], point2[1][g]),
                            (point3[0][g], point3[1][g]), (point4[0][g], point4[1][g]) ]
            int_value1 = success[0][1], success[0][0]
            p10 = (lng[f], lat[f])
            if (is_inside_polygon(points = polygon1, p = p10)):
                # 6. Determine distance between PoI and first_point
                distance1 = haversine((lat[f], lng[f]), int_value1)
                # 7. Calculate time
                time1 = distance1/speed
                # Exposure is time * concentration
                tot_exp = (Concentration*time1) + tot_exp
                print('Cumulative exposure :' , tot_exp)
            # 5. Determine if second_point is inside polygon
            int_value1 = success[0][1], success[0][0]
            p20 = (lng[m], lat[m])
            if (is_inside_polygon(points = polygon1, p = p20)):
                # 6bis. Calculate distance between point of intersection and coordinate
                distance1 = haversine((lat[m], lng[m]), int_value1)   
                # 7. Calculate time
                time1 = distance1/speed
                # Exposure is time * concentration
                tot_exp = (Concentration*time1) + tot_exp
                print('Cumulative exposure :' , tot_exp)
        
        if n == 0:
        # Both points are either inside the polygon or both are outside it on one side - no crossing the pixel
            inside = []
            if __name__ == '__main__':
                polygon1 = [ (point1[0][g], point1[1][g]), (point2[0][g], point2[1][g]),
                            (point3[0][g], point3[1][g]), (point4[0][g], point4[1][g]) ]
                p10 = (lng[f], lat[f])
                if (is_inside_polygon(points = polygon1, p = p10)):
                    inside.append(p10)
                   
                p20 = (lng[m], lat[m])
                if (is_inside_polygon(points = polygon1, p = p20)):
                    inside.append(p20)
                   
            if len(inside) == 2:
                Concentration = []
                if alt[f] < Altitude[g]:
                    y = np.interp(alt[f], x_values1, gaussian(x_values1, Altitude[g], sigma_low)*scale)
                    Concentration = (y * molar_mass * DU)
                    # print('VCD and y:', VCD, y)
                    print('n = 0; Concentration:', Concentration, ' µg/m3')
                elif alt[f] > Altitude[g]:
                    y = np.interp(alt[f], x_values2, gaussian(x_values2, Altitude[g], sigma_hi)*scale)
                    Concentration = (y * molar_mass * DU)
                    print('n = 0; Concentration:', Concentration, ' µg/m3')
                    # print('VCD and y:', VCD, y)
                time1 = time_step
                tot_exp = (Concentration*time1) + tot_exp
                print('Cumulative exposure :', tot_exp)

##############################################################################
# Plot Cumulative Exposure graph
##############################################################################
table = pd.read_csv(prm_dir1 + flight_name + '_' + orbit_number + '.csv')
fig, ax = plt.subplots()
plt.plot(table.Seconds/60, table.Cumulative_exposure)
# plt.plot(table.Seconds, cumuexpo)

plt.xlabel('Time (minutes)', fontsize = 13)
plt.ylabel('Cumulative Exposure (µg/m$^3$)', fontsize = 13)
plt.tight_layout()
fig.savefig(prm_dir1 + flight_name + '_' + orbit_number + '.jpeg', format='jpeg', dpi=1200)

