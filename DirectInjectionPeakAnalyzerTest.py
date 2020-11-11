##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: Tues Nov 10, 2020

@author: sarahzeichner
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.integrate import simps

#####################################################################
########################## CONSTANTS ################################
#####################################################################

WINDOW_LENGTH  = 5
SLOPE_THRESHHOLD = 0.008
NAN_REPLACER = 0.0000001

#####################################################################
########################## FUNCTIONS ################################
#####################################################################

def calc_slope(x):
    slope = np.polyfit(range(len(x)), x, 1)[0]
    return slope

def gaussian(x, a, b, c):
    return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))

#####################################################################
########################## MAIN ####################################
#####################################################################

# Generate dummy dataset
x_dummy = np.linspace(start=-15, stop=15, num=100)
y_dummy = norm.pdf(x_dummy, 0, 1)
# Add noise from a Gaussian distribution
noise = 0.01 * np.random.normal(size=y_dummy.size)
y_dummy = y_dummy + noise

#set up the dataframe with x, y, and rolling slopes
data = pd.DataFrame(y_dummy, columns=['y'])
data['slope'] = data.rolling(WINDOW_LENGTH).apply(calc_slope)
data['slope'] = data['slope'].fillna(NAN_REPLACER)
data['abs_slope'] =  abs(data['slope'])
data['x'] = x_dummy

#Find threshhold and choose subset of data based on slope threshhold
data = data[data['abs_slope'] > SLOPE_THRESHHOLD] 

#Fit a gaussian to the data
popt, pcov = curve_fit(gaussian, data['x'], data['y'])
plt.scatter(x_dummy, y_dummy)
plt.plot(data['x'], gaussian(data['x'], *popt), 'g--', 
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

#integrate the polynomial to get area under the curve
area = simps(gaussian(data['x'], *popt))

print(area)
