##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: Wed Sept 1, 2021
@author: sarahzeichner
"""

import DataAnalyzerWithPeakInteg

#(1) Data Folder: Specify what folder contains the files to be analyzed
inputStandardFolder = ""

#####################################################################

#(2) Experimental parameters
#(2a) Specify the order of isotopes as extracted from FT statistic
#Some examples:
#isotopeList = ['UnSub', '13C', '2x13C']
#isotopeList = ['UnSub', '2x13C']
isotopeList = ['UnSub','13C']


#(2b) Specify the time frame of each eluted peak
#e.g., peakTimeFrames = [(8.8,11.0), (8.8,11.0), (8.8,11.0)] #analysis of 3 fragments of one eluted peak
peakTimeFrames = [(8.8,11.0), (8.8,11.0), (8.8,11.0)]

#(2c) Specify time to extract out the background. Length of the array must be same length as the number of fragments
# you are analyzing. (e.g., backgroundNLTimes = [(8.6,8.8),(8.6,8.8),(8.6,8.8)])
backgroundNLTimes = []

#(2d)Specify any ratios to omit if you are analyzing multiple peaks
# e.g., omitRatios = ['13C/D',  'D/13C']
omitRatios = []

#####################################################################
#(3) Set toggles for data analysis

dataCullThreshhold = None #default = None; a target variable, like 'tic', or 'TIC*IT' to use to determine which scans to call. 
cullZeroScansOn = False #True/False
baselineSubstractionOn = False #True/False, when False, no baseline correction is applied to NL scores
gc_elution_on = True #True/False, when False, all scans are taken to calculate ratios
weight_by_NL_on = False
cullAmount = 2  
percentCullThreshhold = 0.10 #Threshhold to cull only above a threshhold for the values where NL/maxUnSub(NL)<threshhold

#####################################################################
 
#(4) Perform the data analysis and return the output! (Output is also exported to CSV in the underlying code, to the inputFolder path)
Output, StatsOutput = DataAnalyzerWithPeakInteg.calc_Folder_Output(inputStandardFolder, cullOn=dataCullThreshhold, cullAmount=cullAmount,\
                                                                    cullZeroScansOn=cullZeroScansOn,\
                                                                    baselineSubstractionOn=baselineSubstractionOn, \
                                                                    gcElutionOn=gc_elution_on, \
                                                                    gcElutionTimes = peakTimeFrames, \
                                                                    backgroundNLTimes = backgroundNLTimes, \
                                                                    isotopeList = isotopeList, \
                                                                    minNL_over_maxNL=percentCullThreshhold, \
                                                                    omitRatios = omitRatios)