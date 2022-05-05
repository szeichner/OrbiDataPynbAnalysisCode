##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: February 23, 2022

@author: sarahzeichner
"""

import DataAnalyzerWithPeakInteg

#(1) Data Folder: Specify what folder contains the files to be analyzed
#inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Murchison AA/April 2021 - Dirty Murchison/April2021_DirtyMurchData/04292021"
#inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Murchison AA/Jan2022 - PristineMurchData/Raw Data/reprocessed-dirtymurch/156"
inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Hayabusa2/raw data/04182022 and 04192022 - 2Da window/m plus 1 data/fluor_pyrene"
#inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Murchison AA/Jan2022 - PristineMurchData/Raw Data/01312022"
#####################################################################

#(2) Experimental parameters
#(2a) Specify the order of isotopes as extracted from FT statistic
#isotopeList = ['UnSub','13C']
isotopeList = ['UnSub', 'D', '2x13C']
#isotopeList = ['UnSub', 'D','2x13C']

#(2b) Specify the time frame of each eluted peak
#peakTimeFrames = [(6.65,6.78),(7.85, 8.12), (10,10.12),(12.37,12.56)] #sample mixture
#peakTimeFrames = [(6.65,6.78),(7.85, 8.12), (10,10.12), (10.21, 10.4),(12.37,12.56)] #standard mixture
#peakTimeFrames = [(11.37, 11.60), (11.37, 11.60), (11.37, 11.60), (11.37, 11.60)] #aspartic labels
#peakTimeFrames = [(8.8,11.0), (8.8,11.0), (8.8,11.0)] #serine fragments, from elise SERC
#peakTimeFrames = [(6.65,6.78)] #alanine, std mixture
#peakTimeFrames = [(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12)] #aspartic
#peakTimeFrames = [(12.37,12.56),(12.37,12.56),(12.37,12.56),(12.37,12.56),(12.37,12.56), (12.37,12.56), (12.37,12.56), (12.37,12.56), (12.37,12.56)] #methionine
#peakTimeFrames = [(5.50, 7.70),(5.50, 7.70),(5.50, 7.70),(5.50, 7.70)] #Chimiak dallas alanine, all peaks--> note the long time frame
#peakTimeFrames = [(32.80, 33.50)] #fluoranthene large window
#peakTimeFrames = [(34.36, 35.12)] #pyrene large windowxylene
#peakTimeFrames = [(32.65, 33.35), (32.65, 33.35), (32.65, 33.35)] #fluoranthene large window
#peakTimeFrames = [(34.15, 34.80), (34.15, 34.80), (34.15, 34.80)] #pyrene large window
#peakTimeFrames = [(8.54, 8.86), (8.54, 8.86)] #p-xylene
#peakTimeFrames = [(68.75, 70.0)]
#peakTimeFrames = [(8.54, 8.86)] #p-xylene
#peakTimeFrames = [(4.56, 6.7)] #chimiak dallas alanine, 184 and 140 peaks

#Murchison project time frames
#peakTimeFrames = [(18.65, 21.7)] #pristine murch, single column aspartic 156
#peakTimeFrames = [(36.75,41.3)] #pristine murch, aspartic acid
#peakTimeFrames = [(19.7,39)] #pristine murch, beta alanine
#peakTimeFrames = [(13.9,26.3)] #pristine murch, alanine
#peakTimeFrames = [(18.7,26.5)] #dirty murch, aspartic acid 156
#peakTimeFrames = [(18.9,24.8)] #dirty murch, aspartic acid 113
#peakTimeFrames = [(8.5,10.2)] #dirty murch, beta alanine

#HB2 project time frames - longer ramp
#peakTimeFrames = [(29.5,29.9)] #napthalene - longer ramp
#peakTimeFrames = [(71.66,72.05)] #phenanthrene - longer ramp
#peakTimeFrames = [(72.05,73.3)] #anthracene - longer ramp
#peakTimeFrames = [(92.15, 93.66)] #fluoranthene - longer ramp
#peakTimeFrames = [(95.35, 97.96)] #pyrene - longer ramp

#peakTimeFrames = [(20.7,22.15)] #napthalene - shorter ramp, a106
#peakTimeFrames = [(20.5,22.2)] #napthalene - shorter ramp, c107
#peakTimeFrames = [(20.5,20.8)] #napthalene - shorter ramp, 2Da window
#peakTimeFrames = [(61.1,62.1)] #phenanthrene - shorter ramp, a106
#peakTimeFrames = [(61.2,62.4)] #phenanthrene - shorter ramp, c107
#peakTimeFrames = [(61.13,61.8)] #phenanthrene - shorter ramp, 2Da window
#peakTimeFrames = [(62.1,63.35)] #anthracene - shorter ramp, a106
#peakTimeFrames = [(62.4,64)] #anthracene - shorter ramp, c107
#peakTimeFrames = [(62.2,62.85)] #anthracene - shorter ramp, 2Da window
#peakTimeFrames = [(81.95,83.32)] #fluoranthene - shorter ramp, a106
#peakTimeFrames = [(81.95,83.32)] #fluoranthene - shorter ramp, c107
#peakTimeFrames = [(82,83.05)] #fluoranthene - shorter ramp, 2Da window
#peakTimeFrames = [(85.3,87)] #pyrene - shorter ramp, a106
#peakTimeFrames = [(85.3,87)] #pyrene - shorter ramp, c107
peakTimeFrames = [(85.25,87)] #pyrene - shorter ramp, 2Da window

#peakTimeFrames = [(20.7,22.15),(61.1,62.1),(62.1,63.35),(81.95,83.32),(85.3,87)] #all peaks - shorter ramp - this does not process 5 together! beware

#(2c) Specify time to extract out the background, similar to what is done in x-calibur when you highlight a part of the background

#HB2 project
#backgroundNLTimes = [(28.8,29.2)] #napthalene
#backgroundNLTimes = [(70,70.5)] #phenanthrene
#backgroundNLTimes = [(72,72.5)] #anthracene
#backgroundNLTimes = [(91,91.5)] #fluoranthene
#backgroundNLTimes = [(95.35,97.96)] #pyrene
backgroundNLTimes = []

#backgroundNLTimes = [(8.6,8.8),(8.6,8.8),(8.6,8.8)]
#backgroundNLTimes = []
#backgroundNLTimes = [(36.45,36.5)] #pristine murch, aspartic
#backgroundNLTimes = [(17.9,18.2)] #pristine murch, aspartic
#backgroundNLTimes = [(18.15,18.45)] #dirty murch, aspartic
#backgroundNLTimes = [(8.2,8.4)] #dirty murch, beta alanine
#backgroundNLTimes = [(13.7, 13.9)] #pristine murch, alanine
#backgroundNLTimes = [(19.4, 19.55)] #pristine murch, beta alanine

#(2d)Specify any ratios to omit
#omitRatios = ['15N/13C', '13C/15N'] #for nitrogen data
#omitRatios = ['13C/D',  'D/13C']
omitRatios = ['2x13C/D','D/2x13C']
#omitRatios = ['2x13C/13C']

#####################################################################
#(3) Set toggles for data analysis

dataCullThreshhold = None #default = None; a target variable, like 'tic', or 'TIC*IT' to use to determine which scans to call. 
cullZeroScansOn = False #True/False
baselineSubstractionOn = False #True/False, when False, no baseline correction is applied to NL scores
gc_elution_on = True #True/False, when False, all scans are taken to calculate ratios
weight_by_NL_on = False
cullAmount = 2  
cullingThreshholdPercentMaxUnsubNLPeak = 0.05 #Threshhold to cull only above a threshhold for the values where NL/maxUnSub(NL)<threshhold
MAX_IT = 3000

#####################################################################
 
#(4) Perform the data analysis and return the output! (Output is also exported to CSV in the underlying code, to the inputFolder path)
Output, StatsOutput = DataAnalyzerWithPeakInteg.calc_Folder_Output(inputStandardFolder, cullOn=dataCullThreshhold, cullAmount=cullAmount,\
                                                                    cullZeroScansOn=cullZeroScansOn,weightByNLOn = weight_by_NL_on,baselineSubstractionOn=baselineSubstractionOn, \
                                                                    gcElutionOn=gc_elution_on, gcElutionTimes = peakTimeFrames, backgroundNLTimes = backgroundNLTimes, \
                                                                    isotopeList = isotopeList, \
                                                                    minNL_over_maxNL=cullingThreshholdPercentMaxUnsubNLPeak, \
                                                                    maxIT = MAX_IT, omitRatios = omitRatios)