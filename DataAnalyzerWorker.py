##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: February 23, 2022

@author: sarahzeichner
"""

import DataAnalyzerWithPeakInteg

#(1) Data Folder: Specify what folder contains the files to be analyzed
#inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Murchison AA/April 2021 - Dirty Murchison/April2021_DirtyMurchData/04292021"
#inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Murchison AA/Jan2022 - PristineMurchData/Raw Data/02032022/alanine"
inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Murchison AA/Jan2022 - PristineMurchData/Raw Data/reprocessed-dirtymurch/156"
#inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Murchison AA/Jan2022 - PristineMurchData/Raw Data/01312022"
#####################################################################

#(2) Experimental parameters
#(2a) Specify the order of isotopes as extracted from FT statistic
isotopeList = ['UnSub','13C']
#isotopeList = ['UnSub', '13C', '2x13C']
#isotopeList = ['UnSub', '2x13C']

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
#peakTimeFrames = [(18.6, 38)]
#peakTimeFrames = [(36.7,41.3)] #pristine murch, aspartic acid
#peakTimeFrames = [(19.45,110)] #pristine murch, beta alanine
#peakTimeFrames = [(13.8,110)] #pristine murch, alanine
peakTimeFrames = [(18.7,23.8)] #dirty murch, aspartic acid
#peakTimeFrames = [(8.15,30)] #dirty murch, beta alanine

#(2c) Specify time to extract out the background, similar to what is done in x-calibur when you highlight a part of the background
#backgroundNLTimes = [(8.6,8.8),(8.6,8.8),(8.6,8.8)]
#backgroundNLTimes = []
#backgroundNLTimes = [(36.45,36.5)] #pristine murch, aspartic
backgroundNLTimes = [(18.15,18.60)] #dirty murch, aspartic
#backgroundNLTimes = [(8.10,8.15)] #dirty murch, beta alanine
#backgroundNLTimes = [(13.75, 13.80)] #pristine murch, alanine
#backgroundNLTimes = [(19.35, 19.40)] #pristine murch, beta alanine

#(2d)Specify any ratios to omit
#omitRatios = ['15N/13C', '13C/15N'] #for nitrogen data
#omitRatios = ['13C/D',  'D/13C']
omitRatios = []

#####################################################################
#(3) Set toggles for data analysis

dataCullThreshhold = None #default = None; a target variable, like 'tic', or 'TIC*IT' to use to determine which scans to call. 
cullZeroScansOn = False #True/False
baselineSubstractionOn = True #True/False, when False, no baseline correction is applied to NL scores
gc_elution_on = True #True/False, when False, all scans are taken to calculate ratios
weight_by_NL_on = True
cullAmount = 2  
cullingThreshholdPercentMaxUnsubNLPeak = 0 #Threshhold to cull only above a threshhold for the values where NL/maxUnSub(NL)<threshhold

#####################################################################
 
#(4) Perform the data analysis and return the output! (Output is also exported to CSV in the underlying code, to the inputFolder path)
Output, StatsOutput = DataAnalyzerWithPeakInteg.calc_Folder_Output(inputStandardFolder, cullOn=dataCullThreshhold, cullAmount=cullAmount,\
                                                                    cullZeroScansOn=cullZeroScansOn,weightByNLOn = weight_by_NL_on,baselineSubstractionOn=baselineSubstractionOn, \
                                                                    gcElutionOn=gc_elution_on, gcElutionTimes = peakTimeFrames, backgroundNLTimes = backgroundNLTimes, \
                                                                    isotopeList = isotopeList, \
                                                                    minNL_over_maxNL=cullingThreshholdPercentMaxUnsubNLPeak, \
                                                                    omitRatios = omitRatios)