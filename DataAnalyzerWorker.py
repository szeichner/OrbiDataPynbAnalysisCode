##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: Thursday, Nov 19 2020

@author: sarahzeichner
"""

import DataAnalyzerWithPeakInteg

#Specify what folder contains the files to be analyzed
inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Direct Injection/DirectInjectionData/022021_pxylene/m+1"

#Specify the order of isotopes as extracted from FT statistic
#isotopeList = ['UnSub','13C']
isotopeList = ['UnSub', 'D']
#isotopeList = ['UnSub', '2x13C']

#Specify whether to cull based on time frames for chromatography/direct elution
gc_elution_on = True

#Specify the time frame of each eluted peak
#peakTimeFrames = [(6.65,6.78),(7.85, 8.12), (10,10.12),(12.37,12.56)] #sample mixture
#peakTimeFrames = [(6.65,6.78),(7.85, 8.12), (10,10.12), (10.21, 10.4),(12.37,12.56)] #standard mixture
#peakTimeFrames = [(11.37, 11.60), (11.37, 11.60), (11.37, 11.60), (11.37, 11.60)] #aspartic labels
#peakTimeFrames = [(14.05,14.65), (14.05,14.65), (14.05,14.65)] #serine fragments, from elise SERC
#peakTimeFrames = [(6.65,6.78)] #alanine, std mixture
#peakTimeFrames = [(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12)] #aspartic
#peakTimeFrames = [(12.37,12.56),(12.37,12.56),(12.37,12.56),(12.37,12.56),(12.37,12.56), (12.37,12.56), (12.37,12.56), (12.37,12.56), (12.37,12.56)] #methionine
#peakTimeFrames = [(5.50, 7.70),(5.50, 7.70),(5.50, 7.70),(5.50, 7.70)] #Chimiak dallas alanine, all peaks--> note the long time frame

#PAH data
#peakTimeFrames = [(32.80, 33.50)] #fluoranthene large window
#peakTimeFrames = [(34.36, 35.12), (34.36, 35.12), (34.36, 35.12)] #pyrene large window\

#p-xylene
#peakTimeFrames = [(8.54, 8.86), (8.54, 8.86)]
peakTimeFrames = [(8.54, 8.86)]

#peakTimeFrames = [(32.65, 33.35), (32.65, 33.35), (32.65, 33.35)] #fluoranthene large window
#peakTimeFrames = [(34.15, 34.80), (34.15, 34.80), (34.15, 34.80)] #pyrene large window


#peakTimeFrames = [(4.56, 6.7)] #chimiak dallas alanine, 184 and 140 peaks

#Specify any ratios to omit
#omitRatios = ['15N/13C', '13C/15N'] #for nitrogen data
#omitRatios = ['13C/D',  'D/13C']
omitRatios = []
 
#Perform the data analysis and return the output! (Output is also exported to CSV in the underlying code, to the inputFolder path)
Output, StatsOutput = DataAnalyzerWithPeakInteg.calc_Folder_Output(inputStandardFolder, cullOn=None, cullZeroScansOn=False, gcElutionOn=gc_elution_on, gcElutionTimes = peakTimeFrames,  cullAmount=2, isotopeList = isotopeList, NL_over_TIC=0.10, omitRatios = omitRatios, fileCsvOutputPath=None)


