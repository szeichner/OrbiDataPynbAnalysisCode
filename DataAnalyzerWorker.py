##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: Thursday, Nov 19 2020

@author: sarahzeichner
"""

import DataAnalyzerWithPeakInteg

#Specify what folder contains the files to be analyzed
inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Direct Injection/DirectInjectionData/direct injection reanalysis/major peaks/standard"

#Specify the order of isotopes as extracted from FT statistic
isotopeList = ['UnSub','15N','13C']
#isotopeList = ['UnSub','13C']

#Specify whether to cull based on time frames for chromatography/direct elution
gc_elution_on = True

#Specify the time frame of each eluted peak
#peakTimeFrames = [(7.85, 8.12), (10,10.12), (12.37,12.56)] #sample mixture
peakTimeFrames = [(7.85, 8.12), (10,10.12), (10.21, 10.4), (12.37,12.56)] #standard mixture
#peakTimeFrames = [(6.65,6.78)] #alanine
#peakTimeFrames = [(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12),(10,10.12)] #aspartic
#peakTimeFrames = [(12.37,12.56),(12.37,12.56),(12.37,12.56),(12.37,12.56),(12.37,12.56), (12.37,12.56)] #methionine


#Specify any ratios to omit
#omitRatios = ['15N/13C', '13C/15N'] #for nitrogen data
omitRatios = []
 
#Perform the data analysis and return the output! (Output is also exported to CSV in the underlying code, to the inputFolder path)
Output, StatsOutput = DataAnalyzerWithPeakInteg.calc_Folder_Output(inputStandardFolder, cullOn=None, cullZeroScansOn=False, gcElutionOn=gc_elution_on, gcElutionTimes = peakTimeFrames,  cullAmount=2, isotopeList = isotopeList, NL_over_TIC=0.10, omitRatios = omitRatios, fileCsvOutputPath=None)


