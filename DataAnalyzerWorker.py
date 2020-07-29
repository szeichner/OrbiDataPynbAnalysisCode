##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: Thurs July 16, 2020

@author: sarahzeichner
"""

import DataAnalyzer

#Change these things to test the different code, or comment out if you're using in conjunction with the python notebook
inputStandardFolder = "/Users/sarahzeichner/Documents/Caltech/Research/Direct Injection/DirectInjectionData/July2020/July172020/Standard"
#inputStandardFile = "/Users/sarahzeichner/Documents/Caltech/Research/Quick Orbitrap Methods/data/June2020/AA_std_2_15_agc_2e4.xlsx"
#isotopeList = ['UnSub','15N','13C']
isotopeList = ['UnSub','15N','13C']
gc_elution_on = True
peakTimeFrames = [(6.65,6.78), (7.85, 8.12), (9.95,10.15), (10.21, 10.40), (12.32,12.47)]
#peakTimeFrames = [(6.65,6.78), (7.86, 8.12), (9.95,10.15), (12.32,12.47)]
#peakTimeFrames = [(9.96,10.15)]
omitRatios = ['15N/13C', '13C/15N']
#omitRatios = []
#peaks = _importPeaksFromFTStatFile(inputStandardFile)
#pandas = _convertToPandasDataFrame(peaks)f
#Merged = _combineSubstituted(pandas, None, gc_elution_on, peakTimeFrames, 2, isotopeList, 0.10, outputPath=None)
#Output = _calcRawFileOutput(Merged, gc_elution_on, isotopeList, omitRatios)
#df = _convertDictToDF(Output)
#Output = _calcFolderOutput(inputStandardFolder, gc_elution_on,  peakTimeFrames,  isotopeList, omitRatios, outputPath)
Output, StatsOutput = DataAnalyzer.calc_Folder_Output(inputStandardFolder, cullOn=None, cullZeroScansOn=False, gcElutionOn=gc_elution_on, weightByNLHeight=True, gcElutionTimes = peakTimeFrames,  cullAmount=2, isotopeList = isotopeList, NL_over_TIC=0.10, omitRatios = omitRatios, fileCsvOutputPath=None)

