##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Modified: Tues July 21, 2020
@author: sarahzeichner

This code has all of the data processing code, to take data after it has been processed by FT statistic (or equivalent) and calculate isotope ratios based on input
"""

import matplotlib
import csv
import os
import numpy as np
import xlrd
import pandas as pd
import math 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import io
from collections import Counter

import MethodFile

def Process_RAW_File_JSON(methodFilePath, jsonPath):
    '''
    Process each JSON with raw file input

    #TODO: incoporate triageModule.py, and readAndRejectJSON.py

    '''
    #methodFile = MethodFile.MethodFile.ReadMethodFile(methodFilePath)


def import_Peaks_From_FTStatFile(inputFileName):
    '''
    Import peaks from FT statistic output file into a workable form, step 1
    
    Inputs:
        inputFileName: The excel file to input from
        
    Outputs:
        A list, containing dictionaries for each mass with a set of peaks in the excel file. The dictionaries have entries for 'tolerance', 'lastScan', 'refMass', and 'scans'. The 'scans' key directs to another list; this has a dictionary for each indvidual scan, giving a bunch of data about that scan. 
    '''

    # Open the worksheet with xlrd
    wb = xlrd.open_workbook(inputFileName)
    ws = wb.sheet_by_index(0)

    # list containing packets containing dicts of microscans for each measured peak
    peaks = []
    onMicroScans = False

    for r in range(ws.nrows)[2:]:
        if 'Tolerance:' in ws.row_values(r):
            peaks.append({})
            # Get the tolerance
            try:
                tol = float(ws.cell_value(r, 1).strip(' ppm'))
            except:
                tol = float(ws.cell_value(r, 1).strip(' mmu'))
            # Get the last scan of microscan packets
            lastScan = int(ws.cell_value(r, 7))
            # Get the ref mass
            refMass = float(ws.cell_value(r, 9))
            peaks[-1] = {'tolerance': tol, 'lastScan': lastScan,
                         'refMass': refMass, 'scans': []}
            continue

        if 'Measured Mass:' and 'Ret. Time:' in ws.row_values(r):
            # Saving rows to know what goes in each column (in case of changes later in the sheet)
            colIndex = ws.row_values(r)
            onMicroScans = True
            continue
        if onMicroScans:
            if 'Aver:' in ws.row_values(r):
                onMicroScans = False
                continue
            measuredMass = ws.cell_value(r, colIndex.index('Measured Mass:'))
            retTime = ws.cell_value(r, colIndex.index('Ret. Time:'))
            scanNumber = ws.cell_value(r, colIndex.index('Scan Number:'))
            absIntensity = ws.cell_value(r, colIndex.index('Abs. Intensity:'))
            integrationTime = ws.cell_value(r, colIndex.index('IT [ms]:'))
            ftResolution = ws.cell_value(r, colIndex.index('FT Resolution:'))
            peakNoise = ws.cell_value(r, colIndex.index('Peak Noise'))
            totalIonCount = ws.cell_value(r, colIndex.index('TIC:'))
            ticTimesIT = ws.cell_value(r, colIndex.index('TIC*IT:'))
            peakResolution = ws.cell_value(
                r, colIndex.index('Peak Resolution'))
            peakBaseline = ws.cell_value(r, colIndex.index('Peak Baseline'))
            if measuredMass == '' and retTime == '' and scanNumber == '':
                continue
            if scanNumber == lastScan:
                onMicroScans = False
            peaks[-1]['scans'].append(({'mass': measuredMass, 'retTime': retTime, 'tic': totalIonCount,
                                        'scanNumber': scanNumber, 'absIntensity': absIntensity, 'integTime': integrationTime,'TIC*IT': ticTimesIT,'ftRes': ftResolution, 'peakNoise': peakNoise, 'peakRes': peakResolution, 'peakBase': peakBaseline}))
    return(peaks)

def convert_To_Pandas_DataFrame(peaks):
    '''
    Import peaks from FT statistic output file into a workable form, step 2
    
    Inputs:
        peaks: The peaks output from _importPeaksFromFTStatistic; a list of dictionaries. 
        
    Outputs: 
        A list, where each element is a pandas dataframe for an individual peak extracted by FTStatistic (i.e. a single line in the FTStat input .txt file). 
    '''
    rtnAllPeakDF = []

    for peak in peaks:
        try:
            columnLabels = list(peak['scans'][0])
            data = np.zeros((len(peak['scans']), len(columnLabels)))
        except:
            print("Could not find peak " + str(peak))
            continue
        # putting all scan data into an array
        for i in range(len(peak['scans'])):
            for j in range(len(columnLabels)):
                data[i, j] = peak['scans'][i][columnLabels[j]]
        # scan numbers as separate array for indices
        scanNumbers = data[:, columnLabels.index('scanNumber')]
        # constructing data frame
        peakDF = pd.DataFrame(data, index=scanNumbers, columns=columnLabels)

        # calculate counts and add to the dataframe
        peakDF = calculate_Counts_And_ShotNoise(peakDF)

        # add it to the return pandas DF
        rtnAllPeakDF.append(peakDF)

    return(rtnAllPeakDF)

def calculate_Counts_And_ShotNoise(peakDF,CN=4.4,z=1,RN=120000,Microscans=1):
    '''
    Calculate counts of each scan peak
    
    Inputs: 
        peakDf: An individual dataframe consisting of a single peak extracted by FTStatistic.
        CN: A factor from the 2017 paper to convert intensities into counts
        RN: A reference resolution, for the same purpose (do not change!)
        z: The charge of the ion, used to convert intensities into counts
        Microscans: The number of scans a cycle is divided into, typically 1.
        
    Outputs: 
        The inputDF, with a column for 'counts' added. 
    '''
    peakDF['counts'] = (peakDF['absIntensity'] /
                        peakDF['peakNoise']) * (CN/z) *(RN/peakDF['ftRes'])**(0.5) * Microscans**(0.5)
    return peakDF

def calc_Append_Ratios(singleDf, allBelowOne = True, isotopeList = ['13C','15N','UnSub']):
    '''
    Calculates both 15N and 13C ratios, writes them such that they are < 1, and adds them to the dataframe.
    Inputs:                               
            singleDF: An individual pandas dataframe, consisting of multiple peaks from FTStat combined into one dataframe by the _combinedSubstituted function.
            allBelowOne: if True, outputs ratios as 'Sub/unSub' or 'unSub/Sub', whichever is below 1. If false, outputs
            all as 'sub/unSub'. 
            isotopeList: A list of isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. This must be the same for each fragment. This is used to determine all ratios of interest, i.e. 13C/UnSub, and label them in the proper order. 
            
    Outputs:
            The dataframe with ratios added. It computes all ratios, because why not. 
    '''
    if allBelowOne:
        for i in range(len(isotopeList)):
            for j in range(len(isotopeList)):
                if j>i:
                    #determine which has more counts, set ratio accordingly
                    if singleDf['counts' + isotopeList[i]].sum() <= singleDf['counts' + isotopeList[j]].sum():
                        singleDf[isotopeList[i] + '/' + isotopeList[j]] = singleDf['counts' + isotopeList[i]] / singleDf['counts' + isotopeList[j]]
                    else:
                        singleDf[isotopeList[j] + '/' + isotopeList[i]] = singleDf['counts' + isotopeList[j]] / singleDf['counts' + isotopeList[i]]
         
    else:
        for i in range(len(isotopeList)):
            for j in range(len(isotopeList)):
                if j>i:
                    singleDf[isotopeList[i] + '/' + isotopeList[j]] = singleDf['counts' + isotopeList[i]] / singleDf['counts' + isotopeList[j]]
    return singleDf

def combine_Substituted_Peaks(peakDF, cullOn = [], cullZeroScansOn = False, gc_elution_on = False, gc_elution_times = [], cullAmount = 2, fragmentIsotopeList = [['13C','15N','UnSub']], NL_over_TIC = 0.10, csv_output_path=None):
    '''
    Merge all extracted peaks from a given fragment into a single dataframe. For example, if I extracted six peaks, the 13C, 15N, and unsubstituted of fragments at 119 and 109, this would input a list of six dataframes (one per peak) and combine them into two dataframes (one per fragment), each including data from the 13C, 15N, and unsubstituted peaks of that fragment.
    
    Inputs: 
        peakDF: A list of dataframes. The list is the output of the _convertToPandasDataFrame function, and containts
        an individual dataframe for each peak extracted with FTStatistic. 
        cullOn: A target variable, like 'tic', or 'TIC*IT' to use to determine which scans to call. 
        cullZeroScansOn: Toggle whether or not you want to cull out zero scan counts.
        gc_elution_on: Set to True if you need to integrate over GC curve, and account for change in ion NL over time
        gc_elution_times: time frames that each peak elutes at, to cull for
        cullAmount: A number of standard deviations from the mean. If an individual scan has the cullOn variable outside of this range, culls the scan; i.e. if cullOn is 'TIC*IT' and cullAmount is 3, culls scans where TIC*IT is more than 3 standard deviations from its mean. 
        NL_over_TIC: specific NL/TIC that designates what a "peak" should look like. default 0.1, currently not implemented
        fragmentIsotopeList: A list of lists, where each interior list corresponds to a peak and gives the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
       
    Outputs: 
        A list of combined dataframes; in the 119/109 example above, it will output a list of two dataframes, [119, 109]
    where each dataframe combines the information for each substituted peak. This allows one to cull on the basis of different
    inputs (i.e. TIC*IT) as well as compute ratios (both of which are done by this function and _calcAppendRatios, which this 
    function calls). 
    '''
    DFList = []
    peakIndex = 0
    
    for fragment in fragmentIsotopeList:
        #First substitution, keep track of TIC*IT etc from here
        sub = fragment[0]
        df1 = peakDF[peakIndex].copy()

        #Rename columns to keep track of them
        df1.rename(columns={'mass':'mass'+sub,'counts':'counts'+sub,'absIntensity':'absIntensity'+sub,
                                    'peakNoise':'peakNoise'+sub},inplace=True)

        #helper variable to assign labels to final dataframe
        isotopeListIndex = 0

        #set start and end indices to  defaults based on what is already in the code,  so if errors, it just takes all  data
        start_index = int(df1.first_valid_index())
        end_index = int(df1.last_valid_index())
        
        #Set up a column to track total NL of peaks of fragment of interest for GC elution
        #and set up parameters for this specific fragment elution
        if gc_elution_on == True:
            thisGCElutionTimeRange = gc_elution_times[int(peakIndex / numberPeaksPerFragment)]
            df1['sumAbsIntensity'] = df1['absIntensity'+sub]
            
       #add additional dataframes
        for additionalPeakIdx in range(len(fragment)-1):
            df2 = peakDF[peakIndex + additionalPeakIdx+1].copy()

            sub = fragment[additionalPeakIdx+1]
            df2.rename(columns={'mass':'mass'+sub,'counts':'counts'+sub,'absIntensity':'absIntensity'+sub,
                            'peakNoise':'peakNoise'+sub},inplace=True)

            #Drop duplicate information
            df2.drop(['retTime','tic','integTime','TIC*IT','ftRes','peakRes','peakBase'],axis=1,inplace=True) 

            # merge into other dataFrame
            df1 = pd.merge_ordered(df1, df2,on='scanNumber',suffixes =(False,False))

            isotopeListIndex += 1
            
            if gc_elution_on == True:
                df1['sumAbsIntensity'] = df1['sumAbsIntensity'] + df2['absIntensity'+sub]

        #Checks each peak for values which were not recorded (e.g. due to low intensity) and fills in zeroes
        #I think this accomplishes the same thing as the zeroFilling in FTStat
        for string in fragment:
            df1.loc[df1['mass' + string].isnull(), 'mass' + string] = 0
            df1.loc[df1['absIntensity' + string].isnull(), 'absIntensity' + string] = 0
            df1.loc[df1['peakNoise' + string].isnull(), 'peakNoise' + string] = 0
            df1.loc[df1['counts' + string].isnull(), 'counts' + string] = 0 

        #Cull zero scans
        if cullZeroScansOn == True:
            df1 = cull_Zero_Scans(df1)

        #Cull based on time frame for GC peaks
        if gc_elution_on == True and gc_elution_times != 0:
            start_index, end_index = cull_On_GC_Peaks(df1, thisGCElutionTimeRange, NL_over_TIC)
        df1 = df1[start_index:end_index]

        #Calculates ratio values and adds them to the dataframe. Weighted averages will be calculated in the next step
        df1 = calc_Append_Ratios(df1,isotopeList = fragment)
        #Given a key in the dataframe, culls scans outside specified multiple of standard deviation from the mean
        if cullOn != None:
            if cullOn not in list(df1):
                raise Exception('Invalid Cull Input')
            maxAllowed = df1[cullOn].mean() + cullAmount * df1[cullOn].std()
            minAllowed = df1[cullOn].mean() - cullAmount * df1[cullOn].std()

            df1 = df1.drop(df1[(df1[cullOn] < minAllowed) | (df1[cullOn] > maxAllowed)].index)

        peakIndex += len(fragment)
        #Adds the combined dataframe to the output list
        DFList.append(df1)
        #Test by writing to CSV, in case you want to check the output
        if csv_output_path != None:
            df1.to_csv(csv_output_path, index=True, header=True)
 
        else:
            pass
    return DFList

def cull_Zero_Scans(df):
    '''
    Inputs:
        df: input dataframe to cull
    Outputs:
        culled df without zero
    '''
    indexNames = df[(df['counts'] == 0)].index
    df = df[~(df == 0).any(axis=1)]
    return df


def cull_On_GC_Peaks(df, gcElutionTimeFrame = (0,0), NL_over_TIC=0.1):
    '''
    Inputs: 
        df: input dataframe to cull
        gcElutionTimeFrame: elution of gc peaks, currently specified by the user
        NL_over_TIC: specific NL/TIC that designates what a "peak" should look like. default 0.1, currently not implemented
    Outputs: 
       culled df based on input elution times for the peaks
    '''
    # get the scan numbers for the retention  time frame
    if gcElutionTimeFrame != (0,0):
        start_index = df.loc[df['retTime'] == gcElutionTimeFrame[0]].index.values.astype(int)[0]
        end_index = df.loc[df['retTime'] == gcElutionTimeFrame[1]].index.values.astype(int)[0]
    
    return start_index, end_index
    
def calc_Raw_File_Output(df, weightByNLHeight=False, isotopeList = ['13C','15N','UnSub'],omitRatios = []):
    '''
    For each ratio of interest, calculates mean, stdev, SErr, RSE, and ShotNoise based on counts. Outputs these in a dictionary which organizes by fragment (i.e different entries for fragments at 119 and 109).
    
    Inputs:
        df: A merged data frames from the _combineSubstituted function corresponding to a single fragment. 
        weightByNLHeight: Specify whether you want to calculate weighted averages by NL score
        isotopeList: A list of isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. This must be the same for each fragment. This is used to determine all ratios of interest, i.e. 13C/UnSub, and label them in the proper order. 
        omitRatios: A list of ratios to ignore. I.e. by default, the script will report 13C/15N ratios, which one may not care about. In this case, the list should be ['13C/15N','15N/13C'], including both versions, to avoid errors. 
         
    Outputs: 
        A dictionary giving mean, stdev, StandardError, relative standard error, and shot noise limit for all peaks.  
    '''
    #Initialize output dictionary 
    rtnDict = {}
      
    #Adds the peak mass to the output dictionary
    key = df.keys()[0]
    massStr = str(round(df[key].median(),1))
        
    rtnDict[massStr] = {}
        
    for i in range(len(isotopeList)):
        for j in range(len(isotopeList)):
            if j>i:
                if isotopeList[i] + '/' + isotopeList[j] in df:
                    header = isotopeList[i] + '/' + isotopeList[j]
                else:
                    try:
                        header = isotopeList[j] + '/' + isotopeList[i]
                    except:
                        raise Exception('Sorry, cannot find ratios for your input isotopeList')

                if header in omitRatios:
                    continue
                else: 
                    if weightByNLHeight==True:
                        #Weight based on NL value for elution 
                        rtnDict[massStr][header] = {}
                        values = df[header]
                        weights = df['absIntensityUnSub']
                        average = np.average(values, weights=weights)
                        rtnDict[massStr][header]['Ratio'] = average
                        rtnDict[massStr][header]['StDev'] = math.sqrt(np.average((values-average)**2, weights=weights))
                        rtnDict[massStr][header]['StError'] = rtnDict[massStr][header]['StDev'] / np.power(len(df),0.5)
                        rtnDict[massStr][header]['RelStError'] = rtnDict[massStr][header]['StError'] / rtnDict[massStr][header]['Ratio']
                        rtnDict[massStr][header]['ShotNoiseLimit by Quadrature'] = (1/df['counts' + isotopeList[i]].sum() +1/df['counts' + isotopeList[j]].sum())**(1/2)

                    else:
                        #perform calculations and add them to the dictionary     
                        rtnDict[massStr][header] = {}
                        rtnDict[massStr][header]['Ratio'] = np.mean(df[header])
                        rtnDict[massStr][header]['StDev'] = np.std(df[header])
                        rtnDict[massStr][header]['StError'] = rtnDict[massStr][header]['StDev'] / np.power(len(df),0.5)
                        rtnDict[massStr][header]['RelStError'] = rtnDict[massStr][header]['StError'] / rtnDict[massStr][header]['Ratio']
                        rtnDict[massStr][header]['ShotNoiseLimit by Quadrature'] = (1/df['counts' + isotopeList[i]].sum() +1/df['counts' + isotopeList[j]].sum())**(1/2)

                        averageTIC = np.mean(df['tic'])
                        valuesTIC = df['tic']
                        rtnDict[massStr][header]['TICVar'] =  math.sqrt(np.mean((valuesTIC-averageTIC)**2))/np.mean(valuesTIC)

                        averageTICIT = np.mean(df['TIC*IT'])
                        valuesTICIT = df['TIC*IT']
                        rtnDict[massStr][header]['TIC*ITMean'] = averageTICIT
                        rtnDict[massStr][header]['TIC*ITVar'] =  math.sqrt(np.mean((valuesTICIT-averageTICIT)**2))/np.mean(valuesTICIT)

    return rtnDict

def calcOutputList(Merged, fragmentIsotopeList, fragmentMostAbundant):
    '''
    For all peaks in the input file, calculates results via calc_Raw_File_Output and adds these results to a list. Outputs the final list. 
    
    Inputs:
        Merged: The list containing all merged data frames from the _combineSubstituted function. 
        fragmentIsotopeList: A list of lists, where each interior list corresponds to a peak and gives the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        fragmentMostAbundant: A list, where each entry is the most abundant isotope in a fragment. The order of fragments should correspond to the order given in "fragmentIsotopeList".  
    
    Outputs:
        A list of dictionaries. Each dictionary has a single key value pair, where the key is the identity of the fragment and the value is a dictionary. The value dictionary has keys of isotope ratios (e.g. "D/13C") keyed to dictionaries giving information about that ratio measurement. 
        
    Future: Maybe rethink this as outputting a dictionary rather than a list, which may be cleaner? But outputting as a list keeps the same ordering as the original Merged list, which I like. 
    '''
    
    outputList = []
    for fIdx, fragment in enumerate(fragmentIsotopeList):
        mostAbundant = fragmentMostAbundant[fIdx]

        perms = []
        for x in fragment:
            for y in fragment:
                perms.append(x + '/' + y)

        omitRatios = [x for x in perms if mostAbundant not in x.split('/')]
        output = calc_Raw_File_Output(Merged[fIdx],isotopeList = fragment,omitRatios = omitRatios)

        outputList.append(output)
        
    return outputList

def calc_Folder_Output(folderPath, cullOn=None, cullZeroScansOn=False, gcElutionOn=False, weightByNLHeight=False, gcElutionTimes = [],  cullAmount=2, fragmentIsotopeList = [['13C','15N','UnSub']], fragmentMostAbundant = ['UnSub'], NL_over_TIC=0.10, fileCsvOutputPath=None):
    '''
    For each raw file in a folder, calculate mean, stdev, SErr, RSE, and ShotNoise based on counts. Outputs these in a dictionary which organizes by fragment (i.e different entries for fragments at 119 and 109).  
    Inputs:
        folderPath: Path that all the .xslx raw files are in. Files must be in this format to be processed.
        cullOn: cull specific range of scans
        cullZeroScansOn: toggle to eliminate any scans with zero counts
        gcElutionOn: Specify whether you expect elution to change over time, so that you can calculate weighted averages
        weightByNLNeight: Toggle "TRUE" to weight R values by height of NL score
        gcElutionTimes: Time frames to cull the GC peaks for
        cullAmount: If you pass in a range of scans for "cullOn", this is the number of stddevs beyond which scans are culled. default is 2.
        fragmentIsotopeList: A list of lists, where each interior list corresponds to a peak and gives the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        fragmentMostAbundant: A list, where each entry is the most abundant isotope in a fragment. The order of fragments should correspond to the order given in "fragmentIsotopeList".  
        NL_over_TIC: Currently not used. Idea would be to cull any scans that the peak is below this percent of total TIC.
        fileCSVOutputPath: path name if you want to output each file as you process
        
    Outputs: 
        Output is a tuple:
        A dataframe giving mean, stdev, standardError, relative standard error, and shot noise limit for all peaks. 
        A dataframe with calculated statistics (average, stddev, stderror and rel std error)
        (Both the dataframes are also exported as csvs to the original input folder)
    '''

    ratio = "Ratio"
    stdev = "StdDev"
    rtnAllFilesDF = []
    header = ["FileNumber", "Fragment", "IsotopeRatio", "Average", "StdDev", "StdError", "RelStdError","TICVar","TIC*ITVar","TIC*ITMean"]
    #get all the file names in the folder with the same end 
    fileNames = [x for x in os.listdir(folderPath) if x.endswith(".xlsx")]
    peakNumber = 0

    #Process through each raw file added and calculate statistics for fragments of interest
    for i in range(len(fileNames)):
        thisFileName = str(folderPath + '/' + fileNames[i])
        print(thisFileName) #for debugging
        thesePeaks = import_Peaks_From_FTStatFile(thisFileName)
        thisPandas = convert_To_Pandas_DataFrame(thesePeaks)
        thisMergedDF = combine_Substituted_Peaks(peakDF=thisPandas,cullOn=cullOn, cullZeroScansOn = cullZeroScansOn, gc_elution_on=gcElutionOn, gc_elution_times=gcElutionTimes, cullAmount=cullAmount, fragmentIsotopeList=fragmentIsotopeList, NL_over_TIC=NL_over_TIC, csv_output_path=fileCsvOutputPath)
        thisOutputList = calcOutputList(thisMergedDF, fragmentIsotopeList, fragmentMostAbundant)

        for fragmentOutput in thisOutputList:
            #there will be just one key, value pair in the "fragmentOutput" dictionary, keying the name of that fragment (e.g. "133") to information about the ratios for that fragment
            fragment = list(fragmentOutput.keys())[0]
            isotopeRatios = fragmentOutput[fragment]

            for ratioName, ratioData in isotopeRatios.items():
                thisPeak = fragment
                thisRatio = ratioName
                #add subkey to each separate df for isotope specific 
                thisRVal = ratioData["Ratio"]
                thisStdDev = ratioData["StDev"]
                thisStError = ratioData["StError"] 
                thisRelStError = ratioData["RelStError"]
                thisTICVar = ratioData["TICVar"] 
                thisTICITVar = ratioData["TIC*ITVar"]
                thisTICITMean = ratioData["TIC*ITMean"]
                thisRow = [thisFileName, thisPeak, thisRatio, thisRVal, thisStdDev, thisStError,thisRelStError,thisTICVar,thisTICITVar,thisTICITMean] 
                rtnAllFilesDF.append(thisRow)

    rtnAllFilesDF = pd.DataFrame(rtnAllFilesDF)
    # set the header row as the df header
    rtnAllFilesDF.columns = header 

    #sort by fragment and isotope ratio, output to csv
    rtnAllFilesDF = rtnAllFilesDF.sort_values(by=['Fragment', 'IsotopeRatio'], axis=0, ascending=True)
    rtnAllFilesDF.to_csv(str(folderPath + '/' + "all_data_output.csv"), index = False, header=True)

    #we can now calculate average, stdev, relstdev for each fragment across replicate measurements 
    if len(fileNames)>1: #only calculate  stats if there is more than one file
        avgDF = rtnAllFilesDF.groupby(['Fragment', 'IsotopeRatio'])["Average"].mean()
        countDF = rtnAllFilesDF.groupby(['Fragment', 'IsotopeRatio'])["Average"].count()
        stdDF = rtnAllFilesDF.groupby(['Fragment', 'IsotopeRatio'])["Average"].std()
        sqrtCountDF = np.power(countDF, 0.5)
        stdErrorDF = np.divide(stdDF, sqrtCountDF)
        relStdErrorDF = np.divide(stdErrorDF, avgDF)

    statsDF = pd.DataFrame([avgDF, countDF, stdDF, stdErrorDF, relStdErrorDF], index=["Avg R Val", "N", "StdDev", "StdError", "RelStdError"]) 

    #statsDF.rename(index={0:'IsotopeRatio',1:'Average R Val',2:'n', 3:'StdDev', 4:'StdError', 5:'RelStdError'}, inplace=True)
    statsDF.to_csv(str(folderPath + '/' + "stats_output.csv"), index = True, header=True)
    #output results to csv
    return rtnAllFilesDF, statsDF
           
