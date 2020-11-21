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


def import_Peaks_From_FTStatFile(inputFileName):
    '''
    Import peaks from FT statistic output file into a workable form, s`tep 1
    
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

def combine_Substituted_Peaks(peakDF, cullOn = [], cullZeroScansOn = False, gc_elution_on = False, gc_elution_times = [], cullAmount = 2, isotopeList = ['13C','15N','UnSub'], NL_over_TIC = 0.10, csv_output_path=None):
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
        isotopeList: A list of isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. This must be the same for each fragment. This is used to determine all ratios of interest, i.e. 13C/UnSub, and label them in the proper order. 
       
    Outputs: 
        A list of combined dataframes; in the 119/109 example above, it will output a list of two dataframes, [119, 109]
    where each dataframe combines the information for each substituted peak. This allows one to cull on the basis of different
    inputs (i.e. TIC*IT) as well as compute ratios (both of which are done by this function and _calcAppendRatios, which this 
    function calls). 
    '''
    DFList = []
    numberPeaksPerFragment = len(isotopeList)

    for peakIndex in range(len(peakDF)):
        if peakIndex % numberPeaksPerFragment == 0:
            df1 = peakDF[peakIndex].copy()
            
            #set start and end indices to  defaults based on what is already in the code,  so if errors, it just takes all  data
            start_index = int(df1.first_valid_index())
            end_index = int(df1.last_valid_index())
            
            #Rename columns to keep track of them
            sub = isotopeList[0]
            df1.rename(columns={'mass':'mass'+sub,'counts':'counts'+sub,'absIntensity':'absIntensity'+sub,
                                'peakNoise':'peakNoise'+sub},inplace=True)

            #Set up a column to track total NL of peaks of fragment of interest for GC elution
            #and set up parameters for this specific fragment elution
            if gc_elution_on == True: #TODO: fix boolean processing
                thisGCElutionTimeRange = gc_elution_times[int(peakIndex / numberPeaksPerFragment)]
                df1['sumAbsIntensity'] = df1['absIntensity'+sub]
                       
            #helper variable to assign labels to final dataframe
            isotopeListIndex = 0

            #add additional dataframes
            for additionalDfIndex in range(numberPeaksPerFragment-1):
                df2 = peakDF[peakIndex + additionalDfIndex+1].copy()
                
                sub = isotopeList[additionalDfIndex+1]
                df2.rename(columns={'mass':'mass'+sub,'counts':'counts'+sub,'absIntensity':'absIntensity'+sub,
                                'peakNoise':'peakNoise'+sub},inplace=True)
                
                if gc_elution_on == True:
                    df1['sumAbsIntensity'] = df1['sumAbsIntensity'] + df2['absIntensity'+sub]
                
                #Drop duplicate information
                df2.drop(['retTime','tic','integTime','TIC*IT','ftRes','peakRes','peakBase'],axis=1,inplace=True) 
                          
                # merge 13C and 15N dataframes
                df1 = pd.merge_ordered(df1, df2,on='scanNumber',suffixes =(False,False))
                
                isotopeListIndex += 1

            #Checks each peak for values which were not recorded (e.g. due to low intensity) and fills in zeroes
            #I think this accomplishes the same thing as the zeroFilling in FTStat
            for string in isotopeList:
                df1.loc[df1['mass' + string].isnull(), 'mass' + string] = 0
                df1.loc[df1['absIntensity' + string].isnull(), 'absIntensity' + string] = 0
                df1.loc[df1['peakNoise' + string].isnull(), 'peakNoise' + string] = 0
                df1.loc[df1['counts' + string].isnull(), 'counts' + string] = 0 

            massStr = str(df1['mass'+isotopeList[0]].tolist()[0])

            #Cull zero scans
            if cullZeroScansOn == True:
                df1 = cull_Zero_Scans(df1)
            
            #Cull based on time frame for GC peaks
            if gc_elution_on == True and gc_elution_times != 0:
                start_index, end_index = cull_On_GC_Peaks(df1, thisGCElutionTimeRange, NL_over_TIC)
            df1 = df1[start_index:end_index]

            #Calculates ratio values and adds them to the dataframe. Weighted averages will be calculated in the next step
            df1 = calc_Append_Ratios(df1,isotopeList = isotopeList)
            #Given a key in the dataframe, culls scans outside specified multiple of standard deviation from the mean
            if cullOn != None:
                if cullOn not in list(df1):
                    raise Exception('Invalid Cull Input')
                maxAllowed = df1[cullOn].mean() + cullAmount * df1[cullOn].std()
                minAllowed = df1[cullOn].mean() - cullAmount * df1[cullOn].std()
                
                df1 = df1.drop(df1[(df1[cullOn] < minAllowed) | (df1[cullOn] > maxAllowed)].index)
            
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


    
def calc_Raw_File_Output(dfList, weightByNLHeight=False, isotopeList = ['13C','15N','UnSub'],omitRatios = []):
    '''
    For each ratio of interest, calculates mean, stdev, SErr, RSE, and ShotNoise based on counts. Outputs these in a dictionary which organizes by fragment (i.e different entries for fragments at 119 and 109).
    
    Inputs:
        dfList: A list of merged data frames from the _combineSubstituted function. Each dataframe constitutes one fragment.
        weightByNLHeight: Specify whether you want to calculate weighted averages by NL score
        isotopeList: A list of isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. This must be the same for each fragment. This is used to determine all ratios of interest, i.e. 13C/UnSub, and label them in the proper order. 
        omitRatios: A list of ratios to ignore. I.e. by default, the script will report 13C/15N ratios, which one may not care about. In this case, the list should be ['13C/15N','15N/13C'], including both versions, to avoid errors. 
         
    Outputs: 
        A dictionary giving mean, stdev, StandardError, relative standard error, and shot noise limit for all peaks.  
    '''
    #Initialize output dictionary 
    rtnDict = {}
      
    # iterate through each fragment
    for fragmentIndex in range(len(dfList)):
        
        #Adds the peak mass to the output dictionary
        key = round(dfList[fragmentIndex]['massUnSub'].mean())
        massStr = str(key)
        
        rtnDict[massStr] = {}
        
        for i in range(len(isotopeList)):
            for j in range(len(isotopeList)):
                if j>i:
                    if isotopeList[i] + '/' + isotopeList[j] in dfList[fragmentIndex]:
                        header = isotopeList[i] + '/' + isotopeList[j]
                    else:
                        try:
                            header = isotopeList[j] + '/' + isotopeList[i]
                        except:
                            raise Exception('Sorry, cannot find ratios for your input isotopeList')
                    
                    if header in omitRatios:
                        print("Ratios omitted:" + header)
                        continue
                    else: 
                        if weightByNLHeight==True:
                            #Weight based on NL value for elution 
                            rtnDict[massStr][header] = {}
                            values = dfList[fragmentIndex][header]
                            weights = dfList[fragmentIndex]['absIntensityUnSub']
                            average = np.average(values, weights=weights)
                            rtnDict[massStr][header]['Ratio'] = average
                            rtnDict[massStr][header]['StDev'] = math.sqrt(np.average((values-average)**2, weights=weights))
                            rtnDict[massStr][header]['StError'] = rtnDict[massStr][header]['StDev'] / np.power(len(dfList[fragmentIndex]),0.5)
                            rtnDict[massStr][header]['RelStError'] = rtnDict[massStr][header]['StError'] / rtnDict[massStr][header]['Ratio']
                            a = dfList[fragmentIndex]['counts' + isotopeList[i]].sum()
                            b = dfList[fragmentIndex]['counts' + isotopeList[j]].sum()
                            shotNoiseByQuad = np.power((1./a + 1./b),0.5)
                            rtnDict[massStr][header]['ShotNoiseLimit by Quadrature'] = shotNoiseByQuad
                            averageTIC = np.mean(dfList[fragmentIndex]['tic'])
                            valuesTIC = dfList[fragmentIndex]['tic']
                            rtnDict[massStr][header]['TICVar'] =  math.sqrt(np.mean((valuesTIC-averageTIC)**2))/np.mean(valuesTIC)
                            
                            averageTICIT = np.mean(dfList[fragmentIndex]['TIC*IT'])
                            valuesTICIT = dfList[fragmentIndex]['TIC*IT']
                            rtnDict[massStr][header]['TIC*ITMean'] = averageTICIT
                            rtnDict[massStr][header]['TIC*ITVar'] =  math.sqrt(np.mean((valuesTICIT-averageTICIT)**2))/np.mean(valuesTICIT)

    
                        else:
                            #perform calculations and add them to the dictionary     
                            rtnDict[massStr][header] = {}
                            rtnDict[massStr][header]['Ratio'] = np.mean(dfList[fragmentIndex][header])
                            rtnDict[massStr][header]['StDev'] = np.std(dfList[fragmentIndex][header])
                            rtnDict[massStr][header]['StError'] = rtnDict[massStr][header]['StDev'] / np.power(len(dfList[fragmentIndex]),0.5)
                            rtnDict[massStr][header]['RelStError'] = rtnDict[massStr][header]['StError'] / rtnDict[massStr][header]['Ratio']
                            a = dfList[fragmentIndex]['counts' + isotopeList[i]].sum()
                            b = dfList[fragmentIndex]['counts' + isotopeList[j]].sum()
                            shotNoiseByQuad = np.power((1./a + 1./b),0.5)
                            rtnDict[massStr][header]['ShotNoiseLimit by Quadrature'] = shotNoiseByQuad
                            
                            averageTIC = np.mean(dfList[fragmentIndex]['tic'])
                            valuesTIC = dfList[fragmentIndex]['tic']
                            rtnDict[massStr][header]['TICVar'] =  math.sqrt(np.mean((valuesTIC-averageTIC)**2))/np.mean(valuesTIC)
                            
                            averageTICIT = np.mean(dfList[fragmentIndex]['TIC*IT'])
                            valuesTICIT = dfList[fragmentIndex]['TIC*IT']
                            rtnDict[massStr][header]['TIC*ITMean'] = averageTICIT
                            rtnDict[massStr][header]['TIC*ITVar'] =  math.sqrt(np.mean((valuesTICIT-averageTICIT)**2))/np.mean(valuesTICIT)

    return rtnDict

def calc_Folder_Output(folderPath, cullOn=None, cullZeroScansOn=False, gcElutionOn=False, weightByNLHeight=False, gcElutionTimes = [],  cullAmount=2, isotopeList = ['UnSub', '15N', '13C'], NL_over_TIC=0.10, omitRatios = [], fileCsvOutputPath=None):
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
        isotopeList: A list of isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. This must be the same for each fragment. This is used to determine all ratios of interest, i.e. 13C/UnSub, and label them in the proper order. 
        NL_over_TIC: Currently not used. Idea would be to cull any scans that the peak is below this percent of total TIC.
        omitRatios: A list of ratios to ignore. I.e. by default, the script will report 13C/15N ratios, which one may not care about. In this case, the list should be ['13C/15N','15N/13C'], including both versions, to avoid errors.
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
    header = ["FileNumber", "Fragment", "IsotopeRatio", "Average", "StdDev", "StdError", "RelStdError","TICVar","TIC*ITVar","TIC*ITMean", 'ShotNoise']
    #get all the file names in the folder with the same end 
    fileNames = [x for x in os.listdir(folderPath) if x.endswith(".xlsx")]
    peakNumber = 0

    #Process through each raw file added and calculate statistics for fragments of interest
    for i in range(len(fileNames)):
        thisFileName = str(folderPath + '/' + fileNames[i])
        print(thisFileName) #for debugging
        thesePeaks = import_Peaks_From_FTStatFile(thisFileName)
        thisPandas = convert_To_Pandas_DataFrame(thesePeaks)
        thisMergedDF = combine_Substituted_Peaks(peakDF=thisPandas,cullOn=cullOn, cullZeroScansOn = cullZeroScansOn, gc_elution_on=gcElutionOn, gc_elution_times=gcElutionTimes, cullAmount=cullAmount, isotopeList=isotopeList, NL_over_TIC=NL_over_TIC, csv_output_path=fileCsvOutputPath)
        thisOutput = calc_Raw_File_Output(thisMergedDF, weightByNLHeight, isotopeList, omitRatios)
        keys = list(thisOutput.keys())
        peakNumber = len(keys)

        for peak in range(peakNumber):
            #key is each peak within the dictionary
            isotopeRatios = list(thisOutput[keys[peak]].keys())
            #subkey is each isotope ratio info for each peak
            for isotopeRatio in range(len(isotopeRatios)):
                thisPeak = keys[peak]
                thisRatio = isotopeRatios[isotopeRatio]
                #add subkey to each separate df for isotope specific 
                thisRVal = thisOutput[thisPeak][thisRatio]["Ratio"]
                thisStdDev = thisOutput[thisPeak][thisRatio]["StDev"]
                thisStError = thisOutput[thisPeak][thisRatio]["StError"] 
                thisRelStError = thisOutput[thisPeak][thisRatio]["RelStError"]
                thisTICVar = thisOutput[thisPeak][thisRatio]["TICVar"] 
                thisTICITVar = thisOutput[thisPeak][thisRatio]["TIC*ITVar"]
                thisTICITMean = thisOutput[thisPeak][thisRatio]["TIC*ITMean"]
                thisShotNoise = thisOutput[thisPeak][thisRatio]["ShotNoiseLimit by Quadrature"]
                thisRow = [thisFileName, thisPeak, thisRatio, thisRVal, thisStdDev, thisStError,thisRelStError,thisTICVar,thisTICITVar,thisTICITMean, thisShotNoise] 
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
           
def plot_Output(output,isotopeList = ['13C','15N','UnSub'],omitRatios = [],numCols = 2,widthMultiple = 4, heightMultiple = 4):
   #TODO: Fix this for the gc weighted average calculation
    '''
    Constructs a series of output plots for easy visualization
    
    Inputs: 
        output: The output dictionary from _calcOutput
        isotopeList: A list of isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. This must be the same for each fragment. This is used to determine all ratios of interest, i.e. 13C/UnSub, and label them in the proper order. 
        omitRatios: A list of ratios to ignore. I.e. by default, the script will report 13C/15N ratios, which one may not care about. In this case, the list should be ['13C/15N','15N/13C'], including both versions, to avoid errors. 
        numCols: The number of columns in the output plot
        widthMultiple: A factor which determines how wide the output plot is; higher is wider. This can be adjusted to make the plot appear nice, and may have to be changed based on how many ratios one is extracting. 
        heightMultiple: AS widthmultiple, but for height. 
        
    Outputs:
        None. Constructs and displays a plot visualizing the output data. 
    '''
    
    #Generate unique headers from list. Each must have two versions, normal and inverse. List of tuples
    headerList = []
    for i in range(len(isotopeList)):
        for j in range(len(isotopeList)):
            if j>i:
                normal = isotopeList[i] + '/' + isotopeList[j]
                inverse = isotopeList[j] + '/' + isotopeList[i]
                
                if normal not in omitRatios and inverse not in omitRatios:
                    headerList.append((normal,inverse))
    
    #Determine the number of plots based on the number of headers. Initialize figure
    numPlots = len(headerList) * 3
    numRows = numPlots // numCols + (numPlots % numCols > 0)
    fig, axes = plt.subplots(numRows, numCols, figsize = (numCols*widthMultiple, numRows*heightMultiple))
    row, col = 0, 0
    
    #For each header, iterate by fragments and pull information from output dictionary. 
    for header in headerList:
        PlotDict = {'Mass':[],'ShotNoise':[],'Error':[],'ErrorShotNoiseRat':[],'Ratio':[],'Type':[]}
        for fragment in output.items():
            if header[0] in fragment[1]:
                currentHeader = header[0]
            elif header[1] in fragment[1]:
                currentHeader = header[1]
            else:
                raise Exception('Current header is not in output' + str(header) + 'Try adding it to omitRatios')
                
            PlotDict['Mass'].append(fragment[0])
            PlotDict['ShotNoise'].append(fragment[1][currentHeader]['ShotNoiseLimit by Quadrature'])
            PlotDict['Error'].append(fragment[1][currentHeader]['RelStError'])
            PlotDict['ErrorShotNoiseRat'] = [a / b for a, b in zip(PlotDict['Error'], PlotDict['ShotNoise'])]
            PlotDict['Ratio'].append(fragment[1][currentHeader]['Ratio'])
            PlotDict['Type'].append(currentHeader)

        #Determine whether normal or inverse ratios are more common, and change all ratios to this format
        mode = Counter(PlotDict['Type'])
        modeStr = mode.most_common(1)[0][0]

        flippedRatios = []
        for index in range(len(PlotDict['Ratio'])):
            if PlotDict['Type'][index] == modeStr:
                flippedRatios.append(PlotDict['Ratio'][index])
            else:
                flippedRatios.append(1/PlotDict['Ratio'][index])
                print('Flipping Ratio' + str(PlotDict['Mass'][index]))
                
        xs = np.arange(len(PlotDict['Mass']))

        #Makes three plots for each isotope
        axes[row, col].set_xticks(xs)
        axes[row, col].set_xticklabels(PlotDict['Mass'])
        axes[row, col].scatter(xs, PlotDict['ShotNoise'], facecolors = 'none', edgecolors = 'k', label = 'Shot Noise')
        axes[row, col].scatter(xs, PlotDict['Error'], c='k', marker="o", label='Error')
        axes[row, col].legend()
        axes[row, col].set_ylim(0, axes[row, col].get_ylim()[1]*0.6)
        axes[row, col].set_title('Relative Standard Error vs Shot Noise, ' + modeStr)
        axes[row, col].set_xlabel('Fragment')
        axes[row, col].set_ylabel('Value')
        maxRSE = max(PlotDict['Error'])
        axes[row, col].set_ylim(0, 1.2 * maxRSE)

        if col == numCols - 1:
            row += 1
            col = 0
        else: 
            col += 1

        axes[row, col].set_xticks(xs)
        axes[row, col].set_xticklabels(PlotDict['Mass'])
        axes[row, col].scatter(xs, PlotDict['ErrorShotNoiseRat'], c= 'k', marker = 'o', label = 'Ratios')
        axes[row, col].axhline(y=2,color = 'r',linestyle = '--', label = 'Accuracy target')
        axes[row, col].legend()
        axes[row, col].set_title('Ratio RelStandardError to Shot Noise, ' + modeStr)
        axes[row, col].set_xlabel('Fragment')
        axes[row, col].set_ylabel('Ratio')

        if col == numCols - 1:
            row += 1
            col = 0
        else: 
            col += 1

        axes[row, col].set_xticks(xs)
        axes[row, col].set_xticklabels(PlotDict['Mass'])
        axes[row, col].scatter(xs, flippedRatios, c= 'k', marker = 'o', label = 'Ratios')
        axes[row, col].set_title('Ratios of ' + modeStr)
        axes[row, col].set_xlabel('Fragment')
        axes[row, col].set_ylabel('Ratio')

        if col == numCols - 1:
            row += 1
            col = 0
        else: 
            col += 1

    plt.tight_layout()
