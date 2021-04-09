import itertools
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
  
def plotTicStats(Merged, figsize = (16,10)):
    '''
    A kind of basic diagnostic output plot, especially useful for LC measurements which should have a constant TIC with time. This shows "TIC", "IT", and "TIC*IT" as a function of retTime, then plots all as functions of each other. It also prints the mean and STDEV of each. The most important "good/bad" criteria we've seen is that STDEV of TIC*IT should be <= 0.10 for good quality data. 
    
    Written idiosyncratically compared to other functions in this document, using gridspec rather than plt.subplots. If one wants to plot more basic measures as functions of each other, one can add more tuples to "plotList". 
    
    Note that with multiple fragments, this will only plot the information from the first fragment. Typically the baseline measures will all be the same so this is ok. But in a weird case where it is not true, this should be evaulated. 
    
    Inputs:
        Merged: The list containing all merged data frames from the _combineSubstituted function. 
        figsize: Size of output figure
    
    Outputs:
        None. Displays figure. 
    
    '''
    #Intiialize figure
    fig = plt.figure(figsize=figsize)

    #Initializes an outer grid. Sets of subplots for each peak will be placed into this outer grid. 
    numberOfPeaks = 1
    outer_grid = gridspec.GridSpec(1, 1, wspace=0.5, hspace=0.25)

    #The base case gives the tic plots for the whole measurement, then the loop plots each fragment
    for i in range(numberOfPeaks):
        if i == 0:
            inner_grid = gridspec.GridSpecFromSubplotSpec(2, 3,
                subplot_spec=outer_grid[i], wspace=0.25, hspace=0.40)

            df = Merged[0]
            plotlist = [('retTime','tic'),('retTime','TIC*IT'),('retTime','integTime'),('integTime','tic'),('tic','TIC*IT'),('integTime','TIC*IT')]

            for i, item in enumerate(plotlist):
                ax = plt.Subplot(fig, inner_grid[plotlist.index(item)])
                if i <= 2:
                    ax.plot(df[item[0]].tolist(),df[item[1]].tolist())
                else:
                    ax.scatter(df[item[0]].tolist(),df[item[1]].tolist())
                    
                ax.set_xlabel(item[0])
                ax.set_ylabel(item[1])
                ax.set_title(item[1] + ' versus ' + item[0])
                fig.add_subplot(ax)

        else:
            pass

    print('TIC Mean and STD')
    print('{:.2e}'.format(Merged[0]['tic'].mean()))
    print(round(Merged[0]['tic'].std()/Merged[0]['tic'].mean(),2))
    print('IT Mean and STD')
    print('{:.2e}'.format(Merged[0]['integTime'].mean()))
    print(round(Merged[0]['integTime'].std()/Merged[0]['integTime'].mean(),2))
    print('TIC*IT Mean and STD')
    print('{:.2e}'.format(Merged[0]['TIC*IT'].mean()))
    print(round(Merged[0]['TIC*IT'].std()/Merged[0]['TIC*IT'].mean(),2))
    
def plotCounts(Merged, isotopeList, omitRatios=[], figsize = (14,16), CPmS = False, printZeros = False):
    '''
    For each fragment, computes the counts observed per scan for each isotope of that fragment. For each fragment, creates two output plots, one showing counts per scan on an automatically scaled y-axis, the other on a y-axis scaled from 0 to 200. The second plot allows one to observe low abundance isotopes in more detail. 

    Optionally, one may instead compute and display counts per millisecond. This measure may be more useful to assess the benefits of a cleaning procedure, for example. 

    Inputs:
        Merged: The list containing all merged data frames from the _combineSubstituted function. 
        isotopeList: A list of the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        figsize: Size of output figure
        CPmS: If true, calculate and output "counts per millisecond" rather than counts per scan. 
        printZeros: prints how many zero scans there are for each peak. 

    Outputs:
        None. Displays Figure. 
    '''
    l = len(isotopeList)
    markerList = itertools.cycle(('s', '^', 'o')) 
    nrows = l//2 + l % 2
    
    #Prevent weirdness with one dimensional arrays
    oneD = False
    if nrows == 1:
        oneD = True
        
    fig, ax = plt.subplots(nrows = nrows, ncols = 4,figsize = figsize)

    #TODO: fix the looping of the fragment mass names and addition to the plot
    for fragmentIndex in range(len(Merged)):

        #Adds the peak mass to the output dictionary
        key = round(Merged[fragmentIndex]['massUnSub'].mean())
        massStr = str(key)

        for i in range(len(isotopeList)):
            for j in range(len(isotopeList)):
                if j>i:
                    if isotopeList[i] + '/' + isotopeList[j] in Merged[fragmentIndex]:
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
                        markerStyle = next(markerList)
                        unsubCounts = np.array(Merged[fragmentIndex]['counts' + isotopeList[i]])
                        subCounts = np.array(Merged[fragmentIndex]['counts' + isotopeList[j]])
                        if CPmS == True:
                            unsubCounts /= np.array(Merged[fragmentIndex]['integTime'])
                            subCounts /= np.array(Merged[fragmentIndex]['integTime'])
                            ylabel = "Counts Per Millisecond"  
                        else:
                            ylabel = "Counts Per Scan"
                
                        fig.suptitle(ylabel)
                
                        if oneD == True:
                            accessSubplot = 2 * (fragmentIndex % 2)
                        else: 
                            accessSubplot = (fragmentIndex // 2, 2* (fragmentIndex % 2))
            
                        cAx = ax[accessSubplot]
                        cAx.scatter(Merged[fragmentIndex]['scanNumber'],unsubCounts,label = massStr+ " " + str(isotopeList[i]),marker = markerStyle,alpha = 0.4)
                        cAx.scatter(Merged[fragmentIndex]['scanNumber'],subCounts,label = massStr+ " " + str(isotopeList[j]),marker = markerStyle,alpha = 0.4)
                        
                        cAx.legend()
                        cAx.set_title(massStr)
                        cAx.set_ylabel(ylabel)
                        cAx.set_xlabel("Scan Number")
            
                        if oneD == True:
                            accessSubplot = 2 * (fragmentIndex % 2) + 1
                        else: 
                            accessSubplot = (fragmentIndex // 2, 2* (fragmentIndex % 2) + 1)
                
                        cAx = ax[accessSubplot]
                        cAx.scatter(Merged[fragmentIndex]['scanNumber'],unsubCounts,label = massStr+ " " + str(isotopeList[i]),marker = markerStyle,alpha = 0.4)
                        cAx.scatter(Merged[fragmentIndex]['scanNumber'],subCounts,label = massStr+ " " + str(isotopeList[j]),marker = markerStyle,alpha = 0.4)
                        cAx.set_ylim(0,200)
            
                        if printZeros:
                            print(str(isotopeList[j])+ ' Zero Counts: ' + str(len(unsubCounts) - np.count_nonzero(unsubCounts)))

            plt.tight_layout()

def plotOutputList(outputList, isotopeList,omitRatios, figsize = (18,16)):
    '''
    For each fragment, creates two plots. The first plot displays the precision of each isotope ratio of the fragment. The second shows the ratio of precision to shot noise for each isotope ratio. A reasonable cut-off for "good" quality data is RSE/SN < 2, so we draw a red line at 2. 
    
    Inputs:
        outputList: A dictionary, where the key is the identity of the fragment and the value is a list of dictionaryies of run properties 
        isotopeList: A list of the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        figsize: Size of output figure
        
    Outputs:
        None. Displays plot. 
    '''
    l = len(isotopeList)
    nrows = l//2 + l % 2
    
    #Prevent weirdness with one dimensional arrays
    oneD = False
    if nrows == 1:
        oneD = True
    fig, ax = plt.subplots(nrows = nrows, ncols = 4,figsize = figsize)

    keys = list(outputList.keys())
    for fragmentIndex in range(len(outputList)):
        #Adds the peak mass to the output dictionary
        massStr = keys[fragmentIndex]
        RSESN = []
        RSE = []
        #xlabels = []
        for i in range(len(isotopeList)):
            for j in range(len(isotopeList)):
                if j>i:
                    if isotopeList[i] + '/' + isotopeList[j] in outputList[massStr]:
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
                        RSE.append(outputList[massStr][header]['RelStError'])
                        RSESN.append(outputList[massStr][header]['RelStError']/outputList[massStr][header]['ShotNoiseLimit by Quadrature'])
       
                        if oneD == True:
                            accessSubplot = 2* (fragmentIndex % 2)
                        else: 
                            accessSubplot = (fragmentIndex // 2, 2* (fragmentIndex % 2))

                        cAx = ax[accessSubplot]
                        cAx.scatter(range(len(RSE)),RSE)
                        cAx.set_xticks(range(len(RSE)))
                        #cAx.set_xticklabels(xlabels,rotation = 40)
                        cAx.set_title('RelStError ' + massStr)
        
                        if oneD == True:
                            accessSubplot = 2* (fragmentIndex % 2) + 1
                        else: 
                            accessSubplot = (fragmentIndex // 2, 2* (fragmentIndex % 2) + 1)
        
                        cAx = ax[accessSubplot]
                        cAx.scatter(range(len(RSESN)),RSESN)
                        cAx.set_xticks(range(len(RSESN)))
                        #cAx.set_xticklabels(xlabels,rotation = 40)
                        cAx.hlines(2,0,len(RSESN)-1,colors='r',linestyle = '--')
                        cAx.set_title('RSE/SN ' + massStr)

    plt.tight_layout()
    
def plotIsotopeRatioVersusInstrumentalParams(Merged, specialIsotopeList, numBins=10, xvars=['retTime','tic','integTime','TIC*IT']):
    '''
    For each fragment, creates two plots. The first plot displays the precision of each isotope ratio of the fragment. The second shows the ratio of precision to shot noise for each isotope ratio. A reasonable cut-off for "good" quality data is RSE/SN < 2, so we draw a red line at 2. 
    
    Inputs:
        Merged: The list containing all merged data frames from the _combineSubstituted function. 
        isotopeList: A list of the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        figsize: Size of output figure
        
    Outputs:
        None. Displays plot. 
    '''

    #Initialize figure
    fig = plt.figure(figsize=(12,10))

    #Initializes an outer grid. Sets of subplots for each peak will be placed into this outer grid. 
    numberOfPeaks = 1
    outer_grid = gridspec.GridSpec(1, 1, wspace=0.5, hspace=0.25)

    count = Merged.shape[0]/numBins
    #The base case gives the tic plots for the whole measurement, then the loop plots each fragment
    for i in range(numberOfPeaks):
        if i == 0:
            for isotope in specialIsotopeList:
                for variable in xvars:
                    avgXList = []
                    avgYList = []
                    errorList = []

                for div in range(0,numBins):
                    if isotope + '/UnSub' in Merged:
                        targetVariable = isotope + '/UnSub'
                    else:
                        targetVariable = 'UnSub/' + isotope

                    sortDf = Merged.sort_values(by=[variable])
                    division = [int(div * count//1),int((div+1)*count//1)]
                    restrictDf = sortDf.iloc[list(range(division[0],division[1]))]
                    Ymean = np.mean(restrictDf[targetVariable])
                    YStd = np.std(restrictDf[targetVariable])
                    YStError = YStd / np.power(len(restrictDf),0.5)
                    Xmean = np.mean(restrictDf[variable])

                    avgXList.append(Xmean)
                    avgYList.append(Ymean)
                    errorList.append(YStError)

                plt.subplot(2,2,xvars.index(variable)+1)
                plt.errorbar(avgXList, avgYList, errorList, fmt='o')
                plt.xlabel(variable)
                plt.ylabel(targetVariable)
                plt.title(targetVariable + ' binned by ' + variable)

    plt.show()