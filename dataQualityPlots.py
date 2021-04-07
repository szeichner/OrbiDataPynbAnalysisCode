import itertools
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
  
def ticPlot(Merged, figsize = (16,10)):
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
    
def peakDriftPlot(Merged, fragmentIsotopeList, fragmentMostAbundant, figsize = (14,16), printMasses = False):
    '''
    For each peak in Merged, calculates the theoretical separation one should see between each isotope and the most abundant isotope of that fragment. Then calculates the actual separation one observes. Finds the difference between actual and theoretical. Plots the absolute value of this difference. As this difference gets larger, one should worry more about coalescence. 
    
    A reasonable first-order cutoff is a 0.0001 amu, so a line is plotted here. But this can and should be tuned in the future. Also note that a more useful plot may show the change in observed isotope ratio with respect to AGC, combining information from multiple measurements. We believe the isotope ratio may shift before any change in position is seen. 
    
    Inputs:
        Merged: The list containing all merged data frames from the _combineSubstituted function. 
        fragmentIsotopeList: A list of lists, where each interior list corresponds to a peak and gives the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        fragmentMostAbundant: A list, where each entry is the most abundant isotope in a fragment. The order of fragments should correspond to the order given in "fragmentIsotopeList".  
        figsize: Size of output figure
        printMasses: print the computed masses of each peak.
        
    Outputs:
        None. Displays figure. 
    '''
    #TODO: put these all lowercase, change everything to lowercase before using, to avoid errors.
    subMassDict = {'D':1.00627674587,'15N':0.997034886,'13C':1.003354835,
               'Unsub':0,'33S':0.999387735,'34S':1.995795825,'36S':3.995009525,
              '18O':2.0042449924,'17O':1.0042171364}
    
    l = len(fragmentIsotopeList)
    nrows = l//2 + l % 2
    
    #Prevent weirdness with one dimensional arrays
    oneD = False
    if nrows == 1:
        oneD = True
        
    fig, ax = plt.subplots(nrows = nrows, ncols = 2,figsize = figsize)


    for fIdx, fragment in enumerate(fragmentIsotopeList):
        mostAbundant = fragmentMostAbundant[fIdx]
        massList = []
        for i in fragment:
            mass = Merged[fIdx][Merged[fIdx]['mass' + i]!=0]['mass' + i].mean()
            massList.append(mass)

        diffIdent = []
        diffList = []
        mAIdx = fragment.index(mostAbundant)
        mAMass = massList[mAIdx]

        mAInc = 0
        mASubs = fragment[mAIdx].split('-')
        for sub in mASubs:
            mAInc += subMassDict[sub]


        for j, y in enumerate(massList):
            if j != mAIdx:

                mInc = 0
                mSubs = fragment[j].split('-')
                for sub in mSubs:
                    mInc += subMassDict[sub]

                massDiffActual = mAInc - mInc
                massDiffObserve = mAMass - y
                diffIdent.append(fragment[j]+ ' | ' + fragment[mAIdx])
                diffList.append(np.abs(massDiffObserve - massDiffActual))
            
            if printMasses:
                print(y)
        if oneD == True:
            accessSubplot = fIdx % 2
        else: 
            accessSubplot = (fIdx // 2, fIdx % 2)
            
        cAx = ax[accessSubplot]
        cAx.scatter(diffIdent,diffList)
        cAx.hlines(0.0001,-0.1,len(diffIdent)-0.9,colors='k',linestyle='--')
        cAx.set_title("Peak Drift " + str(round(massList[0],1)))
        cAx.set_xticklabels(diffIdent,rotation = 20)
        cAx.set_ylabel("Peak Drift in amu")

    plt.tight_layout()



def plotOutputList(outputList,  fragmentIsotopeList, fragmentMostAbundant, massStr, figsize = (18,16)):
    '''
    For each fragment, creates two plots. The first plot displays the precision of each isotope ratio of the fragment. The second shows the ratio of precision to shot noise for each isotope ratio. A reasonable cut-off for "good" quality data is RSE/SN < 2, so we draw a red line at 2. 
    
    Inputs:
        outputList: A list of dictionaries. Each dictionary has a single key value pair, where the key is the identity of the fragment and the value is a dictionary. The value dictionary has keys of isotope ratios (e.g. "D/13C") keyed to dictionaries giving information about that ratio measurement. 
        fragmentIsotopeList: A list of lists, where each interior list corresponds to a peak and gives the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        fragmentMostAbundant: A list, where each entry is the most abundant isotope in a fragment. The order of fragments should correspond to the order given in "fragmentIsotopeList".
        massStr: A string giving the mass of this fragment, e.g. "133".
        figsize: Size of output figure
        
    Outputs:
        None. Displays plot. 
    '''
    l = len(fragmentIsotopeList)
    nrows = l//2 + l % 2
    
    #Prevent weirdness with one dimensional arrays
    oneD = False
    if nrows == 1:
        oneD = True
    fig, ax = plt.subplots(nrows = nrows, ncols = 4,figsize = figsize)

    for fIdx, fragment in enumerate(fragmentIsotopeList):
        output = outputList[fIdx]
        RSESN = []
        RSE = []
        xlabels = []

        #Generally only one item, so one outer iteration
        for peak, data in output.items(): 
            for sub, measurement in data.items():

                xlabels.append(sub)
                RSE.append(measurement['RelStError'])
                RSESN.append(measurement['RelStError']/measurement['ShotNoiseLimit by Quadrature'])
        
        if oneD == True:
            accessSubplot = 2* (fIdx % 2)
        else: 
            accessSubplot = (fIdx // 2, 2* (fIdx % 2))

        cAx = ax[accessSubplot]
        cAx.scatter(range(len(RSE)),RSE)
        cAx.set_xticks(range(len(RSE)))
        cAx.set_xticklabels(xlabels,rotation = 40)
        cAx.set_title('RelStError ' + massStr[fIdx])
        
        if oneD == True:
            accessSubplot = 2* (fIdx % 2) + 1
        else: 
            accessSubplot = (fIdx // 2, 2* (fIdx % 2) + 1)
        
        cAx = ax[accessSubplot]
        cAx.scatter(range(len(RSESN)),RSESN)
        cAx.set_xticks(range(len(RSESN)))
        cAx.set_xticklabels(xlabels,rotation = 40)
        cAx.hlines(2,0,len(RSESN)-1,colors='r',linestyle = '--')
        cAx.set_title('RSE/SN ' + massStr[fIdx])

    plt.tight_layout()
    
def countsPlot(Merged, fragmentIsotopeList, fragmentMostAbundant, massStr, figsize = (14,16), CPmS = False, printZeros = False):
    '''
    For each fragment, computes the counts observed per scan for each isotope of that fragment. For each fragment, creates two output plots, one showing counts per scan on an automatically scaled y-axis, the other on a y-axis scaled from 0 to 200. The second plot allows one to observe low abundance isotopes in more detail. 

    Optionally, one may instead compute and display counts per millisecond. This measure may be more useful to assess the benefits of a cleaning procedure, for example. 

    Inputs:
        Merged: The list containing all merged data frames from the _combineSubstituted function. 
        fragmentIsotopeList: A list of lists, where each interior list corresponds to a peak and gives the isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
        fragmentMostAbundant: A list, where each entry is the most abundant isotope in a fragment. The order of fragments should correspond to the order given in "fragmentIsotopeList".
        massStr: A string giving the mass of this fragment, e.g. "133".
        figsize: Size of output figure
        CPmS: If true, calculate and output "counts per millisecond" rather than counts per scan. 
        printZeros: prints how many zero scans there are for each peak. 

    Outputs:
        None. Displays Figure. 
    '''
    l = len(fragmentIsotopeList)
    markerList = itertools.cycle(('s', '^', 'o')) 
    nrows = l//2 + l % 2
    
    #Prevent weirdness with one dimensional arrays
    oneD = False
    if nrows == 1:
        oneD = True
        
    fig, ax = plt.subplots(nrows = nrows, ncols = 4,figsize = figsize)


    for fIdx, fragment in enumerate(fragmentIsotopeList):
        if printZeros:
            print(massStr[fIdx])
        for sub in fragment:
            markerStyle = next(markerList)
            counts = np.array(Merged[fIdx]['counts' + sub])
            if CPmS == True:
                counts /= np.array(Merged[fIdx]['integTime'])
                ylabel = "Counts Per Millisecond"  
            else:
                ylabel = "Counts Per Scan"
                
            fig.suptitle(ylabel)
                
            if oneD == True:
                accessSubplot = 2 * (fIdx % 2)
            else: 
                accessSubplot = (fIdx // 2, 2* (fIdx % 2))
            
            cAx = ax[accessSubplot]
            cAx.scatter(Merged[fIdx].scanNumber,counts,label = sub,marker = markerStyle,alpha = 0.4)
            cAx.legend()
            cAx.set_title(massStr[fIdx])
            cAx.set_ylabel(ylabel)
            cAx.set_xlabel("Scan Number")
            
            if oneD == True:
                accessSubplot = 2 * (fIdx % 2) + 1
            else: 
                accessSubplot = (fIdx // 2, 2* (fIdx % 2) + 1)
                
            cAx = ax[accessSubplot]
            cAx.scatter(Merged[fIdx].scanNumber,counts,label = sub,marker = markerStyle,alpha = 0.4)
            cAx.set_ylim(0,200)
            
            if printZeros:
                print(sub + ' Zero Counts: ' + str(len(counts) - np.count_nonzero(counts)))

    plt.tight_layout()