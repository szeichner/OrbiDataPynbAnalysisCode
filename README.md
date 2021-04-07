# OrbiDataPynbAnalysisCode
Welcome to our "one code to rule them all," i.e., the Eiler Lab's collective coding effort to process Orbitrap data.

The main functionality of the code is in "DataAnalyzerWithPeakInteg.py", which has all the processing functions to convert RAW files into python-friendly data structures, and calculate isotope ratios from these data structures. These results are all output to the original folder/file location that the input RAW file comes from.

You can process data using "DataAnalyzerWorker.py", which has a sample of how you would call a folder with a combination of experiments, and what "toggles" you can turn on and off to process the data in different ways.

'dataQualityPlots.py' has functionality that allows you to plot individual acquisitions and assess the quality of data based on various parameters.

If you would like to run these functions from Python notebooks, there are a few within the "Test Folder" that test LC M+N experiment data, as well GC EI data with multiple input files. There are also files within the test folder that you can use to test code changes on M+N experiments, and GC data.

We hope this helps you process your orbitrap data and discover novel isotope anomalies! We welcome suggestions to improve the code, and implement new helpful functionality that would be generally useful to all users. Thanks!
